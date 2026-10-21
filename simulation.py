"""Actual GROMACS runs, no mocked optimization. All silica atoms are frozen."""
EM_MDP = '''integrator = steep
nsteps = 20000
emtol = 500
emstep = 0.005
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw-switch = 1.0
rvdw = 1.2
pbc = xyz
constraints = none
freezegrps = Silica
freezedim = Y Y Y
comm-mode = None
pcoupl = no
tcoupl = no
nstenergy = 100
nstlog = 100
'''
NVT_MDP = '''integrator = md
dt = 0.001
nsteps = 100000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.3
; Explicit 0.1 nm buffer: frozen, unthermostatted silica has zero DOF.
verlet-buffer-tolerance = -1
coulombtype = PME
rcoulomb = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw-switch = 1.0
rvdw = 1.2
pbc = xyz
; Do not constrain silica bonds: constraint initialization can move frozen H.
constraints = none
constraint-algorithm = lincs
freezegrps = Silica
freezedim = Y Y Y
comm-mode = None
pcoupl = no
tcoupl = v-rescale
tc-grps = Mobile Silica
tau-t = 1.0 -1
ref-t = 300 300
gen-vel = yes
gen-temp = 300
gen-seed = 1947
nstenergy = 1000
nstlog = 1000
nstxout-compressed = 1000
'''
INDEX_SCRIPT = '''import json,sys
from pathlib import Path
manifest=json.loads(Path('manifest.json').read_text())
frozen=manifest['frozenAtoms']; count=int(Path(sys.argv[1]).read_text().splitlines()[1])
fixed=set(frozen)
with open('index.ndx','w') as f:
 for name,ids in [('System',range(1,count+1)),('Silica',frozen),('Mobile',[i for i in range(1,count+1) if i not in fixed])]:
  ids=list(ids); f.write('[ '+name+' ]\\n')
  for start in range(0,len(ids),15): f.write(' '.join(map(str,ids[start:start+15]))+'\\n')
if len(sys.argv)>2:
 before=Path(sys.argv[2]).read_text().splitlines()[2:]
 after=Path(sys.argv[1]).read_text().splitlines()[2:]
 if any(before[i-1][20:44]!=after[i-1][20:44] for i in frozen):
  raise SystemExit('FATAL: silica moved')
 print('Verified: every silica coordinate unchanged')
'''
RUN_SCRIPT = '''#!/usr/bin/env bash
set -euo pipefail
if [ -z "${GMX:-}" ]; then
  if command -v gmx_d >/dev/null; then GMX=gmx_d; else GMX=gmx; fi
fi
python3 -c "import json; assert json.load(open('manifest.json'))['status']=='completed', 'Not a validated minimized export'"
cp surface.gro reference.gro
# Keep the original dry topology reproducible on repeated runs.
cp topol.top solvated.top
python3 - <<'PY'
from pathlib import Path
p=Path('solvated.top'); s=p.read_text(); s=s.replace('[ system ]','#include "tip3p.itp"\\n#include "ions.itp"\\n\\n[ system ]'); p.write_text(s)
PY
"$GMX" solvate -cp surface.gro -cs spc216.gro -o solvated.gro -p solvated.top
python3 update_index.py solvated.gro reference.gro
"$GMX" grompp -f em.mdp -c solvated.gro -p solvated.top -n index.ndx -o solvated-em.tpr
"$GMX" mdrun -deffnm solvated-em -ntmpi 1 -ntomp 4
python3 update_index.py solvated-em.gro reference.gro
python3 - <<'PY'
from pathlib import Path
assert 'converged to Fmax <' in Path('solvated-em.log').read_text(), 'Solvated EM not converged'
PY
"$GMX" grompp -f nvt.mdp -c solvated-em.gro -p solvated.top -n index.ndx -o nvt.tpr
"$GMX" mdrun -deffnm nvt -ntmpi 1 -ntomp 4
python3 update_index.py nvt.gro reference.gro
'''

import io, json, re, shutil, subprocess, threading, uuid, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

POOL=ThreadPoolExecutor(max_workers=1)
JOBS={}
LOCK=threading.Lock()

def engine():
    native=shutil.which('gmx_d') or shutil.which('gmx')
    if native: return [native]
    if shutil.which('wsl'):
        check=subprocess.run(['wsl','-d','Ubuntu','--','bash','-lc','if [ -x /home/watermiu/gromacs-2025.3-double/bin/gmx_d ]; then printf /home/watermiu/gromacs-2025.3-double/bin/gmx_d; else command -v gmx_d || command -v gmx; fi'],capture_output=True,text=True,timeout=20)
        if check.returncode==0 and check.stdout.strip(): return ['wsl','-d','Ubuntu','--',check.stdout.strip()]
    raise RuntimeError('GROMACS unavailable; no minimized result will be claimed')

def run_command(command,directory,log):
    if command[0]=='wsl':
        path='/mnt/'+str(directory)[0].lower()+str(directory)[2:].replace('\\','/')
        command=command[:3]+['--cd',path]+command[3:]
    with open(directory/log,'w') as stream, open(directory/(log+'.stderr'),'w') as errors:
        result=subprocess.run(command,cwd=directory,stdout=stream,stderr=errors,timeout=1800)
    if result.returncode:
        text=(directory/log).read_text(errors='replace')+(directory/(log+'.stderr')).read_text(errors='replace')
        raise RuntimeError(log+': '+text[-3500:])

def submit(parameters):
    from builder import ROOT
    ident=uuid.uuid4().hex
    with LOCK: JOBS[ident]=dict(id=ident,status='queued',message='等待 GROMACS',parameters=parameters,directory=str(ROOT/'output'/'jobs'/ident))
    POOL.submit(work,ident,parameters,ROOT/'output'/'jobs'/ident)
    return ident

def work(ident,parameters,directory):
    from builder import build,archive,gro
    def update(**fields):
        with LOCK:
            JOBS[ident].update(fields)
            if directory.exists(): (directory/'job.json').write_text(json.dumps(JOBS[ident]))
    try:
        update(status='building',message='独立生成坐标与完整拓扑')
        model=build(**parameters); directory.mkdir(parents=True)
        with zipfile.ZipFile(io.BytesIO(archive(model))) as z: z.extractall(directory)
        initial_manifest=json.loads((directory/'manifest.json').read_text())
        shutil.copy2(directory/'surface.gro',directory/'initial.gro')
        command=engine()
        update(status='preprocessing',message='GROMACS 检查 bonded / nonbonded 参数')
        run_command(command+['grompp','-f','em.mdp','-c','surface.gro','-p','topol.top','-n','index.ndx','-o','em.tpr','-pp','processed.top'],directory,'grompp.log')
        update(status='minimizing',message='正在最小化聚合物；silica XYZ 冻结')
        run_command(command+['mdrun','-deffnm','em','-ntmpi','1','-ntomp','4'],directory,'mdrun.log')
        before,_=gro(directory/'initial.gro'); after,_=gro(directory/'em.gro')
        if len(before)!=len(after): raise RuntimeError('Atom count changed during minimization')
        displacement=max(max(abs(a-b) for a,b in zip(before[i-1]['xyz'],after[i-1]['xyz'])) for i in model['frozen'])
        if displacement!=0: raise RuntimeError(f'Silica moved: {displacement} nm; export blocked')
        log=(directory/'em.log').read_text(errors='replace')
        converged='converged to Fmax <' in log
        if not converged: raise RuntimeError('EM did not reach Fmax < 500 kJ/mol/nm. Inspect em.log; no MD-ready export.')
        fmax=re.findall(r'Maximum force\s*=\s*([^\s]+)',log)
        model['status']='completed'; model['minimization']=dict(engine=command[-1],converged=True,maxSilicaDisplacementNm=displacement,fmax=float(fmax[-1]) if fmax else None,coordinateCheckPrecisionNm=.001)
        for a,b in zip(model['atoms'],after): a['xyz']=b['xyz']
        # Verify full-precision frozen TPR/TRR coordinates via gmx dump as well.
        run_command(command+['dump','-s','em.tpr'],directory,'tpr-dump.txt')
        run_command(command+['dump','-f','em.trr'],directory,'trr-dump.txt')
        pattern=r'(?m)^\s*x\[\s*(\d+)\]\s*=\s*\{\s*([^}]+)\}'
        initial={int(i):tuple(float(x) for x in v.split(',')) for i,v in re.findall(pattern,(directory/'tpr-dump.txt').read_text())}
        final={int(i):tuple(float(x) for x in v.split(',')) for i,v in re.findall(pattern,(directory/'trr-dump.txt').read_text())}
        if any(i-1 not in initial or i-1 not in final for i in model['frozen']): raise RuntimeError('Full-precision frozen-coordinate verification unavailable')
        exact=max(max(abs(a-b) for a,b in zip(initial[i-1],final[i-1])) for i in model['frozen'])
        if exact!=0: raise RuntimeError(f'Frozen silica moved in full-precision output: {exact} nm')
        model['minimization']['fullPrecisionMaxSilicaDisplacementNm']=exact
        # Keep the exact parameter snapshot that grompp used, even if source
        # force-field files are edited while this job is running.
        with zipfile.ZipFile(io.BytesIO(archive(model))) as z:
            (directory/'surface.gro').write_bytes(z.read('surface.gro'))
        initial_manifest['status']='completed'
        initial_manifest['minimization']=model['minimization']
        (directory/'manifest.json').write_text(json.dumps(initial_manifest,indent=2))
        (directory/'model.json').write_text(json.dumps(model))
        files=['surface.gro','initial.gro','topol.top','BRUSH.itp','SUB_O34-H4.itp','forcefield.itp','ffbonded.itp','ffnonbonded.itp','cmap.itp','nbfix.itp','ffmissingdihedrals.itp','tip3p.itp','ions.itp','index.ndx','freeze.mdp','em.mdp','nvt.mdp','run_md.sh','update_index.py','manifest.json','README.txt','em.log','grompp.log','grompp.log.stderr','mdrun.log','mdrun.log.stderr','em.tpr','processed.top']
        with zipfile.ZipFile(directory/'surface-md.zip','w',zipfile.ZIP_DEFLATED) as z:
            for name in files: z.write(directory/name,name)
            for folder in ['geometry_templates','parameter_audit']:
                for file in sorted((directory/folder).glob('*.json')): z.write(file,file.relative_to(directory))
        update(status='completed',message='最小化收敛；silica 位移为 0；完整模拟包可下载',result=model['minimization'],directory=str(directory))
    except Exception as e: update(status='failed',message=str(e),directory=str(directory))

def job(ident):
    from builder import ROOT
    if not re.fullmatch('[a-f0-9]{32}',ident): raise ValueError('Invalid job ID')
    with LOCK:
        if ident not in JOBS:
            directory=ROOT/'output'/'jobs'/ident
            if (directory/'job.json').exists():
                record=json.loads((directory/'job.json').read_text())
                if record['status'] not in ('completed','failed'):
                    record.update(status='failed',message='服务重启中断了任务；请重新提交')
                JOBS[ident]=record
            elif (directory/'surface-md.zip').exists() and (directory/'model.json').exists():
                JOBS[ident]=dict(id=ident,status='completed',message='已恢复完成的优化结果',directory=str(directory))
            else: raise ValueError('Unknown job')
        return JOBS[ident].copy()
