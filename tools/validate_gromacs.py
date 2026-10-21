"""Integration checks against the actual GROMACS executable."""
import io,json,sys,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from builder import ROOT,build,archive,gro
from simulation import engine,run_command,NVT_MDP,INDEX_SCRIPT

command=engine(); results=[]
if '--solvated' in sys.argv:
    import shutil
    directory=ROOT/'output/solvated-validation'; directory.mkdir(exist_ok=True)
    source=ROOT/'output/final-CBAA20'
    if not (source/'surface-md.zip').exists():
        import uuid
        from simulation import JOBS,work
        ident=uuid.uuid4().hex
        source=ROOT/'output'/('validation-'+ident)
        JOBS[ident]=dict(id=ident,status='queued')
        work(ident,dict(key='CBAA',density='high',repeats=20,methyl=True),source)
        assert JOBS[ident]['status']=='completed',JOBS[ident]
    with zipfile.ZipFile(source/'surface-md.zip') as z:z.extractall(directory)
    text=(directory/'topol.top').read_text().replace('[ system ]','#include "tip3p.itp"\n#include "ions.itp"\n\n[ system ]')
    (directory/'solvated.top').write_text(text)
    run_command(command+['solvate','-cp','surface.gro','-cs','spc216.gro','-o','solvated.gro','-p','solvated.top'],directory,'solvate.log')
    run_command([sys.executable,'update_index.py','solvated.gro','surface.gro'],directory,'index.log')
    run_command(command+['grompp','-f','em.mdp','-c','solvated.gro','-p','solvated.top','-n','index.ndx','-o','wet-em.tpr'],directory,'wet-grompp.log')
    run_command(command+['mdrun','-deffnm','wet-em','-ntmpi','1','-ntomp','4'],directory,'wet-mdrun.log')
    assert 'converged to Fmax <' in (directory/'wet-em.log').read_text()
    run_command([sys.executable,'update_index.py','wet-em.gro','surface.gro'],directory,'wet-fixed.log')
    (directory/'smoke.mdp').write_text(NVT_MDP.replace('nsteps = 100000','nsteps = 10').replace('nstxout-compressed = 1000','nstxout-compressed = 1'))
    run_command(command+['grompp','-f','smoke.mdp','-c','wet-em.gro','-p','solvated.top','-n','index.ndx','-o','smoke.tpr'],directory,'smoke-grompp.log')
    run_command(command+['mdrun','-deffnm','smoke','-ntmpi','1','-ntomp','4'],directory,'smoke-mdrun.log')
    run_command([sys.executable,'update_index.py','smoke.gro','surface.gro'],directory,'smoke-fixed.log')
    print('Solvated EM converged; 10-step NVT smoke test passed; silica coordinates unchanged.')
else:
    for key in ['CBAA','SBMA']:
        for methyl in [False,True]:
            for n in [1,60]:
                directory=ROOT/'output/preflight'/f'{key}-{n}-{methyl}';directory.mkdir(parents=True,exist_ok=True)
                model=build(key,'low',repeats=n,methyl=methyl)
                with zipfile.ZipFile(io.BytesIO(archive(model))) as z:z.extractall(directory)
                try:
                    run_command(command+['grompp','-f','em.mdp','-c','surface.gro','-p','topol.top','-n','index.ndx','-o','em.tpr'],directory,'grompp.log')
                    results.append(dict(template=key,repeats=n,methyl=methyl,preprocessing='passed'))
                except RuntimeError as e:
                    import re
                    missing=[]
                    for term in re.findall(r"'([0-9 ]+) 9'",str(e)):
                        missing.append([model['topology']['atoms'][int(i)-1][1] for i in term.split()])
                    results.append(dict(template=key,repeats=n,methyl=methyl,preprocessing='failed',missing=missing))
                print(results[-1],flush=True)
    (ROOT/'parameter_audit/preprocessing_checks.json').write_text(json.dumps(results,indent=2))
