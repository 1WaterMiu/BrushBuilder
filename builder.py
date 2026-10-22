"""Fixed-silica brush assembly. Organic geometry is independently generated."""
import io, json, math, re, zipfile, hashlib
from pathlib import Path
from functools import lru_cache

ROOT = Path(__file__).resolve().parent
EXAMPLE = ROOT / 'exmaple'
TEMPLATES = {'CBAA': 'SUP_CBAA3N20', 'SBMA': 'SUP_SBMA20'}

def sections(path):
    result, section = {}, None
    for line in path.read_text().splitlines():
        line = line.split(';')[0].strip()
        if line.startswith('['):
            section = line.strip('[] ').lower()
            result.setdefault(section, [])
        elif line and not line.startswith('#') and section:
            result[section].append(line.split())
    return result

def gro(path):
    lines = path.read_text().splitlines()
    count = int(lines[1])
    atoms = [dict(res=int(l[:5]), residue=l[5:10].strip(), name=l[10:15].strip(), xyz=[float(l[a:a+8]) for a in (20,28,36)]) for l in lines[2:2+count]]
    return atoms, [float(x) for x in lines[2+count].split()]

@lru_cache(maxsize=2)
def template(key):
    if key not in TEMPLATES:
        raise ValueError('Unknown template')
    sub = sections(EXAMPLE/'SUB_O34-H4.itp')
    sup = sections(EXAMPLE/(TEMPLATES[key]+'.itp'))
    atoms, box = gro(EXAMPLE/(key+'.gro'))
    ns, np = len(sub['atoms']), len(sup['atoms'])
    blocks = [atoms[i*ns:(i+1)*ns] for i in range(140)]
    chains = [atoms[140*ns+i*np:140*ns+(i+1)*np] for i in range(20)]
    for block in blocks:
        assert [a['name'] for a in block] == [a[4] for a in sub['atoms']], 'SUB atom order mismatch'
    for chain in chains:
        assert [a['name'] for a in chain] == [a[4] for a in sup['atoms']], 'SUP atom order mismatch'
    return sub, sup, blocks, chains, box

def missing_includes():
    missing = set()
    for name in ('forcefield.itp','ffbonded.itp','ffnonbonded.itp'):
        for include in re.findall(r'^\s*#include\s+"([^"]+)"', (EXAMPLE/name).read_text(), re.M):
            if not (EXAMPLE/include).exists(): missing.add(include)
    return sorted(missing)

def periodic_delta(a,b,box):
    d = [a[i]-b[i] for i in range(3)]
    return [d[i]-round(d[i]/box[i])*box[i] for i in range(3)]

REMOVED_SITES = {2,7,12,17,4,9,14,19}

def select_sites(chains, box, seed=0):
    return [i for i in range(20) if i+1 not in REMOVED_SITES]

def restore_sub(chain, sub, box):
    # All existing silica coordinates survive. Restore missing O34/H4 along Si1-C1.
    source = {a['name']:a for a in chain[:52]}
    origin = source['Si1']['xyz']
    direction = [0,0,1]  # canonical graft normal, independent of source organic coordinates
    # Learn Si1-O34 and O34-H4 vectors from an actual bare SUB, not guessed charges.
    _,_,blocks,_,_ = template('CBAA')
    ref = {a['name']:a['xyz'] for a in blocks[0]}
    axis = periodic_delta(ref['O34'],ref['Si1'],box)
    distance = math.sqrt(sum(v*v for v in axis))
    u = [v/distance for v in axis]
    h = periodic_delta(ref['H4'],ref['O34'],box)
    # Minimal rotation mapping reference Si-O to new Si-O (Rodrigues).
    cross = [u[1]*direction[2]-u[2]*direction[1],u[2]*direction[0]-u[0]*direction[2],u[0]*direction[1]-u[1]*direction[0]]
    cosine = sum(u[i]*direction[i] for i in range(3))
    def cp(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
    if cosine < -0.999999:
        perpendicular = cp(u,[1,0,0] if abs(u[0])<0.9 else [0,1,0])
        norm=math.sqrt(sum(v*v for v in perpendicular)); perpendicular=[v/norm for v in perpendicular]
        rotated=[2*perpendicular[i]*sum(perpendicular[j]*h[j] for j in range(3))-h[i] for i in range(3)]
    else:
        first=cp(cross,h); second=cp(cross,first)
        rotated=[h[i]+first[i]+second[i]/(1+cosine) for i in range(3)]
    oxygen=[origin[i]+distance*direction[i] for i in range(3)]
    source['O34']={'xyz':oxygen}
    source['H4']={'xyz':[oxygen[i]+rotated[i] for i in range(3)]}
    return [dict(res=chain[0]['res'],residue='SUB',name=a[4],xyz=source[a[4]]['xyz'][:]) for a in sub['atoms']]

def build(key='CBAA', density='high', seed=0, repeats=20, methyl=True,head_group=None,linker=None,spacer=None):
    from polymer import make_chain, VERSION
    if key not in TEMPLATES:raise ValueError('Unknown template')
    from chemistry import normalize,reference,metadata,VERSION as HEAD_VERSION
    spec=normalize(key,head_group,linker,spacer)
    ref=reference(spec)
    key='CBAA' if spec[1]=='amide' else 'SBMA'
    chemical=metadata(spec,methyl)
    if density not in ('high','low'): raise ValueError('Density must be high or low')
    if not isinstance(repeats,int) or not 1<=repeats<=60: raise ValueError('Repeat units must be an integer from 1 to 60')
    if not isinstance(methyl,bool): raise ValueError('Methyl must be boolean')
    if not isinstance(seed,int) or not 0<=seed<=19: raise ValueError('Site seed must be 0–19')
    sub,sup,blocks,chains,box=template(key)
    keep=list(range(20)) if density=='high' else select_sites(chains,box,seed)
    converted=[restore_sub(c,sub,box) for i,c in enumerate(chains) if i not in keep]
    generated=[make_chain(key,repeats,methyl,chains[i],spec) for i in keep]
    molecules=[(b,sub) for b in blocks+converted]+generated
    sup=generated[0][1]
    atoms=[]; frozen=[]; charge=0; bonds=[]
    for block,top in molecules:
        offset=len(atoms)
        for a,t in zip(block,top['atoms']):
            idx=len(atoms)+1
            is_silica=t[3] in ('SUB','SUP')
            if is_silica: frozen.append(idx)
            q=None if t[6]=='UNASSIGNED' else float(t[6])
            atoms.append(dict(a,charge=q,type=t[1],frozen=is_silica))
            if q is not None:charge+=q
        bonds.extend([[offset+int(b[0])-1,offset+int(b[1])-1] for b in top['bonds']])
    area=box[0]*box[1]
    box=box[:]; box[2]=max(box[2],max(a['xyz'][2] for a in atoms)+4)
    from cgenff_library import status as parameter_status
    chemical['parameterization']=parameter_status(spec,methyl,repeats)
    complete=all(a['charge'] is not None for a in atoms)
    chemical['parameterization']['chargesComplete']=complete
    return dict(template=key,chemistry=chemical,parameterized=chemical['parameterization']['mdReady'],density=density,repeats=repeats,methyl=methyl,geometryVersion=HEAD_VERSION,sidechainOrientation='alternating reflection; identical sequence for every chain',box=box,area=area,chainCount=len(keep),graftDensity=len(keep)/area,subCount=140+len(converted),atoms=atoms,bonds=bonds,frozen=frozen,charge=charge if complete else None,keptSites=[i+1 for i in keep],removedSites=sorted(REMOVED_SITES) if density=='low' else [],topology=sup,topologyCounts={k:len(v) for k,v in sup.items()},missingIncludes=missing_includes(),status='generated',minimization=None)

def catalog():
    result=[]
    for key in TEMPLATES:
        sub,sup,_,_,box=template(key)
        groups={}
        for atom in sup['atoms']: groups.setdefault(atom[2],[]).append(atom)
        first=list(groups.values())[2]
        ids={int(a[0]) for a in first}
        result.append(dict(id=key,name='CBMAA · spacer 3 · amide (source CBAA3)' if key=='CBAA' else 'SBMA · ester',linker='amide' if key=='CBAA' else 'ester',repeats=len(groups)-2,methyl=True,atomsPerChain=len(sup['atoms']),charge=sum(float(a[6]) for a in sup['atoms']),area=box[0]*box[1],anchorAtoms=sup['atoms'][:74],repeatAtoms=first,repeatBonds=[b for b in sup['bonds'] if int(b[0]) in ids and int(b[1]) in ids],source=TEMPLATES[key]+'.itp'))
    from chemistry import FAMILIES,template_comparison
    return dict(templates=result,headGroups=FAMILIES,linkers=['ester','amide'],spacers=[1,2,3,4],templateComparisons=[template_comparison(k) for k in TEMPLATES],missingIncludes=missing_includes())

def archive(model):
    if not model.get('parameterized',True):raise ValueError('New chemistry has no assigned partial charges or force-field types; MD export is unavailable')
    from polymer import itp_text
    from simulation import EM_MDP, NVT_MDP, RUN_SCRIPT, INDEX_SCRIPT
    memory=io.BytesIO()
    with zipfile.ZipFile(memory,'w',zipfile.ZIP_DEFLATED) as z:
        lines=['BrushBuilder dry surface; silica fixed',str(len(model['atoms']))]
        for i,a in enumerate(model['atoms'],1):
            x,y,v=a['xyz']; lines.append(f"{a['res']%100000:5d}{a['residue']:<5.5}{a['name']:>5.5}{i%100000:5d}{x:8.3f}{y:8.3f}{v:8.3f}")
        lines.append(' '.join(f'{v:.5f}' for v in model['box']))
        z.writestr('surface.gro','\n'.join(lines)+'\n')
        ids=model['frozen']; frozen_set=set(ids); mobile=[i for i in range(1,len(model['atoms'])+1) if i not in frozen_set]
        index=''
        for name,values in [('System',list(range(1,len(model['atoms'])+1))),('Silica',ids),('Mobile',mobile)]:
            index+=f'[ {name} ]\n'+'\n'.join(' '.join(map(str,values[i:i+15])) for i in range(0,len(values),15))+'\n'
        z.writestr('index.ndx',index)
        z.writestr('freeze.mdp','; Merge into EVERY simulation stage; not a complete MDP\nfreezegrps = Silica\nfreezedim = Y Y Y\n; Keep lateral box dimensions fixed. Disable COM removal on Silica.\n; Rebuild Mobile/index groups after solvation.\n')
        topname='BRUSH'
        z.writestr('BRUSH.itp',itp_text(model['topology']))
        z.writestr('topol.top',f'#include "forcefield.itp"\n#include "SUB_O34-H4.itp"\n#include "{topname}.itp"\n\n[ system ]\nBrushBuilder dry surface\n\n[ molecules ]\nSUB {model["subCount"]}\n{topname} {model["chainCount"]}\n')
        files=['forcefield.itp','ffbonded.itp','ffnonbonded.itp','SUB_O34-H4.itp']
        for optional in ['cmap.itp','nbfix.itp','ffmissingdihedrals.itp']:
            if (EXAMPLE/optional).exists(): files.append(optional)
        for file in files: z.write(EXAMPLE/file,file)
        report={k:v for k,v in model.items() if k not in ('atoms','bonds','frozen','topology')}
        report['frozenAtoms']=ids
        report['atomCount']=len(model['atoms']); report['frozenCount']=len(ids)
        report['sourceSHA256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [EXAMPLE/(model['template']+'.gro'),*[EXAMPLE/f for f in files]]}
        report['generatorSHA256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['builder.py','polymer.py','simulation.py','chemistry.py']}
        z.writestr('manifest.json',json.dumps(report,indent=2))
        for p in sorted((ROOT/'geometry_templates').glob('*.json')): z.write(p,'geometry_templates/'+p.name)
        for p in sorted((ROOT/'parameter_audit').glob('*.json')): z.write(p,'parameter_audit/'+p.name)
        for f in ['tip3p.itp','ions.itp']: z.write(EXAMPLE/f,f)
        z.writestr('em.mdp',EM_MDP); z.writestr('nvt.mdp',NVT_MDP)
        z.writestr('run_md.sh',RUN_SCRIPT); z.writestr('update_index.py',INDEX_SCRIPT)
        z.writestr('README.txt','Generated organic coordinates are independent of source GRO coordinates.\nSilica is frozen in all XYZ directions during every stage. Low removes sites 2,7,12,17,4,9,14,19 (12 chains remain).\nSee manifest.json for minimization/convergence and frozen-coordinate verification. Only a completed job is MD-ready.\nRun bash run_md.sh with GROMACS and Python 3 for TIP3P solvation, solvated minimization and a 100 ps NVT example. No protein or added salt; customize duration/conditions for your study. Box area is fixed, no pressure coupling.\nAlpha-H substitution preserves charge using an explicit +0.09 H / -0.09 alpha-C increment and valence-based aliphatic types. This is an assigned force-field patch, not newly QM-validated chemistry.\n')
    return memory.getvalue()
