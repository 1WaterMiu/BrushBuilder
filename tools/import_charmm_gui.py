"""Audit CHARMM-GUI archives and import CGenFF data without template charges.

Atom order is accepted only after element, formal-charge and bond-order checks
against the submitted SDF and a separate RTF connectivity check. Fail closed if
the service reorders atoms; do not guess a correspondence from atom names.
"""
from pathlib import Path
import sys, json, tarfile, hashlib, zipfile
from collections import defaultdict
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.audit_cgenff import parse, coverage
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]/'cgenff'

def normalized(m):
    if m is None: raise ValueError('Unreadable structure')
    for c in m.GetAtoms():
        oo=[a for a in c.GetNeighbors() if a.GetSymbol()=='O' and a.GetDegree()==1]
        if c.GetSymbol()=='C' and len(oo)==2 and sum(a.GetFormalCharge() for a in oo)==-1:
            for o in oo:
                b=m.GetBondBetweenAtoms(c.GetIdx(),o.GetIdx())
                b.SetBondType(Chem.BondType.SINGLE if o.GetFormalCharge()==-1 else Chem.BondType.DOUBLE)
                b.SetIsAromatic(False);o.SetIsAromatic(False)
            c.SetIsAromatic(False)
    Chem.SanitizeMol(m)
    return ([(a.GetSymbol(),a.GetFormalCharge()) for a in m.GetAtoms()],
            sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),b.GetBondTypeAsDouble()) for b in m.GetBonds()))

def import_archive(path):
    stem,job=path.stem.rsplit('_',1)
    info=json.loads((ROOT/'inputs'/f'{stem}.json').read_text())
    if path.suffix=='.zip':
        with zipfile.ZipFile(path) as archive:
            files={n:archive.read(n).decode() for n in archive.namelist() if n.endswith(('.rtf','.prm','.mol','.out'))}
    else:
        with tarfile.open(path) as archive:
            files={m.name:archive.extractfile(m).read().decode() for m in archive.getmembers()
                   if m.isfile() and (m.name.endswith(('.rtf','.prm','.mol','.out')))}
    def get(suffix):
        found=[v for k,v in files.items() if k.endswith(suffix)]
        if len(found)!=1:raise ValueError(f'Ambiguous archive member: {suffix}')
        return found[0]
    rtfs=[v for k,v in files.items() if k.endswith('.rtf') and '/toppar/' not in k and not k.endswith('_g.rtf')]
    if len(rtfs)!=1:raise ValueError('Expected one generated residue')
    rtf=rtfs[0];data=parse(rtf)
    if not all(s in rtf for s in ('program version 4.0','parameter files version 5.0')):
        raise ValueError('Unexpected CGenFF version')
    old=normalized(Chem.MolFromMolFile(str(ROOT/'inputs'/f'{stem}.sdf'),removeHs=False))
    new=normalized(Chem.MolFromMolBlock(get('/drawing_3D.mol'),removeHs=False))
    if old!=new:raise ValueError(f'Ordered chemical graph differs: {stem}')
    names=list(data['atoms']);target=[a['cgenffName'] for a in info['atomMap']]
    if len(names)!=len(target):raise ValueError('Atom count differs')
    if [n[0] for n in names]!=[a[0] for a in old[0]]:raise ValueError('RTF element order differs')
    indices={n:i for i,n in enumerate(names)}
    expected={(a,b) for a,b,_ in old[1]}
    actual={tuple(sorted((indices[a],indices[b]))) for a,b in data['bonds']}
    if actual!=expected:raise ValueError('RTF ordered connectivity differs')
    mapping=dict(zip(names,target))
    data['atoms']={mapping[n]:a for n,a in data['atoms'].items()}
    for section in ('bonds','impropers'):data[section]=[[mapping[a] for a in t] for t in data[section]]
    base=parse('read param\n'+get('/toppar/par_all36_cgenff.prm'))['parameters']
    custom_text=[v for k,v in files.items() if k.endswith('.prm') and '/toppar/' not in k]
    if len(custom_text)!=1:raise ValueError('Expected one custom parameter file')
    custom=parse('read param\n'+custom_text[0])['parameters']
    for section,n in [('BONDS',2),('ANGLES',3),('DIHEDRALS',4),('IMPROPERS',4)]:
        rows=custom.get(section,[])
        keys={min(tuple(r[:n]),tuple(r[:n][::-1])) for r in rows}
        data['parameters'][section]=[r for r in base.get(section,[]) if min(tuple(r[:n]),tuple(r[:n][::-1])) not in keys]+rows
        used={a['type'] for a in data['atoms'].values()}|{'X'}
        data['parameters'][section]=[r for r in data['parameters'][section] if all(t in used for t in r[:n])]
    cov=coverage(data)
    if not cov['complete']:raise ValueError(f'Missing bonded terms: {cov}')
    charge=sum(a['charge'] for a in data['atoms'].values())
    if abs(charge-info['formalCharge'])>1e-5:raise ValueError('Charge sum mismatch')
    units=defaultdict(float)
    for a in info['atomMap']:units[a['unit']]+=data['atoms'][a['cgenffName']]['charge']
    report=dict(file=stem,source='CHARMM-GUI / CGenFF',jobId=job,atomCount=len(names),bondCount=len(actual),
        improperCount=len(data['impropers']),bondedCoverage=cov,atomMapExact=True,bondGraphExact=True,
        orderedChemicalGraphExact=True,charge=charge,unitCharges={k:round(v,6) for k,v in units.items()},
        chargePenalty=data['chargePenalty'],parameterPenalty=data['parameterPenalty'],
        inputSHA256=info['sha256'],resultSHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
        resultArchive=str(path.relative_to(ROOT)),version='5.0',programVersion='4.0',guiToInputAtomMap=mapping)
    # Keep the full download locally; commit a compact, lossless archive of
    # scientific source files rather than duplicate every unrelated force field.
    compact=path.with_suffix('.zip')
    selected={k:v for k,v in files.items() if '/toppar/' not in k or k.endswith('/toppar/par_all36_cgenff.prm')}
    with zipfile.ZipFile(compact,'w',zipfile.ZIP_DEFLATED) as z:
        for name,text in sorted(selected.items()):
            item=zipfile.ZipInfo(name,date_time=(2026,1,1,0,0,0));item.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(item,text)
    report['sourceArchive']=str(compact.relative_to(ROOT))
    report['sourceArchiveSHA256']=hashlib.sha256(compact.read_bytes()).hexdigest()
    data['provenance']=report
    out=ROOT/'parsed';out.mkdir(exist_ok=True)
    destination=out/f'{stem}.json'
    if destination.exists() and json.loads(destination.read_text()).get('provenance',{}).get('source')!='CHARMM-GUI / CGenFF':
        raise ValueError('Refusing to overwrite original CGenFF result')
    temporary=destination.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(data,indent=2),encoding='utf-8')
    temporary.replace(destination)
    # Preserve original output and atom mapping; never label this as an original STR.
    (out/f'{stem}.rtf').write_text(rtf,encoding='utf-8')
    (out/f'{stem}.prm').write_text(custom_text[0],encoding='utf-8')
    return report

if __name__=='__main__':
    candidates={p.stem:p for ext in ('*.zip','*.tgz') for p in sorted((ROOT/'charmm-gui-results').glob(ext))}
    reports=[import_archive(p) for _,p in sorted(candidates.items())]
    (ROOT/'charmm-gui-audit.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
    print(json.dumps([{'file':r['file'],'atoms':r['atomCount'],'coverage':r['bondedCoverage']['complete'],'penalty':r['parameterPenalty']} for r in reports],indent=2))
