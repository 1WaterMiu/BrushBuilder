"""Audit actual web results against explicit-H input graphs and atom maps."""
from pathlib import Path
import json,zipfile,re,hashlib
from collections import defaultdict
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1]/'cgenff'

def parse(text):
    result=dict(atoms={},bonds=[],impropers=[],parameters=defaultdict(list))
    mode='topology';section=None
    for line in text.splitlines():
        raw=line.split('!')[0].split()
        if not raw:continue
        if line.lower().startswith('read param'):mode='parameters'
        if mode=='topology':
            if raw[0]=='RESI':
                result['charge']=float(raw[2])
                result['parameterPenalty']=float(re.search(r'param penalty=\s*([\d.]+)',line)[1])
                result['chargePenalty']=float(re.search(r'charge penalty=\s*([\d.]+)',line)[1])
            if raw[0]=='ATOM':
                result['atoms'][raw[1]]=dict(type=raw[2],charge=float(raw[3]),penalty=float(line.split('!')[1].split()[0]))
            if raw[0]=='BOND':result['bonds'].extend([raw[i:i+2] for i in range(1,len(raw),2)])
            if raw[0]=='IMPR':result['impropers'].extend([raw[i:i+4] for i in range(1,len(raw),4)])
        elif raw[0] in ('BONDS','ANGLES','DIHEDRALS','IMPROPERS','NONBONDED','NBFIX','END'):section=raw[0]
        elif section and section!='END' and not line.startswith('*'):result['parameters'][section].append(raw)
    return result

def coverage(data):
    adj={a:set() for a in data['atoms']}
    for a,b in data['bonds']:adj[a].add(b);adj[b].add(a)
    angles=[(a,b,c) for b,neighbors in adj.items() for a,c in combinations(sorted(neighbors),2)]
    torsions=set()
    for b,c in data['bonds']:
        for a in adj[b]-{c}:
            for d in adj[c]-{b,a}:torsions.add(min((a,b,c,d),(d,c,b,a)))
    terms={'BONDS':data['bonds'],'ANGLES':angles,'DIHEDRALS':sorted(torsions),'IMPROPERS':data['impropers']}
    missing={};counts={}
    for section,rows in terms.items():
        absent=set();counts[section]=len(rows)
        for row in rows:
            types=tuple(data['atoms'][a]['type'] for a in row)
            def matches(p):return any(all(x=='X' or x==y for x,y in zip(p,t)) for t in (types,types[::-1]))
            if not any(matches(p[:len(types)]) for p in data['parameters'].get(section,[])):absent.add(types)
        missing[section]=sorted(absent)
    return dict(counts=counts,missingTypes=missing,complete=not any(missing.values()))

def main():
    reports=[]
    for row in json.loads((ROOT/'manifest.json').read_text()):
        stem=Path(row['file']).stem;zpath=ROOT/'results'/f'{stem}.zip'
        if not zpath.exists():continue
        with zipfile.ZipFile(zpath) as z:
            text=z.read(stem+'.str').decode();err=z.read(stem+'.err').decode()
            parsed=parse(text)
            output=ROOT/'parsed';output.mkdir(exist_ok=True)
            (output/(stem+'.str')).write_text(text,encoding='utf-8')
            (output/(stem+'.err')).write_text(err,encoding='utf-8')
        info=json.loads((ROOT/'inputs'/(stem+'.json')).read_text())
        names={a['cgenffName'] for a in info['atomMap']}
        mol2=(ROOT/'inputs'/row['file']).read_text();lines=mol2.split('@<TRIPOS>BOND')[1].split('@<TRIPOS>')[0].strip().splitlines()
        mapping={i+1:a['cgenffName'] for i,a in enumerate(info['atomMap'])}
        expected={tuple(sorted((mapping[int(x.split()[1])],mapping[int(x.split()[2])]))) for x in lines}
        actual={tuple(sorted(e)) for e in parsed['bonds']}
        atom_match=names==set(parsed['atoms']);bond_match=expected==actual
        if not atom_match or not bond_match:raise ValueError(f'Input/output graph mismatch: {stem}')
        units=defaultdict(float)
        for a in info['atomMap']:units[a['unit']]+=parsed['atoms'][a['cgenffName']]['charge']
        report=dict(file=stem,atomCount=len(names),bondCount=len(actual),improperCount=len(parsed['impropers']),bondedCoverage=coverage(parsed),
          atomMapExact=atom_match,bondGraphExact=bond_match,charge=sum(a['charge'] for a in parsed['atoms'].values()),
          unitCharges={k:round(v,6) for k,v in units.items()},chargePenalty=parsed['chargePenalty'],parameterPenalty=parsed['parameterPenalty'],
          inputSHA256=info['sha256'],resultSHA256=hashlib.sha256(zpath.read_bytes()).hexdigest(),version='5.0')
        parsed['provenance']=report
        (output/(stem+'.json')).write_text(json.dumps(parsed,indent=2),encoding='utf-8')
        reports.append(report)
    (ROOT/'audit.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
    print('Audited',len(reports),'results; exact atom and bond graphs.')
    print('Bonded coverage complete:',sum(r['bondedCoverage']['complete'] for r in reports),'/',len(reports))
    print('High penalty:',[(x['file'],x['chargePenalty'],x['parameterPenalty']) for x in reports if max(x['chargePenalty'],x['parameterPenalty'])>50])

if __name__=='__main__':main()
