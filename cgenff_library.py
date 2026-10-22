"""Read only audited CGenFF results. Never fall back to template organic charges."""
import json
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parent/'cgenff'
INTERFACE={'C1','H4','H5'}

def key(spec,methyl):return f'{spec[0]}_{spec[1]}_{spec[2] or 0}_m{int(methyl)}'

@lru_cache(maxsize=160)
def load(stem):
    path=ROOT/'parsed'/(stem+'.json')
    if not path.exists():return None
    data=json.loads(path.read_text(encoding='utf-8'))
    info=json.loads((ROOT/'inputs'/(stem+'.json')).read_text(encoding='utf-8'))
    data['mapping']={(a['unit'],a['name']):a['cgenffName'] for a in info['atomMap']}
    data['reverse']={v:k for k,v in data['mapping'].items()}
    return data

def status(spec,methyl,n=20):
    base=key(spec,methyl)
    mono=load(base+'_monomer');free=load(base+'_pentamer');anchor=load(base+'_anchor5')
    complete=bool(anchor) and n>=5
    return dict(source='CGenFF / CHARMM-GUI / force field 5.0',monomerAvailable=bool(mono),
        polymerContextAvailable=bool(free),anchorContextAvailable=bool(anchor),
        chargesComplete=complete,mdReady=False,
        reason='GROMACS parameter coverage and minimization validation are pending' if complete else
          ('Short chains require dedicated anchored CGenFF results' if anchor else 'CGenFF results for this combination or its anchor are incomplete'),
        chargePenalty=(anchor or free or mono or {}).get('chargePenalty'),
        parameterPenalty=(anchor or free or mono or {}).get('parameterPenalty'),
        interfacePolicy='Only silica and anchor C1/H4/H5 retain source charges/types; all other organic atoms require CGenFF.',
        chargeNormalization='none')

def apply(topology,groups,spec,methyl):
    n=len(groups)-1;base=key(spec,methyl)
    anchor=load(base+'_anchor5');free=load(base+'_pentamer')
    result=anchor or free
    assigned={};provenance=[]
    for r,group in enumerate(groups):
        for name,i in group.items():
            row=topology['atoms'][i-1]
            if r==0 and name in INTERFACE:continue
            row[1]=row[6]='UNASSIGNED'
            if n<5 or result is None:continue
            if r==0:unit=-1 if anchor else None
            elif r in (1,2):unit=r if anchor else None
            elif r==n:unit=5
            elif r==n-1:unit=4
            else:unit=3
            lookup=(6,'END_H') if name=='H14' else (unit,name)
            atom_name=result['mapping'].get(lookup)
            if atom_name is None:continue
            param=result['atoms'][atom_name]
            row[1]=param['type'];row[6]=f"{param['charge']:.6f}"
            assigned[i]=(result,atom_name)
            provenance.append(dict(atom=i,unit=r,name=name,sourceAtom=atom_name,
                source=result['provenance']['file'],charge=param['charge'],penalty=param['penalty']))
    # Use the exact ordered CGenFF improper definitions within every mapped group.
    # Keep the original silica/interface impropers only; unassigned groups retain
    # structural carbonyl placeholders, never claim parameter completeness.
    new_impropers=[]
    for r,group in enumerate(groups):
        reverse={atom:i for i,(data,atom) in assigned.items() if i in group.values()}
        if not reverse or not result:continue
        for term in result['impropers']:
            if all(atom in reverse for atom in term):new_impropers.append([str(reverse[a]) for a in term]+['2'])
    assigned_groups={int(topology['atoms'][int(t[0])-1][2]) for t in new_impropers}
    topology['dihedrals']=[t for t in topology['dihedrals'] if t[4]!='2' or int(topology['atoms'][int(t[0])-1][2]) not in assigned_groups]+new_impropers
    return provenance
