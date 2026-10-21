"""Topology assembly and independent, versioned internal-coordinate fragments.

No polymer coordinate is read from a GRO file. Source ITP chemistry is retained.
RDKit generates capped organic fragments; identical local coordinates are translated
to every site. GROMACS, not RDKit, performs the final surface energy minimization.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

VERSION = 'independent-fragments-v1'
STEP = 0.2514

def source(key):
    from builder import sections, EXAMPLE, TEMPLATES
    return sections(EXAMPLE/(TEMPLATES[key]+'.itp'))

def fragment(key, kind, methyl=True, terminal=False):
    src=source(key)
    residue='9' if kind=='anchor' else ('29' if terminal else '11')
    old=[a[:] for a in src['atoms'] if a[2]==residue]
    ids={int(a[0]):a[4] for a in old}
    edges=[(ids[int(b[0])],ids[int(b[1])]) for b in src['bonds'] if int(b[0]) in ids and int(b[1]) in ids]
    atoms={a[4]:a for a in old}
    if kind=='repeat' and not methyl:
        # Methyl group is neutral in both source topologies. Replace it by one
        # aliphatic H (+0.09), with compensating -0.09 on the bonded alpha carbon.
        remove={'C18'}|{v if u=='C18' else u for u,v in edges if 'C18' in (u,v) and atoms[v if u=='C18' else u][1].startswith('H')}
        atoms={n:a for n,a in atoms.items() if n not in remove}
        edges=[e for e in edges if not set(e)&remove]
        atoms['C7'][1]='CG321' if terminal else 'CG311'
        atoms['C7'][6]=f'{float(atoms["C7"][6])-.09:.6f}'
        htype='HGA2' if terminal else 'HGA1'
        if terminal: atoms['H14'][1]=htype
        atoms['HA']=['0',htype,residue,'REP','HA','0','0.09','1.008']
        edges.append(('C7','HA'))
    return atoms,edges

def _generate_fragment(key, kind, methyl):
    atoms,edges=fragment(key,kind,methyl)
    names=list(atoms)
    mol=Chem.RWMol()
    for name,a in atoms.items():
        symbol='H' if float(a[7])<2 else ('N' if a[1].startswith('N') else 'O' if a[1].startswith('O') else 'S' if a[1].startswith('S') else 'C')
        atom=Chem.Atom(symbol); atom.SetNoImplicit(True)
        if a[1]=='NG3P0': atom.SetFormalCharge(1)
        mol.AddAtom(atom)
    lookup={n:i for i,n in enumerate(names)}
    for u,v in edges: mol.AddBond(lookup[u],lookup[v],Chem.BondType.SINGLE)
    for name,a in atoms.items():
        i=lookup[name]
        if a[1] in ('CG2O1','CG2O2','CG2O3','SG3O1'):
            oxy=[n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol()=='O']
            if a[1]=='SG3O1': doubles=oxy[:2]; negative=oxy[2:]
            elif a[1]=='CG2O3': doubles=oxy[:1]; negative=oxy[1:]
            else:
                doubles=[j for j in oxy if mol.GetAtomWithIdx(j).GetDegree()==1]; negative=[]
            for j in doubles: mol.GetBondBetweenAtoms(i,j).SetBondType(Chem.BondType.DOUBLE)
            for j in negative: mol.GetAtomWithIdx(j).SetFormalCharge(-1)
    if kind=='repeat':
        fixed={'C6':(0,0,0),'C7':(.889,0,1.257)}
        caps=[('C6',(.629,0,-.890)),('C7',(.260,0,2.147))]
    else:
        fixed={'C1':(0,0,-7),'C5':(.889,0,-1.257)}
        caps=[('C1',(0,0,-8.09)),('C5',(.260,0,-.367))]
    for parent,xyz in caps:
        i=mol.AddAtom(Chem.Atom('H')); mol.AddBond(lookup[parent],i,Chem.BondType.SINGLE)
        fixed[f'cap{i}']=xyz; lookup[f'cap{i}']=i
    mol=mol.GetMol(); Chem.SanitizeMol(mol)
    params=AllChem.ETKDGv3(); params.randomSeed=1947; params.useRandomCoords=True
    params.SetCoordMap({lookup[n]:Point3D(*xyz) for n,xyz in fixed.items()})
    if AllChem.EmbedMolecule(mol,params)!=0: raise ValueError('Independent fragment embedding failed')
    # Keep the four polymerization-port atoms fixed in their canonical frame.
    conf=mol.GetConformer()
    for n,xyz in fixed.items(): conf.SetAtomPosition(lookup[n],Point3D(*xyz))
    ff=AllChem.UFFGetMoleculeForceField(mol)
    for n in fixed: ff.AddFixedPoint(lookup[n])
    ff.Initialize(); ff.Minimize(maxIts=2000)
    return {n:[round(v/10,8) for v in conf.GetAtomPosition(i)] for i,n in enumerate(names)}

@lru_cache(maxsize=8)
def coordinates(key,kind,methyl=True):
    from builder import ROOT, EXAMPLE, TEMPLATES
    # Anchor chemistry is identical; names differ but generated coordinates use
    # the same graph ordering and seed. Never sample source GRO organic atoms.
    digest=hashlib.sha256((EXAMPLE/(TEMPLATES[key]+'.itp')).read_bytes()).hexdigest()[:12]
    path=ROOT/'geometry_templates'/f'{VERSION}-{key}-{kind}-{int(methyl)}-{digest}.json'
    if path.exists(): return json.loads(path.read_text())['coordinates_nm']
    result=_generate_fragment(key,kind,methyl)
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(dict(version=VERSION,method='RDKit ETKDG + UFF; fixed canonical polymerization ports; no GRO polymer coordinates',coordinates_nm=result),indent=2))
    return result

def make_topology(key,n,methyl):
    src=source(key)
    atoms=[a[:] for a in src['atoms'][:52]]
    for a in atoms: a[2]='1'
    groups=[]; edges=[]
    for b in src['bonds']:
        if int(b[0])<=52 and int(b[1])<=52: edges.append((int(b[0]),int(b[1])))
    for r in range(n+1):
        kind='anchor' if r==0 else 'repeat'
        block,bonds=fragment(key,kind,methyl,r==n and r>0)
        mapping={}
        for name,a in block.items():
            i=len(atoms)+1; mapping[name]=i
            a=a[:]; a[0]=a[5]=str(i); a[2]=str(r+2); a[3]='POL' if r==0 else ('END' if r==n else 'REP')
            atoms.append(a)
        edges.extend((mapping[u],mapping[v]) for u,v in bonds)
        edges.append((34,mapping['C1']) if r==0 else ((groups[-1]['C5'] if r==1 else groups[-1]['C7']),mapping['C6']))
        groups.append(mapping)
    adj={i:set() for i in range(1,len(atoms)+1)}
    for i,j in edges: adj[i].add(j); adj[j].add(i)
    angles=[]; torsions=set(); pairs=set()
    for j,neighbors in adj.items():
        for i,k in combinations(sorted(neighbors),2): angles.append([str(i),str(j),str(k),'5'])
    for j,k in edges:
        for i in adj[j]-{k}:
            for l in adj[k]-{j,i}:
                t=(i,j,k,l); torsions.add(min(t,t[::-1]))
    for i,_,_,l in torsions:
        if l not in adj[i] and not (adj[i]&adj[l]): pairs.add(tuple(sorted((i,l))))
    proper=[list(map(str,t))+['9'] for t in sorted(torsions)]
    # Preserve source anchor improper and instantiate carbonyl planarity for
    # every repeat, including newly added or terminal repeat units.
    anchor=groups[0]
    improper=[[str(anchor[x]) for x in ('C4','C5','N1','O35')]+['2']]
    for g in groups[1:]:
        names=('C8','C7','N2','O36') if key=='CBAA' else ('C8','C7','O36','O34')
        improper.append([str(g[x]) for x in names]+['2'])
    return dict(moleculetype=[['BRUSH','3']],atoms=atoms,bonds=[list(map(str,b))+['1'] for b in edges],pairs=[list(map(str,p))+['1'] for p in sorted(pairs)],angles=angles,dihedrals=proper+improper),groups

def make_chain(key,n,methyl,silica):
    topology,groups=make_topology(key,n,methyl)
    local_anchor=coordinates(key,'anchor',True)
    local_repeat=coordinates(key,'repeat',methyl)
    si=np.array(silica[33]['xyz'])
    # Canonical Si-C1 distance 0.185 nm. Backbone grows along +z.
    origin=si+np.array([0,0,.885])
    positions=[a['xyz'][:] for a in silica[:52]]
    for r,g in enumerate(groups):
        local=local_anchor if r==0 else local_repeat
        for name in g:
            if name=='H14':
                # Terminal cap replaces the next C6 bond, retaining tetrahedral direction.
                xyz=np.array(local['C7'])+np.array([-.0629,0,.0890])
            else: xyz=np.array(local[name])
            # Alternating sidechain orientation avoids placing adjacent bulky
            # zwitterions on the same side of an extended backbone. The same
            # deterministic sequence is used in every grafted chain.
            if r>0 and r%2==0: xyz[1]*=-1
            positions.append((origin+xyz+np.array([0,0,STEP*max(0,r-1)])).tolist())
    block=[dict(res=int(a[2]),residue=a[3],name=a[4],xyz=positions[i]) for i,a in enumerate(topology['atoms'])]
    return block,topology

def itp_text(topology):
    return '; Generated connectivity; source charges + documented alpha-H patch\n'+''.join('\n[ '+k+' ]\n'+'\n'.join(' '.join(row) for row in rows)+'\n' for k,rows in topology.items())
