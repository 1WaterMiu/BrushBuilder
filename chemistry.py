"""Explicit chemical graphs. New combinations deliberately have no partial charges."""
from functools import lru_cache
from collections import Counter
import json
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

FAMILIES=['CB','SB','PB','sulfabetaine','PC','TMAO']
VERSION='head-graphs-v1'
SOURCES={
 'SB':dict(url='https://pubs.rsc.org/en/content/articlepdf/2026/lp/d6lp00019c',title='RSC Applied Polymers 2026, Scheme 1',scope='SBA / SBMA / SBAA / SBMAA; standard spacer 3. Other spacers are enumerated homologues.'),
 'CB':dict(url='https://patents.google.com/patent/US20140275614A1/en',title='CBMA, CBAA and CBMAA synthesis definitions',scope='Acrylamide versus methacrylamide distinction. Exact spacer is recorded separately; source ITP supplies the spacer-3 reference graph.'),
 'PB':dict(url='https://pubmed.ncbi.nlm.nih.gov/1880723/',title='Trialkylammonio alkyl hydrogen phosphate phosphobetaines',scope='Phosphate-type N+-(CH2)n-O-P(=O)(O-)-OH head motif, as confirmed by the user. Polymerizable combinations are enumerated derivatives, not all independently literature-verified monomers.'),
 'sulfabetaine':dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6631746/',title='Sulfabetaine monomers and polymers (2019)',scope='Sulfate ester head, distinct from sulfobetaine sulfonate; methacrylate and methacrylamide architectures. Additional homologues are enumerated derivatives.'),
 'PC':dict(url='https://pubchem.ncbi.nlm.nih.gov/compound/128934',title='MPC / 2-methacryloyloxyethyl phosphorylcholine',scope='Exact MPC structure for ester + alpha-methyl. Other linker/backbone combinations are explicitly labelled derivatives.'),
 'TMAO':dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6570511/',title='TMAO-derived acrylamide monomer (2019)',scope='Oxidized dimethylaminopropyl acrylamide; N+-O- direct bond. Other backbones/linkers are derivatives.'),
}

def nomenclature(spec,methyl):
    head,link,n=spec
    suffix=('MA' if methyl else 'A')+('A' if link=='amide' else '')
    polymerizable=('methacrylate' if methyl else 'acrylate') if link=='ester' else ('methacrylamide' if methyl else 'acrylamide')
    if head in ('CB','SB'):name=head+suffix
    elif head=='PC' and link=='ester' and methyl:name='MPC'
    else:name={'PB':'Phosphate-type phosphobetaine','sulfabetaine':'Sulfabetaine','PC':'Phosphorylcholine-derived','TMAO':'TMAO-derived'}[head]+' '+polymerizable
    return dict(name=name,polymerizableGroup=polymerizable,displayName=name+(f' · spacer {n}' if n is not None else ''),source=SOURCES[head])

def normalize(key='CBAA',head_group=None,linker=None,spacer=None):
    head=head_group or ('SB' if key=='SBMA' else 'CB')
    link=linker or ('ester' if key=='SBMA' else 'amide')
    if link=='amid':link='amide'
    if head not in FAMILIES or link not in ('ester','amide'):raise ValueError('Unknown head group or linker')
    if head in ('PC','TMAO'):
        if spacer not in (None,'', 'fixed'):raise ValueError('PC and TMAO have fixed charge separation; no spacer selection')
        n=None
    else:
        raw=3 if spacer is None else spacer
        if isinstance(raw,bool) or str(raw)!=str(int(raw)) or not 1<=int(raw)<=4:raise ValueError('Charge spacer must be an integer from 1 to 4')
        n=int(raw)
    return (head,link,n)

def reference(spec):
    return {('CB','amide',3):'CBAA',('SB','ester',3):'SBMA'}.get(spec)

def repeat_smiles(spec,methyl):
    head,link,n=spec
    # This connector is NOT the adjustable charge spacer. Preserve the two
    # supplied architectures: ester-O-CH2-CH2 and amide-NH-(CH2)3.
    connector='[O:34]CC' if link=='ester' else '[NH:22]CCC'
    if head=='PC':tail='OP(=O)([O-])OCC[N+](C)(C)C'
    elif head=='TMAO':tail='[N+](C)(C)[O-]'
    else:
        end={'CB':'C(=O)[O-]','SB':'S(=O)(=O)[O-]',
             'PB':'OP(=O)([O-])O','sulfabetaine':'OS(=O)(=O)[O-]'}[head]
        tail='[N+](C)(C)'+'C'*n+end
    alpha='[C:7](C)' if methyl else '[CH:7]'
    return '[*:1][CH2:6]'+alpha+'([C:8](=[O:36])'+connector+tail+')[*:2]'

@lru_cache(maxsize=80)
def graph(spec,methyl):
    mol=Chem.AddHs(Chem.MolFromSmiles(repeat_smiles(spec,methyl)))
    special={6:'C6',7:'C7',8:'C8',36:'O36',34:'O34',22:'N2'}
    names={}
    for a in mol.GetAtoms():
        if a.GetAtomicNum()==0:names[a.GetIdx()]='PORT'+str(a.GetAtomMapNum())
        else:names[a.GetIdx()]=special.get(a.GetAtomMapNum(),a.GetSymbol()+str(100+a.GetIdx()))
    return mol,names

def repeat_fragment(spec,methyl,terminal=False):
    mol,names=graph(spec,methyl)
    atoms={names[a.GetIdx()]:['0','UNASSIGNED','0','REP',names[a.GetIdx()],'0','UNASSIGNED',str(a.GetMass())]
           for a in mol.GetAtoms() if a.GetAtomicNum()}
    edges=[(names[b.GetBeginAtomIdx()],names[b.GetEndAtomIdx()]) for b in mol.GetBonds()
           if b.GetBeginAtom().GetAtomicNum() and b.GetEndAtom().GetAtomicNum()]
    if terminal:
        atoms['H14']=['0','UNASSIGNED','0','END','H14','0','UNASSIGNED','1.008']
        edges.append(('C7','H14'))
    return atoms,edges

def capped_molecule(spec,methyl):
    original,names=graph(spec,methyl)
    mol=Chem.RWMol(original)
    for a in mol.GetAtoms():
        if not a.GetAtomicNum():a.SetAtomicNum(1);a.SetNoImplicit(True)
    result=mol.GetMol();Chem.SanitizeMol(result)
    return result,names

def monomer(spec,methyl):
    mol,_=graph(spec,methyl);rw=Chem.RWMol(mol)
    for i in sorted([a.GetIdx() for a in rw.GetAtoms() if not a.GetAtomicNum()],reverse=True):rw.RemoveAtom(i)
    ids={a.GetAtomMapNum():a.GetIdx() for a in rw.GetAtoms() if a.GetAtomMapNum()}
    rw.GetBondBetweenAtoms(ids[6],ids[7]).SetBondType(Chem.BondType.DOUBLE)
    result=rw.GetMol();Chem.SanitizeMol(result)
    for a in result.GetAtoms():a.SetAtomMapNum(0)
    return result

@lru_cache(maxsize=80)
def coordinates(spec,methyl):
    from builder import ROOT
    head,link,n=spec
    path=ROOT/'geometry_templates'/f'{VERSION}-{head}-{link}-{n or 0}-{int(methyl)}.json'
    if path.exists():
        cached=json.loads(path.read_text())
        if cached.get('smiles')==repeat_smiles(spec,methyl):return cached['coordinates_nm']
    mol,names=capped_molecule(spec,methyl)
    lookup={v:k for k,v in names.items()}
    fixed={'C6':(0,0,0),'C7':(.889,0,1.257),'PORT1':(.629,0,-.890),'PORT2':(.260,0,2.147)}
    params=AllChem.ETKDGv3();params.randomSeed=1947;params.useRandomCoords=True
    params.SetCoordMap({lookup[n]:Point3D(*xyz) for n,xyz in fixed.items()})
    if AllChem.EmbedMolecule(mol,params)!=0:raise ValueError('Head-group geometry embedding failed')
    conf=mol.GetConformer()
    for name,xyz in fixed.items():conf.SetAtomPosition(lookup[name],Point3D(*xyz))
    ff=AllChem.UFFGetMoleculeForceField(mol)
    if ff is None:raise ValueError('Geometry force field unavailable')
    for name in fixed:ff.AddFixedPoint(lookup[name])
    ff.Initialize();ff.Minimize(maxIts=2000)
    result={name:[round(v/10,8) for v in conf.GetAtomPosition(i)] for i,name in names.items() if not name.startswith('PORT')}
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(dict(version=VERSION,smiles=repeat_smiles(spec,methyl),coordinates_nm=result),indent=2))
    return result

def metadata(spec,methyl):
    mol,names=graph(spec,methyl)
    vinyl=monomer(spec,methyl)
    heavy=Chem.RemoveHs(vinyl)
    for a in heavy.GetAtoms():a.SetAtomMapNum(0)
    drawer=rdMolDraw2D.MolDraw2DSVG(620,210)
    drawer.drawOptions().clearBackground=False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,heavy)
    drawer.FinishDrawing()
    elements=Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum())
    ref=reference(spec)
    return dict(headGroup=spec[0],linker=spec[1],spacer=spec[2],nomenclature=nomenclature(spec,methyl),
        monomerSmiles=Chem.MolToSmiles(heavy),monomerFormula=rdMolDescriptors.CalcMolFormula(vinyl),
        connectorCH2=2 if spec[1]=='ester' else 3,
        referenceTemplate=ref,parameterized=bool(ref),formalCharge=Chem.GetFormalCharge(mol),
        matchesOriginalItp=bool(ref) and methyl,
        repeatAtomCount=sum(elements.values()),repeatElements=dict(elements),
        repeatSmiles=repeat_smiles(spec,methyl),structureSvg=drawer.GetDrawingText(),
        parameterStatus='source-template' if ref else 'unassigned',
        chargePath={'CB':f'N+–(CH2){spec[2]}–C(=O)O−','SB':f'N+–(CH2){spec[2]}–S(=O)2O−',
          'PB':f'N+–(CH2){spec[2]}–O–P(=O)(O−)–OH',
          'sulfabetaine':f'N+–(CH2){spec[2]}–O–S(=O)2–O−',
          'PC':'P(=O)(O−)–O–CH2–CH2–N+(CH3)3','TMAO':'N+–O−'}[spec[0]])

def source_repeat_graph(key):
    from polymer import fragment
    atoms,edges=fragment(key,'repeat',True)
    rw=Chem.RWMol();ids={}
    for name,a in atoms.items():
        element='H' if float(a[7])<2 else name[0]
        atom=Chem.Atom(element);atom.SetNoImplicit(True)
        if a[1]=='NG3P0':atom.SetFormalCharge(1)
        ids[name]=rw.AddAtom(atom)
    for u,v in edges:rw.AddBond(ids[u],ids[v],Chem.BondType.SINGLE)
    for name,a in atoms.items():
        if a[1] not in ('CG2O1','CG2O2','CG2O3','SG3O1'):continue
        oxy=[x.GetIdx() for x in rw.GetAtomWithIdx(ids[name]).GetNeighbors() if x.GetSymbol()=='O']
        if a[1]=='SG3O1':double,negative=oxy[:2],oxy[2:]
        elif a[1]=='CG2O3':double,negative=oxy[:1],oxy[1:]
        else:double,negative=[j for j in oxy if rw.GetAtomWithIdx(j).GetDegree()==1],[]
        for j in double:rw.GetBondBetweenAtoms(ids[name],j).SetBondType(Chem.BondType.DOUBLE)
        for j in negative:rw.GetAtomWithIdx(j).SetFormalCharge(-1)
    for port,name in [(1,'C6'),(2,'C7')]:
        atom=Chem.Atom(0);i=rw.AddAtom(atom);rw.AddBond(ids[name],i,Chem.BondType.SINGLE)
    mol=rw.GetMol();Chem.SanitizeMol(mol)
    return mol

def template_comparison(key):
    """Exact named atom/bond comparison over the entire N20 grafted molecule."""
    from polymer import source,make_topology
    old=source(key);new,_=make_topology(key,20,True)
    def atom_keys(top):
        residues=list(dict.fromkeys(a[2] for a in top['atoms']))
        return {a[0]:(residues.index(a[2]),a[4]) for a in top['atoms']}
    a,b=atom_keys(old),atom_keys(new)
    old_parameters={a[r[0]]:(r[1],r[6],r[7]) for r in old['atoms']}
    new_parameters={b[r[0]]:(r[1],r[6],r[7]) for r in new['atoms']}
    oldedges={tuple(sorted((a[r[0]],a[r[1]]))) for r in old['bonds']}
    newedges={tuple(sorted((b[r[0]],b[r[1]]))) for r in new['bonds']}
    spec=('CB','amide',3) if key=='CBAA' else ('SB','ester',3)
    generated=Chem.Mol(graph(spec,True)[0])
    for x in generated.GetAtoms():x.SetAtomMapNum(0)
    reference_graph=source_repeat_graph(key)
    identical=Chem.MolToSmiles(generated)==Chem.MolToSmiles(reference_graph)
    return dict(template=key,sourceAtomCount=len(a),generatedAtomCount=len(b),
       sourceBondCount=len(oldedges),generatedBondCount=len(newedges),
       missingAtoms=sorted(set(a.values())-set(b.values())),extraAtoms=sorted(set(b.values())-set(a.values())),
       missingBonds=sorted(oldedges-newedges),extraBonds=sorted(newedges-oldedges),
       explicitHydrogenRepeatGraphIdentical=identical,
       atomTypesChargesMassesIdentical=old_parameters==new_parameters,
       passed=set(a.values())==set(b.values()) and oldedges==newedges and identical and old_parameters==new_parameters)
