"""Reproducible explicit-H inputs; never import template partial charges."""
import sys,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors
from openbabel import openbabel as ob
from chemistry import FAMILIES,monomer,graph,nomenclature

ROOT=Path(__file__).resolve().parents[1]/'cgenff'

def oligomer(spec,methyl,n=5):
    fragment,names=graph(spec,methyl)
    combined=Chem.RWMol();mapping=[];ports=[]
    for r in range(n):
        ids={}
        for a in fragment.GetAtoms():
            if not a.GetAtomicNum():continue
            atom=Chem.Atom(a);atom.SetAtomMapNum(0)
            ids[a.GetIdx()]=combined.AddAtom(atom)
            mapping.append({'unit':r+1,'name':names[a.GetIdx()]})
        for b in fragment.GetBonds():
            u,v=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            if u in ids and v in ids:combined.AddBond(ids[u],ids[v],b.GetBondType())
        ports.append({names[i]:j for i,j in ids.items() if names[i] in ('C6','C7')})
    for r in range(n-1):combined.AddBond(ports[r]['C7'],ports[r+1]['C6'],Chem.BondType.SINGLE)
    for i,name in [(ports[0]['C6'],'START_H'),(ports[-1]['C7'],'END_H')]:
        h=combined.AddAtom(Chem.Atom('H'));combined.AddBond(i,h,Chem.BondType.SINGLE)
        mapping.append({'unit':0 if name=='START_H' else n+1,'name':name})
    mol=combined.GetMol();Chem.SanitizeMol(mol)
    return mol,mapping

def anchored_oligomer(spec,methyl,n=5):
    from polymer import fragment
    atoms,edges=fragment('CBAA','anchor')
    mol,mapping=oligomer(spec,methyl,n)
    rw=Chem.RWMol(mol)
    start=next(i for i,x in enumerate(mapping) if x['name']=='START_H')
    rw.RemoveAtom(start);mapping.pop(start)
    first=next(i for i,x in enumerate(mapping) if x['unit']==1 and x['name']=='C6')
    ids={}
    for name,row in atoms.items():
        atom=Chem.Atom('H' if float(row[7])<2 else name[0]);atom.SetNoImplicit(True)
        ids[name]=rw.AddAtom(atom);mapping.append(dict(unit=-1,name=name))
    for a,b in edges:rw.AddBond(ids[a],ids[b],Chem.BondType.DOUBLE if {a,b}=={'C4','O35'} else Chem.BondType.SINGLE)
    rw.AddBond(ids['C5'],first,Chem.BondType.SINGLE)
    cap=rw.AddAtom(Chem.Atom('H'));rw.AddBond(ids['C1'],cap,Chem.BondType.SINGLE)
    mapping.append(dict(unit=-2,name='INTERFACE_CAP'))
    result=rw.GetMol();Chem.SanitizeMol(result)
    return result,mapping

def write_input(mol,path,mapping):
    mol=Chem.Mol(mol)
    for a in mol.GetAtoms():a.SetAtomMapNum(0)
    params=AllChem.ETKDGv3();params.randomSeed=1947;params.useRandomCoords=True
    if AllChem.EmbedMolecule(mol,params):raise ValueError(f'Embedding failed: {path.name}')
    convergence=AllChem.UFFOptimizeMolecule(mol,maxIters=3000)
    mol.SetProp('_Name',path.stem)
    sdf=Chem.MolToMolBlock(mol)
    path.with_suffix('.sdf').write_text(sdf+'\n$$$$\n',encoding='utf-8')
    conversion=ob.OBConversion();conversion.SetInAndOutFormats('mol','mol2')
    target=ob.OBMol();assert conversion.ReadString(target,sdf)
    # Suppress automatic Gasteiger charges: CGenFF must assign all charges.
    target.SetAutomaticPartialCharge(False)
    for atom in ob.OBMolAtomIter(target):atom.SetPartialCharge(0.)
    lines=conversion.WriteString(target).splitlines();lines[4]='USER_CHARGES'
    section=''; atom_index=0
    for i,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):section=line;continue
        if section=='@<TRIPOS>ATOM' and line.strip():
            fields=line.split();a=mol.GetAtomWithIdx(atom_index)
            fields[1]=a.GetSymbol()+str(atom_index+1)
            fields[7]='MOL';fields[8]=f'{a.GetFormalCharge():.6f}'
            lines[i]=' '.join(fields)
            mapping[atom_index].update(cgenffName=fields[1],element=a.GetSymbol(),formalCharge=a.GetFormalCharge())
            atom_index+=1
    assert atom_index==mol.GetNumAtoms()
    text='\n'.join(lines)+'\n';path.write_text(text,encoding='utf-8')
    # Verify the exported file preserves every atom, H and bond order.
    reverse=ob.OBConversion();reverse.SetInAndOutFormats('mol2','mol')
    recovered=ob.OBMol();assert reverse.ReadString(recovered,text)
    check=Chem.MolFromMolBlock(reverse.WriteString(recovered),sanitize=True,removeHs=False)
    if check is not None:Chem.RemoveStereochemistry(check)
    Chem.RemoveStereochemistry(mol)
    if check is None or Chem.MolToSmiles(check)!=Chem.MolToSmiles(mol):
        raise ValueError(f'MOL2 graph round-trip failed: {path.name}')
    result=dict(atoms=mol.GetNumAtoms(),bonds=mol.GetNumBonds(),formalCharge=Chem.GetFormalCharge(mol),
        formula=rdMolDescriptors.CalcMolFormula(mol),smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
        uffConverged=convergence==0,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),atomMap=mapping)
    path.with_suffix('.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return {k:v for k,v in result.items() if k!='atomMap'}

def main():
    directory=ROOT/'inputs';directory.mkdir(parents=True,exist_ok=True)
    manifest=ROOT/'manifest.json'
    extras=[r for r in json.loads(manifest.read_text()) if r['kind'] not in ('monomer','pentamer')] if manifest.exists() else []
    records=[]
    for head in FAMILIES:
        for link in ('ester','amide'):
            for spacer in ([None] if head in ('PC','TMAO') else range(1,5)):
                for methyl in (False,True):
                    spec=(head,link,spacer);key=f'{head}_{link}_{spacer or 0}_m{int(methyl)}'
                    for kind in ('monomer','pentamer'):
                        path=directory/f'{key}_{kind}.mol2'
                        if kind=='monomer':
                            mol=monomer(spec,methyl)
                            mapping=[{'unit':1,'name':f'atom{i+1}'} for i in range(mol.GetNumAtoms())]
                        else:mol,mapping=oligomer(spec,methyl)
                        if path.exists() and path.with_suffix('.json').exists():
                            record=json.loads(path.with_suffix('.json').read_text());record.pop('atomMap')
                        else:record=write_input(mol,path,mapping)
                        records.append(dict(key=key,kind=kind,spec=spec,methyl=methyl,file=path.name,**record))
                        print(path.name,record['atoms'],flush=True)
                    manifest.write_text(json.dumps(records+extras,indent=2),encoding='utf-8')
    print('Complete:',len(records),'inputs',flush=True)

if __name__=='__main__':main()
