"""Charge-free structure exports, intentionally distinct from MD-ready packages."""
import io,json,zipfile
from rdkit import Chem
from rdkit.Chem import AllChem
from chemistry import graph,capped_molecule,monomer

def sdf(mol,title):
    mol=Chem.Mol(mol)
    for a in mol.GetAtoms():a.SetAtomMapNum(0)
    mol.SetProp('_Name',title)
    params=AllChem.ETKDGv3();params.randomSeed=1947
    if AllChem.EmbedMolecule(mol,params)!=0:raise ValueError('Monomer embedding failed')
    AllChem.UFFOptimizeMolecule(mol,maxIters=1000)
    return Chem.MolToMolBlock(mol)+'\n$$$$\n'

def archive_structure(model):
    meta=model['chemistry'];spec=(meta['headGroup'],meta['linker'],meta['spacer'])
    memory=io.BytesIO()
    with zipfile.ZipFile(memory,'w',zipfile.ZIP_DEFLATED) as z:
        lines=['BrushBuilder STRUCTURE ONLY - not an MD input',str(len(model['atoms']))]
        for i,a in enumerate(model['atoms'],1):
            x,y,v=a['xyz'];lines.append(f"{a['res']%100000:5d}{a['residue']:<5.5}{a['name']:>5.5}{i%100000:5d}{x:8.3f}{y:8.3f}{v:8.3f}")
        lines.append(' '.join(f'{v:.5f}' for v in model['box']))
        z.writestr('surface.gro','\n'.join(lines)+'\n')
        z.writestr('monomer.sdf',sdf(monomer(spec,model['methyl']),'Vinyl monomer; explicit H; formal charges only'))
        z.writestr('repeat-capped.sdf',sdf(capped_molecule(spec,model['methyl'])[0],'Hydrogen-capped repeat; explicit H; formal charges only'))
        z.writestr('monomer.svg',meta['structureSvg'])
        z.writestr('structure.json',json.dumps({k:v for k,v in model.items() if k!='chemistry'},indent=2))
        z.writestr('chemistry.json',json.dumps({k:v for k,v in meta.items() if k!='structureSvg'},indent=2))
        z.writestr('README.txt','STRUCTURE ONLY. Not an MD-ready package.\nNew chemistry has no assigned partial charges or force-field atom types. Formal N+/O- charges in SDF describe chemistry, not MD partial charges.\nmonomer.sdf retains the vinyl double bond; repeat-capped.sdf is a saturated H-capped repeat. Both contain explicit hydrogens. Surface geometry uses polymerized single bonds.\nstructure.json bonds use zero-based global atom indices; topology uses one-based local indices. Angle/torsion entries describe graph connectivity only; they are not newly fitted parameters.\nSilica coordinates remain fixed.\n')
    return memory.getvalue()
