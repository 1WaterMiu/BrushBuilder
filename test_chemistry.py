import io,json,unittest,zipfile
from collections import Counter
from rdkit import Chem
from chemistry import FAMILIES,normalize,graph,repeat_fragment,template_comparison,reference,nomenclature,metadata
from structure_export import archive_structure,monomer
from builder import build,archive,template
from polymer import make_chain

class ChemistryTests(unittest.TestCase):
    def test_full_reference_graphs_include_all_hydrogens(self):
        for key in ['CBAA','SBMA']:
            report=template_comparison(key)
            self.assertTrue(report['passed'],report)

    def test_all_combinations_valence_counts_coordinates_and_caps(self):
        silica=template('CBAA')[3][0]
        checked=0
        for head in FAMILIES:
            for link in ['ester','amide']:
                previous={}
                for n in ([None] if head in ('PC','TMAO') else range(1,5)):
                    spec=normalize(head_group=head,linker=link,spacer=n)
                    for methyl in [False,True]:
                        with self.subTest(spec=spec,methyl=methyl):
                            mol,names=graph(spec,methyl)
                            self.assertEqual(Chem.GetFormalCharge(mol),0)
                            self.assertEqual(sum(a.GetNumImplicitHs() for a in mol.GetAtoms()),0)
                            atoms,bonds=repeat_fragment(spec,methyl)
                            if methyl in previous:self.assertEqual(len(atoms)-previous[methyl],3)
                            previous[methyl]=len(atoms)
                            chain,top=make_chain('CBAA' if link=='amide' else 'SBMA',2,methyl,silica,spec)
                            self.assertEqual(len(chain),52+22+2*len(atoms)+1)
                            self.assertEqual([a['xyz'] for a in chain[:52]],[a['xyz'] for a in silica[:52]])
                            counts=Counter(i for b in top['bonds'] for i in b[:2])
                            for a in top['atoms'][74:]:
                                if a[4].startswith('H'):self.assertEqual(counts[a[0]],1)
                            vm=monomer(spec,methyl)
                            self.assertEqual(vm.GetNumAtoms(),len(atoms))
                            checked+=1
        self.assertEqual(checked,72)

    def test_charge_free_new_combinations_and_structure_export(self):
        model=build(head_group='TMAO',linker='amide',repeats=2,methyl=True,density='low')
        self.assertIsNone(model['charge'])
        self.assertFalse(model['parameterized'])
        for a in model['topology']['atoms'][74:]:self.assertEqual(a[6],'UNASSIGNED')
        with self.assertRaises(ValueError):archive(model)
        with zipfile.ZipFile(io.BytesIO(archive_structure(model))) as z:
            self.assertNotIn('topol.top',z.namelist())
            self.assertIn('monomer.sdf',z.namelist())
            m=Chem.MolFromMolBlock(z.read('monomer.sdf').decode(),removeHs=False)
            self.assertIsNotNone(m)
            self.assertEqual(Chem.GetFormalCharge(m),0)

    def test_spacer_and_linker_validation(self):
        for head in ['PC','TMAO']:
            with self.assertRaises(ValueError):normalize(head_group=head,spacer=3)
        for bad in [0,5,True,2.5]:
            with self.assertRaises(ValueError):normalize(spacer=bad)
        self.assertEqual(normalize(linker='amid')[1],'amide')

    def test_literature_nomenclature_and_mpc_graph(self):
        for link,methyl,name in [('ester',False,'SBA'),('ester',True,'SBMA'),('amide',False,'SBAA'),('amide',True,'SBMAA')]:
            self.assertEqual(nomenclature(('SB',link,3),methyl)['name'],name)
        self.assertEqual(nomenclature(('CB','amide',3),True)['name'],'CBMAA')
        self.assertEqual(nomenclature(('CB','amide',3),False)['name'],'CBAA')
        m=monomer(('PC','ester',None),True)
        pubchem=Chem.MolFromSmiles('CC(=C)C(=O)OCCOP(=O)([O-])OCC[N+](C)(C)C')
        self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(m)),Chem.MolToSmiles(pubchem))
        self.assertEqual(metadata(('PC','ester',None),True)['monomerFormula'],'C11H22NO6P')

if __name__=='__main__':unittest.main()
