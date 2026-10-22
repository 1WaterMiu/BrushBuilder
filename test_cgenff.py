import unittest,json,hashlib,tempfile,shutil
from pathlib import Path
from builder import build,archive
from cgenff_library import ROOT,load

class CGenFFTests(unittest.TestCase):
    def test_compact_gui_archive_reimports_without_full_download(self):
        from tools import import_charmm_gui as importer
        from unittest.mock import patch
        row=json.loads((ROOT/'charmm-gui-audit.json').read_text())[0]
        stem=row['file']
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder);(target/'inputs').mkdir()
            for ext in ('.json','.sdf'):shutil.copyfile(ROOT/'inputs'/(stem+ext),target/'inputs'/(stem+ext))
            archive=target/Path(row['sourceArchive']).name
            shutil.copyfile(ROOT/row['sourceArchive'],archive)
            with patch.object(importer,'ROOT',target):report=importer.import_archive(archive)
            self.assertTrue(report['bondedCoverage']['complete'])
            self.assertEqual(json.loads((target/'parsed'/(stem+'.json')).read_text())['atoms'],load(stem)['atoms'])

    def test_charmm_gui_provenance_and_graphs(self):
        rows=json.loads((ROOT/'charmm-gui-audit.json').read_text())
        self.assertGreaterEqual(len(rows),2)
        for r in rows:
            self.assertTrue(r['orderedChemicalGraphExact'])
            self.assertTrue(r['bondedCoverage']['complete'])
            self.assertAlmostEqual(r['charge'],0,places=8)
            self.assertEqual(hashlib.sha256((ROOT/r['sourceArchive']).read_bytes()).hexdigest(),r['sourceArchiveSHA256'])
            data=load(r['file'])
            self.assertEqual(set(r['guiToInputAtomMap'].values()),set(data['atoms']))

    def test_result_provenance_and_coverage(self):
        rows=json.loads((ROOT/'audit.json').read_text())
        self.assertEqual(len(rows),35)
        for r in rows:
            self.assertTrue(r['atomMapExact'] and r['bondGraphExact'])
            self.assertTrue(r['bondedCoverage']['complete'])
            self.assertAlmostEqual(r['charge'],0,places=8)
            self.assertEqual(hashlib.sha256((ROOT/'results'/(r['file']+'.zip')).read_bytes()).hexdigest(),r['resultSHA256'])

    def test_complete_charge_transfer_and_terminal_impropers(self):
        raw=load('CB_amide_3_m1_anchor5')
        for n in (5,20,60):
            model=build(repeats=n,density='low')
            self.assertTrue(model['chemistry']['parameterization']['chargesComplete'])
            self.assertAlmostEqual(model['charge'],0,places=8)
            top=model['topology']
            self.assertEqual(sum(t[4]=='2' for t in top['dihedrals']),2*n+1)
            for row in top['atoms'][52:]:
                r=int(row[2])-2;name=row[4]
                if r==0 and name in ('C1','H4','H5'):continue
                unit=-1 if r==0 else r if r<3 else 5 if r==n else 4 if r==n-1 else 3
                source=raw['mapping'][(6,'END_H') if name=='H14' else (unit,name)]
                self.assertEqual(float(row[6]),raw['atoms'][source]['charge'])
                self.assertEqual(row[1],raw['atoms'][source]['type'])
            self.assertFalse(model['parameterized'])
            with self.assertRaises(ValueError):archive(model)

    def test_no_template_charge_fallback(self):
        for opts in (dict(key='SBMA'),dict(repeats=2),dict(methyl=False)):
            model=build(**opts)
            self.assertIsNone(model['charge'])
            self.assertFalse(model['parameterized'])
            self.assertTrue(any(a[6]=='UNASSIGNED' for a in model['topology']['atoms'][52:]))

if __name__=='__main__':unittest.main()
