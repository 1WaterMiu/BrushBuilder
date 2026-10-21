import io, json, unittest, zipfile
from builder import *

class BuilderTests(unittest.TestCase):
    def test_high_exact_coordinates(self):
        for key in TEMPLATES:
            m=build(key); sub,sup,blocks,chains,box=template(key)
            original=[a for b in blocks+chains for a in b]
            self.assertEqual([a['xyz'] for a in original],[a['xyz'] for a in m['atoms']])
            self.assertAlmostEqual(m['charge'],0,places=8)
            self.assertEqual(m['chainCount'],20)
            self.assertEqual(len(m['frozen']),140*54+20*52)

    def test_low_preserves_silica_and_restores_charge(self):
        for key in TEMPLATES:
            sub,sup,blocks,chains,box=template(key)
            high=build(key); low=build(key,'low')
            self.assertAlmostEqual(low['graftDensity']*2,high['graftDensity'])
            self.assertEqual(low['subCount'],150)
            self.assertEqual(len(low['frozen']),150*54+10*52)
            self.assertAlmostEqual(low['charge'],0,places=8)
            for chain in chains:
                new=restore_sub(chain,sub,box)
                lookup={a['name']:a for a in new}
                for a in chain[:52]: self.assertEqual(a['xyz'],lookup[a['name']]['xyz'])
                self.assertEqual(len(new),54)
                oh=periodic_delta(lookup['H4']['xyz'],lookup['O34']['xyz'],box)
                self.assertTrue(.08<math.sqrt(sum(v*v for v in oh))<.11)

    def test_patterns_and_invalid_input(self):
        for seed in range(20): self.assertEqual(len(set(build('SBMA','low',seed)['keptSites'])),10)
        for kwargs in [dict(repeats=10),dict(key='unknown'),dict(density='other'),dict(seed=20)]:
            with self.assertRaises(ValueError): build(**kwargs)

    def test_export_consistency(self):
        for key in TEMPLATES:
            m=build(key,'low')
            with zipfile.ZipFile(io.BytesIO(archive(m))) as z:
                lines=z.read('surface.gro').decode().splitlines()
                self.assertEqual(int(lines[1]),len(m['atoms']))
                self.assertEqual(len(lines),len(m['atoms'])+3)
                report=json.loads(z.read('manifest.json'))
                self.assertEqual(report['atomCount'],len(m['atoms']))
                self.assertIn('SUB 150',z.read('topol.top').decode())
                self.assertIn('freezedim = Y Y Y',z.read('freeze.mdp').decode())
                self.assertEqual(z.read(TEMPLATES[key]+'.itp'),(EXAMPLE/(TEMPLATES[key]+'.itp')).read_bytes())
                for name in z.namelist():
                    if name.endswith(('.itp','.top')):
                        includes=re.findall(r'^\s*#include\s+"([^"]+)"',z.read(name).decode(),re.M)
                        for dependency in includes:
                            self.assertIn(dependency,z.namelist(),f'{name} requires {dependency}')

if __name__=='__main__': unittest.main()
