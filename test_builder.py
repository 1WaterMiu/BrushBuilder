import io,json,re,unittest,zipfile
from collections import Counter
from unittest.mock import patch
from builder import *
from polymer import make_topology,make_chain,coordinates

class BuilderTests(unittest.TestCase):
    def test_silica_invariance(self):
        for key in TEMPLATES:
            sub,sup,blocks,chains,box=template(key)
            original=Counter(tuple(a['xyz']) for b in blocks for a in b)
            original.update(tuple(a['xyz']) for c in chains for a in c[:52])
            for density in ['high','low']:
                m=build(key,density,repeats=3,methyl=False)
                fixed=Counter(tuple(m['atoms'][i-1]['xyz']) for i in m['frozen'])
                self.assertFalse(original-fixed)
                self.assertEqual(len(m['frozen']),8600 if density=='high' else 8616)
                self.assertAlmostEqual(m['charge'],0,places=8)
                self.assertEqual(m['chainCount'],20 if density=='high' else 12)
                self.assertEqual(m['removedSites'],[] if density=='high' else [2,4,7,9,12,14,17,19])

    def test_graph_and_terminal_regeneration(self):
        for key in TEMPLATES:
            for methyl in [True,False]:
                last=None
                for n in [1,2,10,20,60]:
                    top,groups=make_topology(key,n,methyl)
                    count=len(top['atoms'])
                    self.assertEqual(len(groups),n+1)
                    self.assertAlmostEqual(sum(float(a[6]) for a in top['atoms']),0,places=8)
                    self.assertEqual(sum(r[4]=='2' for r in top['dihedrals']),n+1)
                    self.assertEqual(sum(a[4]=='H14' for a in top['atoms']),1)
                    self.assertEqual(sum(a[4]=='HA' for a in top['atoms']),0 if methyl else n)
                    self.assertEqual(sum(a[4]=='C18' for a in top['atoms']),n if methyl else 0)
                    adj={i:set() for i in range(1,count+1)}
                    for row in top['bonds']:
                        a,b=map(int,row[:2]);adj[a].add(b);adj[b].add(a)
                    for i,a in enumerate(top['atoms'],1):
                        if a[4]=='C7':self.assertEqual(len(adj[i]),4)
                    for section,arity in [('bonds',2),('pairs',2),('angles',3),('dihedrals',4)]:
                        for row in top[section]:
                            ids=list(map(int,row[:arity]))
                            self.assertTrue(all(1<=x<=count for x in ids))
                            self.assertEqual(len(set(ids)),arity)
                    if last:self.assertGreater(count,last)
                    last=count

    def test_independent_coordinates_and_identical_copies(self):
        for key in TEMPLATES:
            _,_,_,chains,_=template(key)
            mutated=[dict(a,xyz=[100,100,100]) if i>=52 else a for i,a in enumerate(chains[0])]
            for methyl in [True,False]:
                a,top=make_chain(key,4,methyl,chains[0])
                b,_=make_chain(key,4,methyl,mutated)
                c,_=make_chain(key,4,methyl,chains[1])
                self.assertEqual(a,b)
                shift=[chains[1][33]['xyz'][j]-chains[0][33]['xyz'][j] for j in range(3)]
                for left,right in zip(a[52:],c[52:]):
                    for j in range(3):self.assertAlmostEqual(right['xyz'][j]-left['xyz'][j],shift[j],places=12)
                shorter,_=make_chain(key,3,methyl,chains[0])
                # Anchor and first two repeats do not depend on total chain length.
                stop=next(i for i,x in enumerate(shorter) if x['res']==5)
                self.assertEqual(a[:stop],shorter[:stop])

    def test_package(self):
        m=build('SBMA','low',repeats=4,methyl=False)
        with zipfile.ZipFile(io.BytesIO(archive(m))) as z:
            for name in ['BRUSH.itp','em.mdp','nvt.mdp','run_md.sh','update_index.py','tip3p.itp','ions.itp']:
                self.assertIn(name,z.namelist())
            for name in z.namelist():
                if name.endswith(('.itp','.top')):
                    for inc in re.findall(r'^\s*#include\s+"([^"]+)"',z.read(name).decode(),re.M):self.assertIn(inc,z.namelist())
            self.assertIn('BRUSH 12',z.read('topol.top').decode())
            self.assertIn('freezedim = Y Y Y',z.read('em.mdp').decode())
            self.assertIn('freezedim = Y Y Y',z.read('nvt.mdp').decode())
            self.assertEqual(int(z.read('surface.gro').decode().splitlines()[1]),len(m['atoms']))

    def test_invalid(self):
        for kwargs in [dict(repeats=0),dict(repeats=61),dict(repeats=2.5),dict(key='x'),dict(methyl='false')]:
            with self.assertRaises(ValueError):build(**kwargs)

if __name__=='__main__':unittest.main()
