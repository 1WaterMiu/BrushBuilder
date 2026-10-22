"""Compare CHARMM-GUI effective parameters to the original CGenFF web result."""
from pathlib import Path
import sys,json,hashlib,tarfile
from itertools import combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.audit_cgenff import parse
import networkx as nx
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]/'cgenff'
JOB=ROOT/'charmm-gui-comparison'/'charmm-gui-9263401822'

def source_bytes(path):
    with tarfile.open(ROOT/'charmm-gui-comparison'/'CBMAA_9263401822.tgz') as archive:
        return archive.extractfile(str(path.relative_to(JOB.parent)).replace('\\','/')).read()

def source_text(path):return source_bytes(path).decode('utf-8')

def graph(data):
    g=nx.Graph()
    for name,a in data['atoms'].items():g.add_node(name,element=name[0],type=a['type'])
    g.add_edges_from(data['bonds']);return g

def lookup(parameters,section,types):
    n=len(types);hits=[]
    for row in parameters.get(section,[]):
        if len(row)<=n:continue
        if any(all(x=='X' or x==y for x,y in zip(row[:n],t)) for t in (types,types[::-1])):
            hits.append(row)
    if not hits:return None
    score=max(sum(x!='X' for x in row[:n]) for row in hits)
    return sorted(set(tuple(float(v) for v in row[n:]) for row in hits if sum(x!='X' for x in row[:n])==score))

def main():
    original=json.loads((ROOT/'parsed'/'CB_amide_3_m1_monomer.json').read_text())
    rtf=source_text(JOB/'cbmaa/cbmaa.rtf')
    gui=parse(rtf)
    mapping=next(nx.algorithms.isomorphism.GraphMatcher(graph(original),graph(gui),
        node_match=lambda a,b:a==b).isomorphisms_iter())
    differences=[dict(source=a,gui=b,old=original['atoms'][a]['charge'],new=gui['atoms'][b]['charge'])
        for a,b in mapping.items() if original['atoms'][a]['charge']!=gui['atoms'][b]['charge']]
    # Check chemical bond order and formal charges before comparing parameter values.
    original_mol=Chem.MolFromMolFile(str(ROOT/'inputs/CB_amide_3_m1_monomer.sdf'),removeHs=False)
    gui_mol=Chem.MolFromMolBlock(source_text(JOB/'drawing_3D.mol'),removeHs=False)
    for m in (original_mol,gui_mol):Chem.RemoveStereochemistry(m)
    raw_graph_same=Chem.MolToSmiles(original_mol)==Chem.MolToSmiles(gui_mol)
    # Marvin writes carboxylate C-O resonance bonds as 'ar'. Normalize only
    # this explicit terminal O=C-O(-) motif; do not erase other bond orders.
    normalized_carboxylates=0
    for carbon in gui_mol.GetAtoms():
        oxy=[a for a in carbon.GetNeighbors() if a.GetSymbol()=='O' and a.GetDegree()==1]
        if carbon.GetSymbol()!='C' or len(oxy)!=2 or sum(a.GetFormalCharge() for a in oxy)!=-1:continue
        for oxygen in oxy:
            b=gui_mol.GetBondBetweenAtoms(carbon.GetIdx(),oxygen.GetIdx())
            b.SetBondType(Chem.BondType.SINGLE if oxygen.GetFormalCharge()==-1 else Chem.BondType.DOUBLE)
            b.SetIsAromatic(False);oxygen.SetIsAromatic(False)
        carbon.SetIsAromatic(False);normalized_carboxylates+=1
    Chem.SanitizeMol(gui_mol)
    chemical_graph_same=Chem.MolToSmiles(original_mol)==Chem.MolToSmiles(gui_mol)
    main_params=parse('read param\n'+source_text(JOB/'toppar/par_all36_cgenff.prm'))['parameters']
    custom=parse('read param\n'+source_text(JOB/'cbmaa/cbmaa.prm'))['parameters']
    g=graph(original)
    angles=[(a,b,c) for b in g for a,c in combinations(sorted(g[b]),2)]
    proper=set()
    for b,c in g.edges:
        for a in set(g[b])-{c}:
            for d in set(g[c])-{b,a}:proper.add(min((a,b,c,d),(d,c,b,a)))
    terms={'BONDS':original['bonds'],'ANGLES':angles,'DIHEDRALS':proper,'IMPROPERS':original['impropers']}
    stats={};mismatch=[]
    for section,rows in terms.items():
        for atoms in rows:
            types=tuple(original['atoms'][a]['type'] for a in atoms)
            old=lookup(original['parameters'],section,types)
            new=lookup(custom,section,types) or lookup(main_params,section,types)
            if old!=new:mismatch.append(dict(section=section,types=types,old=old,new=new))
        stats[section]=len(rows)
    report=dict(jobId='9263401822',molecule='CBMAA / spacer 3 / alpha methyl / vinyl monomer',
        programVersion='4.0',forceFieldVersion='5.0',versionHeadersMatch=all(x in rtf for x in ('program version 4.0','parameter files version 5.0')),
        rawBondRepresentationIdentical=raw_graph_same,carboxylateResonanceNormalizations=normalized_carboxylates,
        chemicalGraphIdenticalAfterResonanceNormalization=chemical_graph_same,atoms=len(mapping),atomMap=mapping,
        chargeDifferences=differences,maxChargeDifference=max([abs(x['old']-x['new']) for x in differences]+[0.]),
        checkedBondedTerms=stats,bondedDifferences=mismatch,
        allCheckedValuesMatch=chemical_graph_same and not differences and not mismatch,
        quotaObservation='One job succeeded. No daily/weekly job limit was displayed; absence of a quota is not proven.',
        sourceSHA256={str(p.relative_to(JOB)):hashlib.sha256(source_bytes(p)).hexdigest() for p in [JOB/'cbmaa/cbmaa.rtf',JOB/'cbmaa/cbmaa.prm',JOB/'toppar/par_all36_cgenff.prm']})
    path=ROOT/'charmm-gui-comparison'/'comparison.json';path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('atomMap','sourceSHA256')},indent=2))

if __name__=='__main__':main()
