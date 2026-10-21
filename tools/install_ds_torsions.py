"""Task-specific, evidence-supported DS torsions. Does not distribute DS sources."""
from pathlib import Path
from decimal import Decimal
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
DS=Path(sys.argv[1]) if len(sys.argv)>1 else Path('C:/Installment/Discovery Studio/share/forcefield/CHARMm')
tuples=[('CG2O1','CG301','CG321','CG311'),('CG301','CG321','CG311','CG2O1'),('NG2S1','CG321','CG321','CG324'),('CG324','CG321','CG321','CG2O3'),('CG2O1','CG301','CG321','CG321'),('CG301','CG321','CG321','CG2O1'),('CG301','CG321','CG321','CG2O2'),('CG2O2','CG311','CG321','CG311'),('CG311','CG321','CG311','CG2O2')]
mapping={'CG2O1':'C','CG2O2':'C','CG2O3':'C','CG301':'CT','CG311':'CT','CG321':'CT','CG324':'CT','NG2S1':'NP'}
sources=[ROOT/'exmaple/atomtypes.atp',*[DS/n for n in ['AMINOH.RTF','msicharmmh_templates.dat','CHARMM.eqv','PARM.PRM']]]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
rows=[]; inside=False
for line_no,line in enumerate((DS/'PARM.PRM').read_text().splitlines(),1):
    tokens=line.split('!')[0].split()
    if not tokens: continue
    if tokens[0]=='PHI': inside=True; continue
    if tokens[0]=='IMPHI': inside=False
    if inside and len(tokens)==7: rows.append((tokens,line_no,line.strip()))
audit=[]; additions=[]
for original in tuples:
    target=tuple(mapping[t] for t in original)
    matches=[]
    for tokens,line_no,line in rows:
        pat=tuple(tokens[:4])
        if any(all(p=='X' or p==t for p,t in zip(pat,seq)) for seq in (target,target[::-1])):
            matches.append((sum(p!='X' for p in pat),tokens,line_no,line))
    if not matches: raise RuntimeError(f'No proper DS torsion for {target}')
    specificity=max(m[0] for m in matches); selected=[m for m in matches if m[0]==specificity]
    patterns={min(tuple(m[1][:4]),tuple(m[1][:4][::-1])) for m in selected}
    if len(patterns)!=1: raise RuntimeError(f'Ambiguous equally specific patterns: {patterns}')
    for _,tokens,line_no,line in selected:
        k=Decimal(tokens[4])*Decimal('4.184'); phi=Decimal(tokens[6]); mult=int(tokens[5])
        additions.append(' '.join(original)+f' 9 {phi:.6f} {k:.6f} {mult} ; added by jf 20261021; PARM.PRM:{line_no}; DS {line}; inferred environment-specific mapping, not DS-engine assignment')
    audit.append(dict(original=original,ds_chemical_types=target,ds_lookup_types=target,equivalence_substitutions=[],match='exact' if specificity==4 else 'wildcard',source_rows=[dict(line=m[2],row=m[3]) for m in selected],status='supported for current acyclic molecular environments'))
ff=ROOT/'exmaple/ffbonded.itp'; data=ff.read_bytes(); marker=b'; BrushBuilder DS proper torsion supplement'
existing={tuple(line.split(';')[0].split()[:5]) for line in data.decode().splitlines()}
new=[row for row in additions if tuple(row.split()[:5]) not in existing and tuple(row.split()[:4][::-1]+['9']) not in existing]
if new:
    with ff.open('ab') as f: f.write(('\n'+marker.decode()+'\n[ dihedraltypes ]\n'+'\n'.join(new)+'\n').encode())
report=dict(source_hashes=hashes,mapping_evidence=dict(ATP_lines={'CG2O1':171,'CG2O3':173,'CG301':191,'CG311':193,'CG321':196,'CG324':199,'NG2S1':264},DS_MASS_lines={'C':26,'CT':52,'NP':79},DS_template_lines={'C':104,'CT':78,'NP':371}),environment='C1-C3 and backbone/sidechain carbons are acyclic tetrahedral explicit-H carbon; CG324 is CH2 next to quaternary N; CG2O1 is amide carbonyl, CG2O3 carboxylate carbon; NG2S1 is neutral secondary amide N. No ring or aromatic alternatives apply. Source CBAA internal CG311 has no H despite its ATP label: preserved source type, DS CT still matches its actual tetrahedral carbon environment.',entries=audit,limitations='Parameter transfer is inferred from definitions/templates, not running DS; preprocessing/EM do not establish cross-force-field physical equivalence.')
report['mapping_evidence']['ATP_lines']['CG2O2']=172
report['environment']+=' CG2O2 is the ester carbonyl carbon (ATP:172), mapping to DS carbonyl C by the same MASS/template evidence.'
(ROOT/'parameter_audit').mkdir(exist_ok=True)
(ROOT/'parameter_audit/ds_torsions.json').write_text(json.dumps(report,indent=2))
print('\n'.join(additions))
