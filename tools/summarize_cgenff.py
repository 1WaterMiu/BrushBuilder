"""Produce a reproducible inventory without conflating fragments and surfaces."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]/'cgenff'

def main():
    manifest=json.loads((ROOT/'manifest.json').read_text())
    rows=json.loads((ROOT/'audit.json').read_text())+json.loads((ROOT/'charmm-gui-audit.json').read_text())
    records={r['file']:r for r in rows}
    assert len(records)==len(rows),'Duplicate models in audit'
    expected={Path(r['file']).stem for r in manifest}
    assert set(records)<=expected
    assert all(r['bondedCoverage']['complete'] and r['atomMapExact'] and r['bondGraphExact'] for r in rows)
    assert all(abs(r['charge'])<1e-5 for r in rows)
    kinds={k:sum(s.endswith('_'+k) for s in records) for k in ('monomer','pentamer','anchor5')}
    families={f:{k:sum(s.startswith(f+'_') and s.endswith('_'+k) for s in records)
                 for k in kinds} for f in ('CB','SB','PB','sulfabetaine','PC','TMAO')}
    high=[{k:r[k] for k in ('file','chargePenalty','parameterPenalty')} for r in rows
          if max(r['chargePenalty'],r['parameterPenalty'])>50]
    summary=dict(uniqueModels=len(records),counts=kinds,families=families,
        originalCGenFFModels=len(rows)-len(json.loads((ROOT/'charmm-gui-audit.json').read_text())),
        newCHARMMGUIModels=len(json.loads((ROOT/'charmm-gui-audit.json').read_text())),
        missingSubmittedModelInputs=sorted(expected-set(records)),
        allAuditedGraphsAndBondedCoveragePass=True,highPenaltyModels=high,
        missingAnchorContexts=sorted(s[:-8] for s in records if s.endswith('_monomer') and s[:-8]+'_anchor5' not in records),
        mdReadySurfaces=0,
        remainingValidation=['71 chemistry combinations need anchored-context results',
            'N=1..4 need separate anchored short-chain contexts',
            'Assembled-surface force-field compatibility, GROMACS and fixed-silica minimization checks remain',
            'High CGenFF penalties require scientific validation; coverage alone does not establish accuracy'],
        quotaObservation='110 new models plus the comparison and one duplicate succeeded; no quota warning observed. Unlimited use is not established.')
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('highPenaltyModels','missingAnchorContexts')},indent=2))
    print('High-penalty models:',len(high))

if __name__=='__main__':main()
