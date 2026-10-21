# BrushBuilder

Local molecular-surface studio built from the supplied `exmaple/` files. No Python packages or frontend build step required.

```powershell
cd C:\personal\BrushBuilder
python server.py
```

Open http://127.0.0.1:8765. Run `python -m unittest -v` for coordinate, charge, replacement and export checks.

## Supported construction

- Original CBAA3 and SBMA N=20 dry surfaces; high density has 20 chains, low density 10.
- High-density coordinates and polymer topologies remain unchanged. Original water, ions and protein are excluded from the new dry-surface export.
- Low density uses deterministic farthest-point selection on the periodic xy plane. The 20 starting-site choices need not all give distinct subsets.
- Removed graft molecules are replaced by the supplied SUB topology, preserving every existing silica coordinate. Missing O34/H4 geometry is transferred from an original bare SUB and rotated to the former Si1–C1 direction. Si1 changes from SSU +0.825 to SIB +1.1; restored O34 is OSU −0.692 and H4 is HSU +0.417.
- All SUB/SUP silica atoms, including restored OH, are in the `Silica` freeze group. `freeze.mdp` is an integration snippet, not a complete simulation protocol. Keep these atoms fixed in every simulation stage; inspect added OH contacts and minimize mobile atoms. Do not allow box scaling to move the fixed substrate.
- Both supplied templates contain a backbone methyl substituent. Both use an amide-containing surface anchor; the repeat-unit sidechain is amide for CBAA3 and ester for SBMA. These distinct connections are not interchangeable.

## Validation boundaries

The interface shows real atomic coordinates (heavy-atom projection) and source atom charges/bonds. Top view labels actual graft sites. There is no synthetic decorative molecule or mock export.

Only fully defined N=20 chemistry templates are enabled. Other lengths, methyl removal, new heads, or exchanging sidechain linkers need explicitly defined repeat/terminal patches and validation; the application does not infer charges for unseen chemistry.

The supplied dry-surface force-field includes are complete: `cmap.itp`, `nbfix.itp`, and `ffmissingdihedrals.itp` are present and included in exports. `ions.itp` is also supplied. The original solvated-system topologies still require `tip3p.itp`, and `SBMA.top` additionally requires `PROA.itp`; neither is needed for the current dry-surface export. `CBAA.top` retains its original `charmm36-jul2022.ff/` include paths; use the generated `topol.top` for the flat export layout.

GROMACS is not installed in the current development environment. Generated models are **not GROMACS validated**. Complete force-field dependencies, inspect restored OH geometry, solvate, update groups/topology and run preprocessing/minimization before production. Original files are never modified. ZIP manifests record source SHA256 hashes, retained sites, actual density, charge and validation state.

## Layout

- `builder.py`: parsers, template audit, site selection, OH restoration, ZIP exporter.
- `server.py`: loopback-only standard-library HTTP server.
- `web/`: responsive UI, actual coordinate projection, parameter inspector.
- `test_builder.py`: scientific integrity regression tests.

This is a local research tool, not a public production server. The explicit GitHub remote is `https://github.com/1WaterMiu/BrushBuilder`.
