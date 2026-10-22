# BrushBuilder

A local web application for fixed-silica zwitterionic brush construction and real GROMACS minimization.

## Current CGenFF migration status — 2026-10-21

**This section supersedes the legacy parameter and MD descriptions below.** Website builds no longer use supplied template charges for organic atoms, except the explicitly approved Si–C interface (anchor C1/H4/H5). Silica charges and all fixed silica coordinates are preserved. No charge renormalization or manual alpha-H charge patch is applied.

All 72 monomers and their 72 saturated pentamers have explicit-H MOL2/SDF inputs, atom maps and graph round-trip checks in `cgenff/inputs`. Actual CGenFF 5.0 web results are available for **17 monomers + 17 pentamers + 1 anchored pentamer = 35 molecules**. These cover all 16 CB combinations and SB/ester/spacer 1/no methyl. The anchored result is CB/amide/spacer 3/alpha methyl. All 35 output atom/bond graphs match their inputs and all enumerated bonds, angles, proper torsions and declared impropers have parameters in the downloaded streams. Penalties above 50 are retained in `cgenff/audit.json`; parameter completeness does not establish quantitative accuracy.

The original CGenFF website reported its free-job limit. The user subsequently authorized testing CHARMM-GUI and continuing there if equivalent. The CBMAA spacer-3 monomer comparison passed: program 4.0, force field 5.0, all 42 atom types and charges identical, and identical effective parameters for 41 bonds, 75 angles, 102 proper torsions and 2 impropers. CHARMM-GUI's custom parameters were combined with its matching base force field for this comparison; raw file lengths are not a valid comparison. The only chemical representation normalization was carboxylate resonance. See `cgenff/charmm-gui-comparison/comparison.json`.

Additional CHARMM-GUI results and their current counts are recorded in `cgenff/charmm-gui-audit.json`. No quota warning has appeared during these submissions; this does not establish unlimited usage. `tools/import_charmm_gui.py` verifies the ordered explicit-H chemical graph, independent RTF connectivity, charge sum, version, and bonded coverage before importing renamed atoms into the library. Original full TGZ downloads are retained locally; compact ZIPs preserve the scientific source files and the matching base parameter file with provenance hashes. New organic charges are never filled from supplied templates. `tools/prepare_cgenff.py` regenerates inputs, and `tools/audit_cgenff.py` audits the original CGenFF-site results separately.

For CB/amide/spacer 3/alpha methyl with N >= 5, all organic charges are mapped from the anchored pentamer, distinguishing first, interior and last units; chain charge remains zero. Other combinations retain only available interior/end assignments and mark unresolved atoms `UNASSIGNED`. Short chains need separate anchored contexts. Available CGenFF results supply ordered impropers, including CB carboxylates missing from the legacy construction. Structure ZIPs include available streams and atom maps.

MD submission/export remains blocked until assembled-surface bonded/nonbonded compatibility and GROMACS checks are complete. Legacy MD jobs using template charges cannot be downloaded as new CGenFF results. The matching official CHARMM36/CGenFF 5.0 archive is retained for the next compatibility check; it has not replaced the supplied force field. Earlier GROMACS tests below validate only the legacy parameter set.

### Completed CHARMM-GUI batch

All 110 remaining monomer/pentamer models have now been retrieved and audited. The library contains **72 monomers + 72 pentamers + 1 anchored pentamer = 145 unique models**. All have exact input/output graph checks and complete bonded coverage. The comparison job and one duplicate TMAO pentamer are excluded from this total. See `cgenff/summary.json` for per-family counts, the 25 high-penalty models, and outstanding validation. The 71 missing anchored contexts and the short-chain contexts are separate from completion of this 144-model monomer/pentamer batch. No assembled surface is yet marked MD-ready under the new parameter set.

Run `python tools/compare_charmm_gui.py` to reproduce the comparison from its retained original archive and `python tools/summarize_cgenff.py` to check the combined inventory.

## Start

Run in C:\personal\BrushBuilder:

    python -m pip install -r requirements.txt
    python server.py

Open http://127.0.0.1:8765. The GROMACS runner detects native gmx_d/gmx or the installed Ubuntu WSL executable. Double precision is preferred because frozen substrate interactions can create a large constant energy contribution.

## Head groups and names

Head group, linker, charge spacer and alpha methyl are independent selectors. CB, SB, phosphate-type PB and sulfabetaine support 1–4 methylenes between charged groups. PC has fixed O–CH2–CH2–N connectivity; TMAO-derived heads have a direct N–O bond. There are 36 head/linker/spacer combinations, or 72 including methyl states. Ester uses an O–(CH2)2 connector; amide uses NH–(CH2)3. These connectors are distinct from the adjustable charge spacer.

Naming is derived from the actual graph: SBA = acrylate, SBMA = methacrylate, SBAA = acrylamide, SBMAA = methacrylamide. CB follows CBA/CBMA/CBAA/CBMAA. The supplied file named SUP_CBAA3N20.itp contains alpha methyl, so its structure is displayed as **CBMAA · spacer 3**, retaining the original filename as provenance. No atom is removed to make the structure agree with the filename. Numerical spacer labels explicitly count CH2 groups and are not assumed to reproduce all naming conventions used in papers.

Sources and scope are embedded in chemistry.py and the UI. [RSC Scheme 1](https://pubs.rsc.org/en/content/articlepdf/2026/lp/d6lp00019c) distinguishes the four SB monomers. [CB synthesis definitions](https://patents.google.com/patent/US20140275614A1/en) distinguish CBAA and CBMAA. [MPC](https://pubchem.ncbi.nlm.nih.gov/compound/128934) is graph-matched to its published SMILES. [Sulfabetaine](https://pmc.ncbi.nlm.nih.gov/articles/PMC6631746/) uses sulfate ester, not sulfonate. [TMAO-derived acrylamide](https://pmc.ncbi.nlm.nih.gov/articles/PMC6570511/) uses a propyl connector and direct N–O bond. [Phosphate-type phosphobetaine](https://pubmed.ncbi.nlm.nih.gov/1880723/) supports the user-confirmed N+–(CH2)n–O–P(=O)(O−)–OH head motif. Enumerated derivatives are not claimed to be separately synthesized/literature-verified monomers; PB is this explicit phosphate variant, not an ambiguous universal PB structure.

New combinations contain **no assigned partial charges or force-field types**. Their preview and structure ZIP are available; MD submission/export is blocked. Structure ZIP contains full explicit-H vinyl monomer SDF, capped repeat SDF, surface GRO and connectivity JSON. Formal charges in SDF represent chemical valence and are not MD partial charges. Existing reference parameters and the previously implemented alpha-H patch remain unchanged. Structure and MD downloads are separate.

N20 template checks compare all named atoms, every bond and explicit-H repeat graphs: source CBAA 915 atoms/916 bonds; SBMA 855 atoms/856 bonds; zero missing or extra atoms/bonds. Atom counts include the silica/anchor portion of each grafted molecule. Coordinates are independently generated and therefore are not expected to match source polymer coordinates.

## Construction controls

- CBAA3 / amide and SBMA / ester source chemical definitions.
- N = 1–60 repeating units, including correct chain termination.
- Alpha methyl checkbox: checked retains source methyl. Unchecked substitutes alpha H with +0.09 charge, compensates the bonded carbon by −0.09, and selects CG311/CG321 plus HGA1/HGA2 by terminal valence. This explicit charge-conserving patch is not a new QM fit.
- High: 20 chains, approximately 0.2938 chains/nm².
- Low: exactly remove original graft sites 2, 7, 12, 17, 4, 9, 14, 19; retain 12 chains (60% of high), approximately 0.1763 chains/nm². The removed sites become original SUB units with O34/H4 and SIB Si1 charge +1.1.
- Numeric site-selection UI removed.

## Independent geometry

Organic coordinates are generated from chemical graphs, not copied from either source GRO. Deterministic RDKit ETKDG/UFF fragment generation produces fixed polymerization ports and independent anchor / amide / ester-sidechain geometry. The persistent JSON templates in geometry_templates/ are reused for every length and graft site.

UFF is only an initial fragment-coordinate tool. Surface energy minimization uses the supplied CHARMM-family force field in GROMACS.

Connectivity, quantities and source charges of the existing Si–(CH2)3–anchor remain unchanged. Repeat sidechain orientation alternates by reflection to avoid identical neighboring bulky heads overlapping; this deterministic sequence is copied to every chain. No claim of equilibrated brush structure is made for the initial coordinates. EM naturally changes polymer coordinates afterward.

The original silica coordinates remain strictly unchanged. Restored OH geometry is built from bare SUB geometry with the canonical +z graft normal. All silica atoms, including OH, are frozen in XYZ during minimization and MD. Only box height may grow to accommodate longer chains; the surface area is fixed.

## Topology and parameters

Bond connectivity is assembled from the supplied ITPs. Angles, proper torsions and 1–4 pairs are enumerated from the resulting graph; anchor and repeat carbonyl improper terms are generated for every unit. Chain-end H and charges are regenerated. CBAA's original internal CG311 type is preserved even though its original methylated connectivity has four heavy neighbors; this source inconsistency is recorded in the parameter audit.

Missing proper torsions are supplied through the ds2gmx workflow, with environment-specific mapping to DS C / CT / NP and lookup in PARM.PRM. All nine supplemented tuples currently select the X CT CT X proper wildcard, converted as k = 0.15 × 4.184 = 0.6276 kJ/mol, phase 0°, multiplicity 3. They are appended to exmaple/ffbonded.itp with "added by jf 20261021", source line and original row. parameter_audit/ds_torsions.json records hashes, mapping evidence and limitations. No proprietary DS source files are bundled.

## Minimize and export

1. Change structure, N, density or alpha methyl; preview and topology counts update.
2. Click Generate and minimize. The background job runs grompp without maxwarn and steepest descent with Fmax < 500 kJ/mol/nm.
3. Every frozen coordinate is compared before/after using GRO and GROMACS dump of TPR/TRR. Failed convergence or any detected silica displacement blocks the MD-ready download.
4. Download the complete ZIP after success. It contains optimized and initial GRO, generated BRUSH.itp, topol.top, all included force-field files, TIP3P/ions, index, EM/NVT MDPs, logs, TPR, processed topology, coordinate templates, parameter provenance and manifest.
5. In the extracted directory run bash run_md.sh with GROMACS and Python 3. Prefer GMX=/path/to/gmx_d bash run_md.sh. This creates a TIP3P water system, performs solvated EM, then runs a 100 ps NVT example. All stages preserve frozen silica and verify its coordinates. This is pure water with no added protein/salt; scientific conditions and production duration should be configured for the study.

Only validated completed jobs are downloadable as MD-ready. Jobs and diagnostics are stored in output/jobs/ (Git ignored). Closing the page does not stop a submitted job. Completed jobs can be recovered after server restart; interrupted jobs must be resubmitted.

NVT uses a separate Silica temperature group with zero degrees of freedom and disabled thermostat, an explicit 0.1 nm neighbor-list buffer, and no automatic bond constraints. Applying h-bond constraints can alter frozen silica H during constraint initialization; rigid water still uses its topology SETTLE definition.

## Verification

    python -m unittest -v
    python tools/validate_gromacs.py
    python tools/validate_gromacs.py --solvated

Unit tests cover silica invariance, exact low-density site removal, charge conservation, terminal valence, graph indices, identical translated chain copies, independence from source organic GRO coordinates, and ZIP dependency closure. GROMACS boundary checks cover both chemistries, methyl states and N=1/60. The solvated test checks actual EM and a short NVT smoke run, not production sampling.

The local development server binds loopback only. Source data are retained under the user-provided directory name exmaple/.
