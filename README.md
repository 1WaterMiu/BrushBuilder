# BrushBuilder

A local web application for fixed-silica zwitterionic brush construction and real GROMACS minimization.

## Start

Run in C:\personal\BrushBuilder:

    python -m pip install -r requirements.txt
    python server.py

Open http://127.0.0.1:8765. The GROMACS runner detects native gmx_d/gmx or the installed Ubuntu WSL executable. Double precision is preferred because frozen substrate interactions can create a large constant energy contribution.

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
