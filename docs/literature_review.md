# Literature review

These five sources shaped the project. The PDFs are **not** redistributed in
this repository because they are copyrighted; follow the links or DOIs instead.

| # | Reference | Access |
|---|---|---|
| 1 | Tang J., Xi P., Zhang B., Hu B., *A finite element parametric modeling technique of aircraft wing structures*, Chinese Journal of Aeronautics 26(5), 1202–1210, 2013 | [doi:10.1016/j.cja.2013.07.019](https://doi.org/10.1016/j.cja.2013.07.019) (open access) |
| 2 | Chen H. H., Chang K. C., Tzong T., Cebeci T., *Aeroelastic Analysis of Aircraft: Wing and Wing/Fuselage Configurations*, McDonnell Douglas report MDC 97K0164 for NASA Ames (contract NAS2-14091), 1997 | [NTRS 19980003841](https://ntrs.nasa.gov/citations/19980003841) |
| 3 | Noor Abdi A., Kapesa T. G., *Optimization of Composite Structures of an Aircraft Wing*, Int. J. of Modern Research in Engineering and Technology 8(7), 2023 | ijmret.org, ISSN 2456-5628 |
| 4 | Peruru S. P., Abbisetti S. B., *Design and Finite Element Analysis of Aircraft Wing Using Ribs and Spars*, IRJET 4(6), 2133–2139, 2017 | irjet.net, e-ISSN 2395-0056 |
| 5 | Vinay Kumar P., Rahul Raj I., Snehith Reddy M., Siva Prasad N., *Design and Finite Element Analysis of Aircraft Wing using Ribs and Spars*, Turkish J. of Computer and Mathematics Education 12(8), 3224–3230, 2021 | turcomat.org |

## What each source contributes

**1. Parametric FE modelling (Tang et al., 2013).** Building an FE model of a
wing by hand (geometry, mid-surfaces, mesh, properties) is slow and error-prone.
The authors encode design rules and expert knowledge in templates. A *skeleton
model* (the layout of spars, ribs and stringers) drives a *geometric mesh model*,
which drives the *FE model*, so the whole model regenerates when a parameter
changes. It targets preliminary design, where many layouts must be analysed
quickly.
*Relevance:* the project's CAD/FE hand-off was the most time-consuming part.
Parametrising the layout (rib pitch, spar positions, thicknesses) is exactly how
[`validation/wingfe/wing.py`](../validation/wingfe/wing.py) rebuilds the model.

**2. Aeroelastic analysis with a loosely coupled interface (Chen et al., 1997).**
Aerodynamic pressures from a flow solver are converted into structural nodal
forces by virtual work. The FE deflections are then mapped back onto the
aerodynamic surface, and the loop repeats until loads and shape converge. For a
transport wing at cruise, the elastic deformation noticeably changes the pressure
distribution and shock position compared with a rigid wing, and the coupled
results agree better with wind-tunnel data.
*Relevance:* it motivates deriving the structural load from an aerodynamic
analysis rather than assuming a pressure. It also shows that modal and aeroelastic
checks, not just static stress, belong in wing design.

**3. Composite wing optimisation (Noor Abdi & Kapesa, 2023).** A broad overview
of optimising composite wing structures together with aerodynamic shape
(twist and thickness distributions, pressure distributions) for weight and
performance. It gives limited quantitative detail about the structural model.
*Relevance:* background for the aluminium-versus-composite comparison and for
treating weight as an objective.

**4. Wing with ribs and spars in Pro/E and Ansys (Peruru & Abbisetti, 2017).** A
trainer-aircraft wing with skin, 15 ribs, an I-section front spar and a
C-section rear spar, analysed in Ansys 14.5 (static, fatigue, modal). It compares
Al 6061-T8, S2-glass and carbon/epoxy at several speeds.
*Relevance:* this is the template for the project workflow (CAD → Workbench →
static/fatigue/modal → material comparison). The project's CFRP properties
(E = 70 GPa, ν = 0.3, ρ = 1.6 g/cc, isotropic) come from this paper's material
table. Treating carbon/epoxy as isotropic with the same modulus as aluminium is
the simplification that made the project's CFRP comparison a pure density change.

**5. Wing with ribs, spars and winglets (Vinay Kumar et al., 2021).** A CATIA
wing (15 ribs, I- and C-section spars) with winglets at 25° and 45°. It is
analysed in Ansys for S2-glass, Kevlar-49 and boron fibre. Boron fibre gave the
lowest deformation and stress of the three materials. The winglet angle made
little difference: in the paper's own table, 25° gives the lowest deformation
(0.000444 vs 0.000466) and 45° a marginally lower stress (12.236 vs 12.254,
0.15 % lower).
*Relevance:* same workflow as paper 4. It also shows how a configuration change
(winglets) can be studied within the same FE framework.

## Research gaps the project set out to address

1. Studies tend to treat aerodynamics and structure separately. Few derive the
   structural load from an aerodynamic analysis of the same wing (papers 2, 4, 5).
2. Many student-level studies stop at static stress; modal and fatigue checks are
   often missing (papers 4, 5).
3. Composite comparisons usually use isotropic "equivalent" properties rather
   than laminate (ply-level) modelling (papers 4, 5).
4. Building the FE model by hand is slow and hard to iterate (paper 1).

## Objectives (as stated for the project)

1. Design a 3D wing with skin, spars and ribs in CAD.
2. Build an FE model in Ansys Workbench.
3. Run static structural analysis for stress and deformation.
4. Run modal analysis for natural frequencies and mode shapes.
5. Run fatigue analysis to estimate life.
6. Compare aluminium with CFRP and evaluate strength-to-weight.

How far the original analysis met these objectives, and how the re-analysis in
[`validation/`](../validation/README.md) closes the remaining gaps, is covered in
[model_audit.md](model_audit.md).
