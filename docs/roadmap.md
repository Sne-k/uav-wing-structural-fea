# Roadmap: what is left to do

The re-analysis in [`validation/`](../validation/README.md) closes the main gaps
found in the [audit](model_audit.md). This list covers what is still open,
most valuable first.

## Needs Ansys (Workbench / Mechanical / ACP)

1. **Rebuild the Workbench model with the corrections**, so the submitted
   screenshots match the verified numbers. The corrected set-up has already
   been re-solved in Ansys MAPDL on the saved mesh
   ([validation/ansys_rerun](../validation/ansys_rerun/README.md)). A corrected
   Workbench project has also been built and solved by script
   (`validation/ansys_rerun/build_corrected_workbench.py`; images in
   `results/corrected/`). What is left is the cleaner GUI version:
   * Fixed Support on the skin root edge and both spar root faces, or use
     *Share Topology* in SpaceClaim so the root cap is part of the skin.
   * Pressure *Normal To* the **upper skin only**, not "Components Z" on the
     whole skin. Split the upper surface into 4–6 spanwise bands carrying the
     lifting-line load in `validation/results/B_spanwise_lift_ultimate.csv`, or
     apply the same forces with an APDL command snippet.
   * Quad-dominant mesh on the skin, shell offset *Bottom* (inward), weak
     springs off.
   * Target values to reproduce: `validation/results/summary.json` (Part A:
     2.70 mm, 22.5 Hz).
2. **Composite layup in Ansys ACP** (T300/5208, the Part D layups) with ply-wise
   Tsai-Wu or Hashin, checked against Part D.
3. **Eigenvalue buckling** in Workbench (Static → Eigenvalue Buckling), checked
   against Part C.

## Engineering scope

4. **V-n diagram and load cases**: manoeuvre (+n, −n), gust (Pratt formula at
   V_C and V_D), and landing/handling cases. At the moment only +4 g limit /
   +6 g ultimate is analysed.
5. **Root fitting and spar-cap design**: the root joint (lugs, bolts or a
   carry-through spar) is where real wings fail. Bearing and shear-out checks.
6. **Stringers or a sandwich skin** if thinner skins are wanted: buckling
   governs the lighter aluminium options.
7. **Aerodynamic-structural coupling**: map XFLR5 or CFD pressures onto the
   structure and iterate on the deformed shape (loosely coupled, as in Chen et
   al.). Not significant at these deflections, but good to demonstrate.
8. **Flutter**: the typical-section estimate (Part F) is a first check. A
   modal-based p-k analysis with a doublet-lattice or strip-theory model of the
   whole wing would confirm the margin.
9. **Fatigue spectrum**: replace the constant-amplitude cycles with a mission
   spectrum (ground-air-ground plus gust exceedances), use Miner's rule, and
   use S-N data for the actual alloy and joint detail.

## Validation

10. **Physical test**: a simple static test (sandbag or whiffletree loading
    of a built wing, with tip deflection measured by a dial gauge) against the
    predicted deflection and strain. That turns the numbers into validated
    results.
