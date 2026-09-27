# Ansys Workbench project

`wing_fea.wbpj` + `wing_fea_files/` is the final project, saved with **Ansys
2026 R1 Student**. Open the `.wbpj` in Workbench; keep the `.wbpj` and the
`_files` folder together and with the same base name.

| System | Content | Solver files |
|---|---|---|
| A | Original static system. It holds the shared Engineering Data; not solved | `dp0/global/MECH/SYS.mechdb` |
| B | Static Structural, aluminium, plus the Fatigue Tool | `dp0/SYS-2/MECH/` |
| C | Modal, aluminium (6 modes) | `dp0/SYS-3/MECH/` |
| D | Static Structural, "CFRP" (copy of B) | `dp0/SYS-4/MECH/` |
| E | Modal, "CFRP" | `dp0/SYS-5/MECH/` |

Each `MECH` folder keeps what the solver saw and produced:

* `ds.dat`: the complete MAPDL input deck (mesh, materials, contacts, loads)
* `solve.out`: solver log, including mass, frequencies and participation factors
* `file.rst`: results

The project was saved as `complete.wbpj` and renamed here together with its
`_files` folder. The renamed project opens in Workbench 2026 R1 Student
(checked in batch mode on a copy; the systems listed were "dummy", "static and
fatigue aluminium", "Modal aluminum", "CFRP static" and "CFRP modal"). The project file
also still contains stale absolute paths into the original project folder: the
geometry file `FEA Proeject (~recovered).step`, and a results file of the earlier
"unfailed attempt" project. The geometry is cached inside the project
(`dp0/SYS/DM/SYS.scdocx`). To update the geometry, re-point the geometry cell to
[`../cad/wing_final.step`](../cad/wing_final.step).

## Corrected project (`ansys/corrected/`)

`python validation/ansys_rerun/build_corrected_workbench.py` builds
`ansys/corrected/wing_fea_corrected.wbpj`, a copy of this project with the audit's
corrections, and re-solves all four systems. Result images are in
[`../results/corrected/`](../results/corrected/README.md).

The copy is about 110 MB; the script rebuilds it from scratch. If a system shows "Solve Required", re-solve just that system with
`--resolve al_static` (or `al_modal`, `cfrp_static`, `cfrp_modal`).

### Doing the same by hand in Workbench

1. Open `wing_fea.wbpj`, and double-click **Model** of a system to open Mechanical.
2. Right-click **Pressure** → **Suppress**.
3. Right-click the analysis (e.g. *Static Structural*) → **Insert → Commands**, and
   paste `validation/ansys_rerun/snippets/root_clamp.inp`. In static analyses, add a
   second **Commands** object with `design_load.inp`.
4. **Analysis Settings**: *Weak Springs* → **Off**. In modal analyses, set
   *Max Modes to Find* to **12** and add Total Deformation results for the extra
   modes.
5. **Solve**. Expected for aluminium:
   * design load: 0.81 mm and about 5.4 MPa;
   * modal: 22.8 Hz first bending and 221 Hz first torsion.

If MAPDL stops with "The memory (-m) size requested … is not currently
available", close other programs, or reduce the cores under **Home → Solve
Process Settings**. That is what happened on the first automated run.

The two Commands objects select nodes and elements by position, so they work on
any mesh of this geometry. The cleaner purely-GUI alternative:
* a **Fixed Support** scoped to the skin root edges and both spar root faces (or
  *Share Topology* in SpaceClaim, so the root cap is part of the skin);
* a **Pressure** on the upper skin, split into spanwise bands.

## Set-up as saved

| Item | Setting |
|---|---|
| Materials | "Aluminium aAlloy": E 70 GPa, ν 0.33, ρ 2700, UTS 310 MPa, S-N curve (10³…10⁷ cycles, 300…90 MPa). "cpfr": E 70 GPa, ν 0.3, ρ 1600, yield 500, UTS 600 MPa |
| Skin / caps | SHELL181, 8 mm / 0.1 mm, mid-surface offset |
| Ribs, spars | SOLID186 |
| Connections | 16 bonded contacts |
| Support | Fixed Support on the root cap (Body44) only |
| Load | Pressure by components, Z = 800 Pa, on all skin faces (resultant 434.6 N) |
| Analysis | linear, weak springs on; modal: 6 modes |
| Fatigue | stress-life, zero-based, Goodman, von Mises, Kf = 1 |

**Read [../docs/model_audit.md](../docs/model_audit.md) before using these
results**: the support and load definitions drive the reported peak stress,
deflection and frequencies.

## Inspecting the decks without Ansys

```bash
python tools/inspect_ds_dat.py ansys/wing_fea_files/dp0/SYS-2/MECH/ds.dat
python validation/read_ansys_results.py        # needs: pip install ansys-mapdl-reader
python validation/ansys_rerun/rerun_in_mapdl.py # re-solve in MAPDL, as saved and corrected (needs Ansys)
```
