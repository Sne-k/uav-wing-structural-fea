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
`_files` folder. Workbench normally accepts that, but it has **not been tested**,
because Ansys was not installed when this repository was made. The project file
also still contains stale absolute paths into the original project folder: the
geometry file `FEA Proeject (~recovered).step`, and a results file of the earlier
"unfailed attempt" project. The geometry is cached inside the project
(`dp0/SYS/DM/SYS.scdocx`). To update the geometry, re-point the geometry cell to
[`../cad/wing_final.step`](../cad/wing_final.step).

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
```
