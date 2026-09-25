"""Fusion 360 script: import an airfoil .dat file (Selig format) as a fitted spline.

Install: UTILITIES > Add-Ins > Scripts and Add-Ins > + (My Scripts) > Create,
paste this file in, then run it and pick the .dat file (e.g. cad/naca2412.dat).

The coordinates are normalised (chord = 1), so scale the sketch to the chord you
need afterwards (0.25 m at the root in this project) and close the trailing edge
with a short line so the profile can be lofted.

Notes from building this project:
  * A header line or a UTF-8 BOM in the .dat file breaks float parsing. Only
    lines with two numbers are read here, and the BOM is stripped.
  * sketchFittedSplines.add() needs an adsk.core.ObjectCollection, not a list.
"""
import traceback

import adsk.core
import adsk.fusion


def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface

        dialog = ui.createFileDialog()
        dialog.title = "Select airfoil DAT file (Selig format)"
        dialog.filter = "*.dat"
        if dialog.showOpen() != adsk.core.DialogResults.DialogOK:
            return

        points = adsk.core.ObjectCollection.create()
        with open(dialog.filename, "r", encoding="utf-8-sig") as f:
            for line in f:
                values = line.split()
                if len(values) != 2:
                    continue
                try:
                    x, y = float(values[0]), float(values[1])
                except ValueError:
                    continue  # header line such as "NACA 2412"
                points.add(adsk.core.Point3D.create(x, y, 0))

        design = adsk.fusion.Design.cast(app.activeProduct)
        root = design.rootComponent
        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.sketchCurves.sketchFittedSplines.add(points)

        ui.messageBox(f"Airfoil imported: {points.count} points")
    except Exception:
        if ui:
            ui.messageBox("Failed:\n{}".format(traceback.format_exc()))
