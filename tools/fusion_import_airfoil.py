"""Fusion 360 script: import an airfoil .dat file (Selig format) as a fitted spline.

Install: UTILITIES > Add-Ins > Scripts and Add-Ins > + (My Scripts) > Create,
paste this file in, then run it and pick the .dat file (e.g. cad/naca2412.dat).

The .dat coordinates are normalised (chord = 1). Fusion's scripting API works in
centimetres, so unscaled points would give a 1 cm chord (to reach 0.25 m you
would have to scale by 25, not 0.25). The script therefore scales the points to
CHORD_M directly. Afterwards, close the trailing edge with a short line so the
profile can be lofted.

Notes from building this project:
  * A header line or a UTF-8 BOM in the .dat file breaks float parsing. Only
    lines with two numbers are read here, and the BOM is stripped.
  * sketchFittedSplines.add() needs an adsk.core.ObjectCollection, not a list.
"""
import traceback

import adsk.core
import adsk.fusion

CHORD_M = 0.25          # chord of the imported profile [m] (root chord of this project)
API_UNITS_PER_M = 100.0  # the Fusion API uses centimetres


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
                s = CHORD_M * API_UNITS_PER_M
                points.add(adsk.core.Point3D.create(x * s, y * s, 0))

        design = adsk.fusion.Design.cast(app.activeProduct)
        root = design.rootComponent
        sketch = root.sketches.add(root.xYConstructionPlane)
        sketch.sketchCurves.sketchFittedSplines.add(points)

        ui.messageBox(f"Airfoil imported: {points.count} points")
    except Exception:
        if ui:
            ui.messageBox("Failed:\n{}".format(traceback.format_exc()))
