"""Build a corrected copy of the Workbench project and re-solve it (Ansys 2026 R1).

1. Writes two mesh-independent APDL command snippets to validation/ansys_rerun/snippets/:
     root_clamp.inp    clamp every node in the root plane (y = 0): skin edge, spar
                       root faces and root cap
     design_load.inp   design ultimate lift (5 kg UAV, 6 g: 147.15 N per semi-span),
                       lifting-line spanwise shape, centre of pressure at 25 % chord,
                       applied as vertical nodal forces on the upper-skin elements
   The same files can be pasted into Mechanical by hand (Insert > Commands).
2. Copies ansys/wing_fea.wbpj to ansys/corrected/wing_fea_corrected.wbpj and, in every
   system: suppresses the old 800 Pa pressure, adds the snippets, switches weak
   springs off, asks for 12 modes in the modal systems, re-solves, and exports result
   images and values (via a Mechanical script run in Workbench batch mode).
3. Collects the images and values in results/corrected/.

Run:  python validation/ansys_rerun/build_corrected_workbench.py   (~10-20 min)
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "validation"))
from wingfe.aero import lifting_line  # noqa: E402
from wingfe.wing import SPAN, chord  # noqa: E402

AWP = os.environ.get("AWP_ROOT261", r"D:\apps\ANSYS\ANSYS Inc\ANSYS Student\v261")
RUNWB2 = os.path.join(AWP, "Framework", "bin", "Win64", "RunWB2.exe")
SNIP = os.path.join(HERE, "snippets")
L_ULT = 6.0 * 5.0 * 9.81 / 2
SYSTEMS = (("SYS 2", "al_static"), ("SYS 1", "al_modal"), ("SYS 3", "cfrp_static"), ("SYS 4", "cfrp_modal"))

CLAMP = """! Correction 1: clamp the whole root section.
! Every node in the root plane y = 0 (skin root edge, spar root faces, root cap) is fixed.
! Units: m (the solver input is written in SI).
nsel,s,loc,y,-1e-5,1e-5
d,all,all
allsel,all
"""


def design_load_snippet():
    ys, shape, _, _ = lifting_line(chord, SPAN)
    grid = np.linspace(0.0, SPAN, 41)
    vals = np.interp(grid, ys, shape)
    arr = "\n".join(f"ysh({i + 1})={v:.6f}" for i, v in enumerate(vals))
    c0, c1 = chord(0.0), chord(SPAN)
    return f"""! Correction 2: design ultimate load on the upper skin.
! Lift for a 5 kg UAV at n = 6 (4 g limit x 1.5): {L_ULT:.2f} N per semi-span.
! Spanwise shape from lifting-line theory (tabulated below, per metre, integral = 1),
! chordwise dp ~ 3 (1 - x/c)^2 (centre of pressure at 25 % chord).
! Applied as vertical nodal forces on the skin elements above the camber line.
! Mesh-independent: it loops over the skin elements of whatever mesh is solved.
! Replaces the 800 Pa pressure (suppress that object). Units: m, N.
ltot={L_ULT:.4f}
*dim,ysh,array,41
{arr}
esel,s,ename,,181
esel,u,cent,y,-1e-5,1e-5
esel,u,cent,y,{SPAN}-1e-5,{SPAN}+1e-5
*get,ne,elem,0,count
raw=0
*do,pass,1,2
*if,pass,eq,2,then
scl=ltot/raw
fcum,add
*endif
el=0
*do,k,1,ne
el=elnext(el)
yc=centry(el)
c={c0}+({c1}-{c0})*yc/{SPAN}
xi=centrx(el)/c
*if,xi,lt,0,then
xi=0
*endif
*if,xi,gt,1,then
xi=1
*endif
*if,xi,lt,0.4,then
zcam=0.02/0.16*(0.8*xi-xi*xi)*c
*else
zcam=0.02/0.36*(0.2+0.8*xi-xi*xi)*c
*endif
*if,centrz(el),gt,zcam,then
n1=nelem(el,1)
n2=nelem(el,2)
n3=nelem(el,3)
n4=nelem(el,4)
*if,n4,eq,n3,then
az=0.5*abs((nx(n2)-nx(n1))*(ny(n3)-ny(n1))-(nx(n3)-nx(n1))*(ny(n2)-ny(n1)))
nn=3
*else
az=0.5*abs((nx(n3)-nx(n1))*(ny(n4)-ny(n2))-(nx(n4)-nx(n2))*(ny(n3)-ny(n1)))
nn=4
*endif
t=yc/{SPAN / 40}
i0=nint(t-0.5)+1
*if,i0,lt,1,then
i0=1
*endif
*if,i0,gt,40,then
i0=40
*endif
fr=t-(i0-1)
lsh=ysh(i0)*(1-fr)+ysh(i0+1)*fr
fe=lsh/c*3*(1-xi)**2*az
*if,pass,eq,1,then
raw=raw+fe
*else
fn=fe*scl/nn
f,n1,fz,fn
f,n2,fz,fn
f,n3,fz,fn
*if,nn,eq,4,then
f,n4,fz,fn
*endif
*endif
*endif
*enddo
*enddo
fcum,repl
allsel,all
"""


MECH_SCRIPT = r'''
import os, traceback
log = open(os.path.join(OUT_DIR, "mech_log_%s.txt" % SYS_TAG), "w")
vals = open(os.path.join(OUT_DIR, "values_%s.txt" % SYS_TAG), "w")

def export(obj, name):
    try:
        obj.Activate()
        Graphics.Camera.SetSpecificViewOrientation(ViewOrientationType.Iso)
        Graphics.Camera.SetFit()
        s = Ansys.Mechanical.Graphics.GraphicsImageExportSettings()
        s.Resolution = GraphicsResolutionType.EnhancedResolution
        s.Background = GraphicsBackgroundType.White
        s.Width = 1600
        s.Height = 900
        s.CurrentGraphicsDisplay = False
        Graphics.ExportImage(os.path.join(OUT_DIR, "%s_%s.png" % (SYS_TAG, name)), GraphicsImageExportFormat.PNG, s)
        log.write("image ok %s\n" % name)
    except Exception:
        log.write("image FAILED %s\n%s\n" % (name, traceback.format_exc()))

def apply_corrections(a, static):
    for c in a.Children:
        if c.GetType().Name == "Pressure":
            c.Suppressed = True
            log.write("suppressed %s\n" % c.Name)
    try:
        a.AnalysisSettings.WeakSprings = WeakSpringsType.Off
    except Exception:
        log.write("weak springs setting not available\n")
    cs = a.AddCommandSnippet()
    cs.Name = "Correction 1 - clamp whole root section"
    cs.Input = open(os.path.join(SNIPPET_DIR, "root_clamp.inp")).read()
    if static:
        cs2 = a.AddCommandSnippet()
        cs2.Name = "Correction 2 - design load on upper skin"
        cs2.Input = open(os.path.join(SNIPPET_DIR, "design_load.inp")).read()
    else:
        a.AnalysisSettings.MaximumModesToFind = 12
        have = set()
        for r in a.Solution.Children:
            if r.GetType().Name == "TotalDeformation":
                have.add(int(r.Mode))
        for m in range(1, 13):
            if m not in have:
                r = a.Solution.AddTotalDeformation()
                r.Mode = m

def solve_with_retry(a):
    a.Solve(True)
    if str(a.Solution.Status) != "Done":
        # e.g. MAPDL could not reserve its memory: retry with fewer, shared-memory cores
        log.write("first solve ended in state %s; retrying on 2 cores\n" % a.Solution.Status)
        try:
            cfg = ExtAPI.Application.SolveConfigurations["My Computer"]
            cfg.SolveProcessSettings.MaxNumberOfCores = 2
            cfg.SolveProcessSettings.DistributeSolution = False
        except Exception:
            log.write("could not change the solve settings\n%s\n" % traceback.format_exc())
        a.Solve(True)
    log.write("solved, state %s\n" % a.Solution.Status)

try:
    a = ExtAPI.DataModel.Project.Model.Analyses[0]
    static = str(a.AnalysisType) == "Static"
    if not RESOLVE_ONLY:
        apply_corrections(a, static)
    solve_with_retry(a)
    results = []
    for r in a.Solution.Children:
        results.append(r)
        if r.GetType().Name == "FatigueTool":
            results.extend(list(r.Children))
    for r in results:
        tn = r.GetType().Name
        if tn in ("SolutionInformation", "FatigueTool"):
            continue
        line = "%s | %s" % (r.Name, tn)
        for prop in ("Maximum", "Minimum", "ReportedFrequency"):
            try:
                q = getattr(r, prop)
                line += " | %s %s %s" % (prop, q.Value, q.Unit)
            except Exception:
                pass
        try:
            line += " | mode %s" % r.Mode
        except Exception:
            pass
        vals.write(line + "\n")
        export(r, r.Name.replace(" ", "_"))
except Exception:
    log.write("FAILED\n%s\n" % traceback.format_exc())
log.close()
vals.close()
'''


def resolve(tags):
    """Re-solve (and re-export) selected systems of an existing ansys/corrected project."""
    proj = os.path.join(REPO, "ansys", "corrected", "wing_fea_corrected.wbpj")
    work = tempfile.mkdtemp(prefix="wbresolve_")
    out = os.path.join(work, "out")
    os.makedirs(out)
    with open(os.path.join(work, "mech.py"), "w") as f:
        f.write(MECH_SCRIPT)
    esc = lambda p: p.replace("\\", "\\\\")  # noqa: E731
    systems = [s for s in SYSTEMS if s[1] in tags]
    journal = f'''# encoding: utf-8
SetScriptVersion(Version="26.1.174")
log = open(r"{os.path.join(out, 'wb_log_resolve.txt')}", "w")
try:
    Open(FilePath=r"{proj}")
    src = open(r"{os.path.join(work, 'mech.py')}").read()
    for name, tag in {systems!r}:
        model = GetSystem(Name=name).GetContainer(ComponentName="Model")
        model.Edit()
        hdr = 'OUT_DIR = r"{esc(out)}"\\nSNIPPET_DIR = r"{esc(SNIP)}"\\nSYS_TAG = "%s"\\nRESOLVE_ONLY = True\\n' % tag
        model.SendCommand(Language="Python", Command=hdr + src)
        model.Exit()
        log.write("done %s\\n" % name)
    Save(Overwrite=True)
    log.write("saved\\n")
except Exception as e:
    log.write("FAILED: %s\\n" % e)
log.close()
'''
    jpath = os.path.join(work, "resolve.wbjn")
    with open(jpath, "w") as f:
        f.write(journal)
    subprocess.run([RUNWB2, "-B", "-R", jpath], cwd=work, capture_output=True, text=True, timeout=5400)
    dest_res = os.path.join(REPO, "results", "corrected")
    for fn in sorted(os.listdir(out)):
        shutil.copy(os.path.join(out, fn), os.path.join(dest_res, fn))
        if fn.endswith(".txt"):
            print(f"--- {fn}\n" + open(os.path.join(out, fn), errors="ignore").read(), flush=True)


def main():
    if not os.path.exists(RUNWB2):
        sys.exit(f"RunWB2 not found at {RUNWB2}")
    if "--resolve" in sys.argv:
        resolve(sys.argv[sys.argv.index("--resolve") + 1:])
        return
    os.makedirs(SNIP, exist_ok=True)
    with open(os.path.join(SNIP, "root_clamp.inp"), "w", newline="\n") as f:
        f.write(CLAMP)
    with open(os.path.join(SNIP, "design_load.inp"), "w", newline="\n") as f:
        f.write(design_load_snippet())

    work = tempfile.mkdtemp(prefix="wbfix_")
    out = os.path.join(work, "out")
    os.makedirs(out)
    shutil.copy(os.path.join(REPO, "ansys", "wing_fea.wbpj"), os.path.join(work, "wing_fea_corrected.wbpj"))
    shutil.copytree(os.path.join(REPO, "ansys", "wing_fea_files"), os.path.join(work, "wing_fea_corrected_files"))
    with open(os.path.join(work, "mech.py"), "w") as f:
        f.write(MECH_SCRIPT)
    esc = lambda p: p.replace("\\", "\\\\")  # noqa: E731
    journal = f'''# encoding: utf-8
SetScriptVersion(Version="26.1.174")
log = open(r"{os.path.join(out, 'wb_log.txt')}", "w")
try:
    Open(FilePath=r"{os.path.join(work, 'wing_fea_corrected.wbpj')}")
    src = open(r"{os.path.join(work, 'mech.py')}").read()
    for name, tag in {list(SYSTEMS)!r}:
        model = GetSystem(Name=name).GetContainer(ComponentName="Model")
        model.Edit()
        hdr = 'OUT_DIR = r"{esc(out)}"\\nSNIPPET_DIR = r"{esc(SNIP)}"\\nSYS_TAG = "%s"\\nRESOLVE_ONLY = False\\n' % tag
        model.SendCommand(Language="Python", Command=hdr + src)
        model.Exit()
        log.write("done %s\\n" % name)
    Save(Overwrite=True)
    log.write("saved\\n")
except Exception as e:
    log.write("FAILED: %s\\n" % e)
log.close()
'''
    jpath = os.path.join(work, "apply.wbjn")
    with open(jpath, "w") as f:
        f.write(journal)
    print("working folder:", work, flush=True)
    subprocess.run([RUNWB2, "-B", "-R", jpath], cwd=work, capture_output=True, text=True, timeout=5400)

    dest_res = os.path.join(REPO, "results", "corrected")
    os.makedirs(dest_res, exist_ok=True)
    for fn in sorted(os.listdir(out)):
        shutil.copy(os.path.join(out, fn), os.path.join(dest_res, fn))
        if fn.endswith(".txt"):
            print(f"--- {fn}\n" + open(os.path.join(out, fn), errors="ignore").read(), flush=True)
    dest_proj = os.path.join(REPO, "ansys", "corrected")
    if os.path.exists(os.path.join(work, "wing_fea_corrected.wbpj")):
        shutil.rmtree(dest_proj, ignore_errors=True)
        os.makedirs(dest_proj)
        shutil.copy(os.path.join(work, "wing_fea_corrected.wbpj"), dest_proj)
        shutil.copytree(os.path.join(work, "wing_fea_corrected_files"),
                        os.path.join(dest_proj, "wing_fea_corrected_files"))
        print("corrected project copied to", dest_proj)


if __name__ == "__main__":
    main()
