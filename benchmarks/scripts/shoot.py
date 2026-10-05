"""Rebuild figure pages from the template and screenshot them with headless Edge.
usage: python shoot.py fig1:0,0 fig1:90,0 ...   (writes figN_Y_X.png)"""
import os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(here, "fig_data_full.json")):
    t = open(os.path.join(here, "fig1.html"), encoding="utf-8").read()
    a = t.index("const D = ") + len("const D = ")
    b = t.index(', FIG = "', a)
    open(os.path.join(here, "fig_data_full.json"), "w", encoding="utf-8").write(t[a:b])
data = open(os.path.join(here, "fig_data_full.json"), encoding="utf-8").read()
tmpl = open(os.path.join(here, "fig_template.html"), encoding="utf-8").read()
edge = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
for spec in sys.argv[1:]:
    fig, rot = spec.split(":")
    page = os.path.join(here, f"{fig}.html")
    open(page, "w", encoding="utf-8").write(tmpl.replace("__DATA__", data).replace("__FIG__", fig))
    out = os.path.join(here, f"{fig}_{rot.replace(',', '_')}.png")
    url = "file:///" + page.replace("\\", "/") + "#" + rot
    try:
        subprocess.run([edge, "--headless=new", "--use-angle=swiftshader", "--enable-unsafe-swiftshader",
                        f"--user-data-dir={os.path.join(here, 'edge_profile')}", "--no-first-run",
                        "--hide-scrollbars", "--window-size=1600,1200", "--virtual-time-budget=15000",
                        f"--screenshot={out}", url], timeout=120, capture_output=True)
    except Exception as e:
        print("ERR", spec, e, flush=True)
    print(out, os.path.exists(out), flush=True)
