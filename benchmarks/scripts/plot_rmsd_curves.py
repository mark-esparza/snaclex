"""Figure 5: cumulative redocking success curves (needs numpy + matplotlib)."""
import json, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE=os.path.dirname(os.path.abspath(__file__))
W=os.path.join(HERE,"..","results")+os.sep
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11})
fig,axes=plt.subplots(1,2,figsize=(10,4.2),sharey=True)
th=np.linspace(0,10,201)
for ax,(s,name) in zip(axes,[("astex_diverse_set","Astex Diverse (n = 85)"),("posebusters_benchmark_set","PoseBusters Benchmark (n = 428)")]):
    ok=json.load(open(W+f"{s}_with_sym.json"))
    single=np.array([[r["runs"][k]["rmsd_sym_A"] for k in range(3)] for r in ok])
    best=np.array([min(r["runs"],key=lambda x:x["score"])["rmsd_sym_A"] for r in ok])
    anyb=single.min(axis=1)
    per_seed=np.array([[100*(single[:,k]<=t).mean() for t in th] for k in range(3)])
    ax.fill_between(th,per_seed.min(0),per_seed.max(0),color="#9ecae1",alpha=.6,lw=0,label="Single run (range of 3 seeds)")
    ax.plot(th,per_seed.mean(0),color="#3182bd",lw=1.6)
    ax.plot(th,[100*(best<=t).mean() for t in th],color="#e6550d",lw=2,label="Best score of 3 runs")
    ax.plot(th,[100*(anyb<=t).mean() for t in th],color="#636363",lw=1.4,ls="--",label="Lowest RMSD of 3 runs")
    ax.axvline(2,color="k",lw=.8,ls=":")
    ax.set_title(name,fontsize=12); ax.set_xlabel("Symmetry corrected heavy atom RMSD (Å)")
    ax.set_xlim(0,10); ax.set_ylim(0,100); ax.grid(alpha=.25)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("Complexes within threshold (%)")
axes[1].legend(loc="lower right",frameon=False,fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(HERE,"..","figures","fig5_benchmark.png"),dpi=300); print("ok")
