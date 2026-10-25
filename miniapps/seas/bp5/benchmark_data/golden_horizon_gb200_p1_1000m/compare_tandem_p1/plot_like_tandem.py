"""Plots for cmp_like_tandem.py output: Tandem p=1 drawn first, CPU next,
GPU last (on top).  Usage: plot_like_tandem.py GPU_DIR CPU_DIR TANDEM_DIR CMP_DIR OUT.pdf PNG_DIR"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

GPU, CPU, TDM, CMP, OUT, PNG = sys.argv[1:7]
os.makedirs(PNG, exist_ok=True)
YR = 365.25 * 86400.0
STATIONS = ["strk+00dp+10", "strk+16dp+10", "strk-16dp+10", "strk-24dp+10", "strk+00dp+22",
            "strk+00dp+00", "strk+16dp+00", "strk-16dp+00", "strk+36dp+00", "strk-36dp+00"]
STY = {  # draw order: Tandem, CPU, GPU (GPU on top)
    "Tandem": dict(color="#2f5fa7", lw=3.0, alpha=0.45, ls="-", zorder=1, label="Tandem p=1"),
    "CPU":    dict(color="black", lw=1.4, alpha=1.0, ls="-", zorder=2, label="CPU 1074099 (Vista, 144 ranks)"),
    "GPU":    dict(color="#d0201a", lw=1.0, alpha=1.0, ls="--", zorder=3, label="GPU 69011 (Horizon, 1 GB200)"),
}
ORDER = ["Tandem", "CPU", "GPU"]
PAIRS = {"gpu_vs_cpu": ("GPU - CPU", "#d0201a"), "gpu_vs_tandem": ("GPU - Tandem", "#e08a1e"),
         "cpu_vs_tandem": ("CPU - Tandem", "#555555")}

def tname(s):
    return f"bp5qd_tandem_p1_x2_{int(s[4:7])}_x3_{int(s[9:12])}.txt"
def load(code, s):
    if code == "Tandem":
        a = np.array([[float(v) for v in ln.split()] for ln in open(os.path.join(TDM, tname(s)))
                      if ln.strip() and ln.strip()[0] not in "#t"])
    else:
        a = np.loadtxt(os.path.join(GPU if code == "GPU" else CPU, f"bp5_cycle_p1_fltst_{s}.txt"), comments="#")
    return a
D = {c: {s: load(c, s) for s in STATIONS} for c in ORDER}
S = json.load(open(os.path.join(CMP, "summary.json")))
P = S["pairs"]
ON = {"GPU": [e["GPU"] for e in P["gpu_vs_cpu"]["event_onsets"]],
      "CPU": [e["CPU"] for e in P["gpu_vs_cpu"]["event_onsets"]],
      "Tandem": [e["Tandem"] for e in P["gpu_vs_tandem"]["event_onsets"]]}
COLS = {1: "slip_strike (m)", 3: "log10 V_strike (m/s)", 5: "tau_strike (MPa)", 7: "log10 state (s)"}
pngs = []

def save(fig, pdf, name):
    pdf.savefig(fig)
    fig.savefig(os.path.join(PNG, name + ".png"), dpi=105)
    pngs.append(name); plt.close(fig)

with PdfPages(OUT) as pdf:
    # page: event onset differences + per-window interseismic maxima at the 6 center/depth stations
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    ev = np.arange(1, 8)
    w = 0.26
    for i, (tag, (lab, col)) in enumerate(PAIRS.items()):
        d = [e["diff_s"] for e in P[tag]["event_onsets"]]
        axs[0, 0].bar(ev + (i - 1) * w, d, w, color=col, label=lab)
    axs[0, 0].axhline(0, color="k", lw=0.6)
    axs[0, 0].set(xlabel="event", ylabel="onset difference (s)", title="Event onsets (earliest station crossing of log10 V = -3)")
    axs[0, 0].legend(fontsize=8); axs[0, 0].grid(alpha=0.3)
    DEPTH6 = ["strk+00dp+00", "strk+00dp+10", "strk+00dp+22", "strk-16dp+10", "strk+16dp+10", "strk-24dp+10"]
    wins = sorted({r["window"] for r in P["gpu_vs_cpu"]["interseismic"]})
    for ax, f, ylab in ((axs[0, 1], "max_abs_dtau_MPa", "max |d tau| (MPa)"),
                        (axs[1, 0], "max_abs_dslip_m", "max |d slip| (m)"),
                        (axs[1, 1], "max_abs_dlog10V", "max |d log10 V|")):
        for tag, (lab, col) in PAIRS.items():
            y = [max(r[f] for r in P[tag]["interseismic"] if r["window"] == wn and r["station"] in DEPTH6) for wn in wins]
            ax.semilogy(range(len(wins)), y, "o-", color=col, label=lab)
        ax.set_xticks(range(len(wins))); ax.set_xticklabels(wins)
        ax.set(xlabel="interseismic window (onset+1 yr to next onset-1 day)", ylabel=ylab,
               title="Between events, 6 center/depth stations")
        ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
    fig.suptitle("BP5-QD 1000 m p=1: GPU 69011 vs CPU 1074099 vs Tandem p=1 (compare_tandem_p1 method, 365.25-day years)", fontsize=10)
    fig.tight_layout(); save(fig, pdf, "summary")

    # event zooms, each code aligned on its own event onset
    for s in ("strk+00dp+10", "strk+36dp+00"):
        fig, axs = plt.subplots(2, 4, figsize=(11, 8.5))
        for j, k in enumerate((1, 6)):
            for code in ORDER:
                a = D[code][s]; t0 = ON[code][k] * YR
                m = (a[:, 0] > t0 - 60) if s.endswith("dp+10") else (a[:, 0] > t0 - 600)
                m &= a[:, 0] < t0 + (600 if s.endswith("dp+10") else 4000)
                for q, c in enumerate((1, 3, 5, 7)):
                    axs[j, q].plot(a[m, 0] - t0, a[m, c], **{**STY[code], "label": code})
                    axs[j, q].set(xlabel=f"t - own onset of event {k+1} (s)", ylabel=COLS[c]); axs[j, q].grid(alpha=0.3)
            axs[j, 0].legend(fontsize=7)
        fig.suptitle(f"{s}: events 2 and 7, each code aligned on its own onset (GPU drawn on top)", fontsize=10)
        fig.tight_layout(); save(fig, pdf, f"events_{s}")

    # one page per station: full cycle + per-window between-event maxima
    for s in STATIONS:
        fig, axs = plt.subplots(3, 2, figsize=(11, 8.5))
        for ax, c in zip(axs.flat[:4], (1, 3, 5, 7)):
            for code in ORDER:
                a = D[code][s]
                ax.plot(a[:, 0] / YR, a[:, c], **STY[code])
            ax.set(xlabel="time (yr)", ylabel=COLS[c]); ax.grid(alpha=0.3)
        axs[0, 0].legend(fontsize=7, loc="upper left")
        for ax, f, ylab in ((axs[2, 0], "max_abs_dtau_MPa", "between events max |d tau| (MPa)"),
                            (axs[2, 1], "max_abs_dslip_m", "between events max |d slip| (m)")):
            for tag, (lab, col) in PAIRS.items():
                rows = [r for r in P[tag]["interseismic"] if r["station"] == s]
                ax.semilogy(range(len(rows)), [r[f] for r in rows], "o-", color=col, label=lab)
                ax.set_xticks(range(len(rows))); ax.set_xticklabels([r["window"] for r in rows])
            ax.set(xlabel="window", ylabel=ylab); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7)
        fig.suptitle(f"{s}: Tandem p=1 (blue, bottom), CPU 1074099 (black), GPU 69011 (red dashed, top). "
                     "Bottom row: compare_tandem_p1 between-event maxima per window", fontsize=8.5)
        fig.tight_layout(); save(fig, pdf, f"station_{s}")
print("wrote", OUT, "pngs:", " ".join(pngs))
