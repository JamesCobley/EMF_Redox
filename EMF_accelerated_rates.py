#@title t5% vs [H2O2] for EMF-accelerated k2 (baseline 10 M^-1 s^-1)
# Computes time to 5% oxidation across [H2O2] = 1 nM → 1 µM
# for rate-constant enhancements (1× … 1000×). Saves figure + CSVs.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Settings ----------
k2_baseline = 10.0  # M^-1 s^-1 (field-free)
folds = [1, 3, 10, 30, 100, 300, 1000]  # EMF-driven accelerations to test
n_points = 300
H2O2_min_nM, H2O2_max_uM = 1.0, 1.0  # sweep range
f_target = 0.05  # 5% oxidation target

# ---------- Build concentration grid ----------
H2O2_grid_M = np.logspace(np.log10(H2O2_min_nM*1e-9),
                          np.log10(H2O2_max_uM*1e-6),
                          n_points)

# ---------- Compute t5% ----------
ln_term = -np.log(1.0 - f_target)  # ~0.051293
rows = []
for fold in folds:
    k2_eff = k2_baseline * fold  # M^-1 s^-1
    kprime = k2_eff * H2O2_grid_M
    t5_sec = ln_term / kprime
    for cM, tsec in zip(H2O2_grid_M, t5_sec):
        rows.append({
            "H2O2_uM": cM * 1e6,
            "fold_enhancement": fold,
            "k2_eff_Minv_sinv": k2_eff,
            "t5_seconds": tsec,
            "t5_minutes": tsec / 60.0
        })

df = pd.DataFrame(rows)

# ---------- Save CSV with full grid ----------
df.to_csv("t5_sweep_EMF_enhancements.csv", index=False)

# ---------- Plot ----------
plt.figure(figsize=(7,5), dpi=200)
for fold in folds:
    sub = df[df["fold_enhancement"] == fold]
    plt.loglog(sub["H2O2_uM"], sub["t5_minutes"], label=f"{fold}×")

plt.xlabel("[H$_2$O$_2$] (µM)")
plt.ylabel("Time to 5% oxidation (min)")
plt.title("t$_{5\\%}$ vs [H$_2$O$_2$] for EMF-accelerated k$_2$ (baseline 10 M$^{-1}$ s$^{-1}$)")
plt.grid(which='both', ls='--', lw=0.4, alpha=0.6)
plt.legend(title="Rate enhancement", ncol=2, fontsize=9)
plt.tight_layout()
plt.savefig("t5_vs_H2O2_by_enhancement.png", dpi=600, bbox_inches="tight", transparent=True)
plt.show()

# ---------- Small summary table at three doses (optional) ----------
def t5_at(conc_uM):
    cM = conc_uM * 1e-6
    out = []
    for fold in folds:
        k2_eff = k2_baseline * fold
        t5_min = ln_term / (k2_eff * cM) / 60.0
        out.append((fold, t5_min))
    return pd.DataFrame(out, columns=["fold_x", f"t5_min_at_{conc_uM}_uM"])

summary = t5_at(0.01).merge(t5_at(0.1), on="fold_x").merge(t5_at(1.0), on="fold_x")
summary.to_csv("t5_summary_selected_doses.csv", index=False)
summary
