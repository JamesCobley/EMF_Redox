#@title Reverse calculator: [H2O2] or k2 needed for 5% oxidation in 5 min
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Targets ---
f_target = 0.05          #@param {type:"number"}  # target fractional oxidation (5%)
t_min = 5                #@param {type:"number"}  # target time (minutes)
t_sec = t_min * 60

# --- Mode A: given k2, solve for required [H2O2] ---
k2 = 10.0                #@param {type:"number"}  # M^-1 s^-1

# --- Mode B: given [H2O2], solve for required k2 ---
H2O2_uM_for_k2 = 1.0     #@param {type:"number"}  # µM

# --- Context (HeLa volume for molecule counts) ---
hela_volume_pL = 3.0     #@param {type:"number"}  # HeLa ~3 pL
NA = 6.022e23

# ---------- Closed-form solutions ----------
# From f = 1 - exp(-k2 * [H2O2] * t)
# (a) required [H2O2] = -ln(1-f) / (k2 * t)
H2O2_req_M = -np.log(1 - f_target) / (k2 * t_sec) if k2 > 0 else np.nan
H2O2_req_uM = H2O2_req_M * 1e6

# (b) required k2 = -ln(1-f) / (t * [H2O2])
H2O2_given_M = H2O2_uM_for_k2 * 1e-6
k2_req = -np.log(1 - f_target) / (t_sec * H2O2_given_M) if H2O2_given_M > 0 else np.nan

# Molecule-count context at the required [H2O2]
cell_L = hela_volume_pL * 1e-12
mol_per_cell_req = H2O2_req_M * cell_L * NA

print("=== Reverse calculations for 5% in 5 min ===")
print(f"(a) Given k2 = {k2:.3f} M^-1 s^-1 → required [H2O2] = {H2O2_req_uM:.3f} µM")
print(f"    (~{mol_per_cell_req:,.0f} H2O2 molecules present in a {hela_volume_pL} pL cell at that concentration)")
print(f"(b) Given [H2O2] = {H2O2_uM_for_k2:.3f} µM → required k2 = {k2_req:.1f} M^-1 s^-1")

# ---------- Optional sweeps for intuition ----------
# Sweep 1: required [H2O2] vs k2
k2_sweep = np.logspace(-1, 3, 200)          # 0.1 to 1000 M^-1 s^-1
H2O2_req_sweep_uM = -np.log(1 - f_target) / (k2_sweep * t_sec) * 1e6

plt.figure(figsize=(6,4))
plt.loglog(k2_sweep, H2O2_req_sweep_uM)
plt.xlabel("k$_2$ (M$^{-1}$ s$^{-1}$)")
plt.ylabel("Required [H$_2$O$_2$] (µM)")
plt.title("Required [H$_2$O$_2$] for 5% in 5 min")
plt.tight_layout()
plt.show()

# Sweep 2: required k2 vs [H2O2]
H2O2_sweep_uM = np.logspace(-3, 1, 200)     # 1 nM to 10 µM
k2_req_sweep = -np.log(1 - f_target) / (t_sec * H2O2_sweep_uM * 1e-6)

plt.figure(figsize=(6,4))
plt.loglog(H2O2_sweep_uM, k2_req_sweep)
plt.xlabel("[H$_2$O$_2$] (µM)")
plt.ylabel("Required k$_2$ (M$^{-1}$ s$^{-1}$)")
plt.title("Required k$_2$ for 5% in 5 min")
plt.tight_layout()
plt.show()
