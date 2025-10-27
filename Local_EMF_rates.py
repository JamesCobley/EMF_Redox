# Field-aware Eyring parameter sweep (reproducible, Colab-ready)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- Constants ----------
kB  = 1.380649e-23      # J/K
NA  = 6.02214076e23     # 1/mol
R   = 8.314462618       # J/(mol K)
T   = 310.0             # K
eps0 = 8.8541878128e-12 # F/m
DEBYE_TO_C_M = 3.33564e-30
KBT_J = kB*T

# ---------- Helper functions ----------
def orient_gain(mu_eff_D, E, theta_c_deg):
    """G_orient(ξ, θc)/P0 with ξ = μ_eff E / (kB T)."""
    mu_eff = mu_eff_D*DEBYE_TO_C_M
    xi = (mu_eff*E)/KBT_J
    th = np.deg2rad(theta_c_deg)
    num = np.exp(xi) - np.exp(xi*np.cos(th))
    den = np.exp(xi) - np.exp(-xi)
    den = np.where(np.isclose(den,0.0), 1e-300, den)
    P = num/den
    P0 = (1 - np.cos(th))/2.0
    return P/P0

def barrier_gain(delta_mu_D, E, delta_alpha_A3=0.0):
    """exp[(Δμ E + 0.5 Δα E^2)/RT]; Δα input in Å^3, converted to SI."""
    dmu = delta_mu_D*DEBYE_TO_C_M
    dalpha_SI = 4.0*np.pi*eps0*(delta_alpha_A3*1e-30)
    exponent = (dmu*E + 0.5*dalpha_SI*E**2)/(R*T)
    return np.exp(exponent)

def fold_accel(E, delta_mu_D, mu_eff_D, theta_c_deg, delta_alpha_A3=0.0, F_eps=1.0):
    """Total acceleration = G_orient * barrier_gain; E scaled by dielectric focusing F_eps."""
    E_loc = F_eps*E
    return orient_gain(mu_eff_D, E_loc, theta_c_deg)*barrier_gain(delta_mu_D, E_loc, delta_alpha_A3)

# ---------- Sweep settings ----------
E_vals        = np.logspace(7, np.log10(3e8), 200)   # 1e7–3e8 V/m
delta_mu_list = [10.0, 15.0, 20.0]                   # Debye
theta_list    = [20.0, 30.0, 45.0]                   # degrees
mu_eff_list   = [1.0, 2.0, 3.0]                      # Debye
F_eps_list    = [1.0, 2.0, 3.0, 5.0]                 # dielectric focusing multipliers
delta_alpha_A3 = 0.0                                  # try 40.0 to test Δα contribution
k0 = 10.0                                             # baseline M^-1 s^-1

# ---------- Run sweep + save CSV ----------
rows = []
for E in E_vals:
    for dm in delta_mu_list:
        for th in theta_list:
            for me in mu_eff_list:
                for Fe in F_eps_list:
                    acc = fold_accel(E, dm, me, th, delta_alpha_A3, Fe)
                    rows.append({
                        "E_Vm": E, "DeltaMu_D": dm, "Theta_deg": th, "Mu_eff_D": me,
                        "F_eps": Fe, "FoldAccel": float(acc), "k_eff_Minv_sinv": float(k0*acc)
                    })
df = pd.DataFrame(rows)
df.to_csv("/content/field_sweep_results.csv", index=False)

# ---------- Figure A: Fold-acceleration vs E (vary Δμ; θc=30°, μ_eff=2D) ----------
plt.figure(figsize=(6,4), dpi=200)
for dm in delta_mu_list:
    acc = fold_accel(E_vals, dm, mu_eff_D=2.0, theta_c_deg=30.0, delta_alpha_A3=delta_alpha_A3, F_eps=1.0)
    plt.loglog(E_vals, acc, label=f"Δμ={dm} D, θc=30°, μ_eff=2 D")
plt.xlabel("Local field, E (V m$^{-1}$)")
plt.ylabel("Fold-acceleration, $k(E)/k_0$")
plt.grid(which='both', ls='--', lw=0.4, alpha=0.6)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig("/content/FigA_fold_accel_vs_E.png", dpi=600, bbox_inches="tight", transparent=True)
plt.show()

# ---------- Figure B: Dielectric focusing curves (Δμ=15 D, θc=30°, μ_eff=2D) ----------
plt.figure(figsize=(6,4), dpi=200)
for Fe in F_eps_list:
    acc = fold_accel(E_vals, delta_mu_D=15.0, mu_eff_D=2.0, theta_c_deg=30.0,
                     delta_alpha_A3=delta_alpha_A3, F_eps=Fe)
    plt.loglog(E_vals, acc, label=f"Fε={Fe}")
plt.xlabel("Local field, E (V m$^{-1}$)")
plt.ylabel("Fold-acceleration, $k(E)/k_0$")
plt.title("Dielectric focusing multiplies E")
plt.grid(which='both', ls='--', lw=0.4, alpha=0.6)
plt.legend(fontsize=8, title="Amplification")
plt.tight_layout()
plt.savefig("/content/FigB_focusing.png", dpi=600, bbox_inches="tight", transparent=True)
plt.show()

# ---------- Representative printed cases ----------
def report_case(E, dm, me, th, Fe):
    acc = fold_accel(E, dm, me, th, delta_alpha_A3, Fe)
    print(f"E={E:.2e} V/m, Δμ={dm} D, μ_eff={me} D, θc={th}°, Fε={Fe} -> "
          f"×{acc:,.1f}, k≈{k0*acc:,.1f} M^-1 s^-1")

print("=== Representative cases ===")
report_case(1e8, 15.0, 2.0, 30.0, 1.0)
report_case(1e8, 20.0, 3.0, 30.0, 1.0)
report_case(2e8, 15.0, 2.0, 30.0, 2.0)
report_case(3e8, 20.0, 3.0, 20.0, 3.0)

print("\nSaved CSV: /content/field_sweep_results.csv")
print("Saved Fig A: /content/FigA_fold_accel_vs_E.png")
print("Saved Fig B: /content/FigB_focusing.png")
