# Local electric field near a plasma membrane + nearby charged residues
# Assumptions: distances 1–5 nm from membrane; screened Coulomb; membrane surface potential;
# dielectric focusing for low-ε protein environment around the reactive site.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- constants ----
eps0 = 8.8541878128e-12; kB = 1.380649e-23; T = 310.0
NA = 6.02214076e23; e = 1.602176634e-19

# ---- user parameters (edit freely) ----
eps_r_bulk   = 80.0         # water
eps_r_shell  = 10.0         # local low-ε shell around cysteine
F_eps        = eps_r_bulk/eps_r_shell   # dielectric focusing factor (normal component)

ionic_strength_M = 0.15      # physiological ionic strength (M)
phi0 = -0.070                # membrane surface potential (V), e.g., -70 mV
z_nm = np.linspace(1.0, 5.0, 41)  # distance along membrane normal (nm)

# Positively charged residues around cysteine (positions in nm relative to cysteine)
residue_positions_nm = np.array([
    [0.6,  0.0, 0.2],
    [-0.3, 0.5, -0.1],
    [0.0, -0.6, 0.3],
])
residue_charges_e = np.array([+1, +1, +1])  # Lys/Arg

# ---- Debye length ----
def debye_length_m(eps_r, I_M, T_K=310.0):
    # λ_D = sqrt( (ε_r ε0 kT) / (2 N_A e^2 I * 1000) )
    return np.sqrt((eps_r*eps0*kB*T_K)/(2*NA*(e**2)*I_M*1000.0))

lambda_D = debye_length_m(eps_r_bulk, ionic_strength_M, T)

# ---- fields ----
def E_membrane(z_m, phi0, lambda_D):
    # Normal field from diffuse layer: E(z) ≈ (phi0/λ_D) * exp(-z/λ_D)
    return (phi0/lambda_D) * np.exp(-z_m/lambda_D)

def E_point_charge_vec(r_vec_m, q_C, eps_r, lambda_D):
    r = np.linalg.norm(r_vec_m, axis=-1, keepdims=True)
    r = np.maximum(r, 1e-12)
    kappa = 1.0/lambda_D
    pref = (1.0/(4.0*np.pi*eps0*eps_r)) * q_C * np.exp(-kappa*r) * (1.0 + kappa*r) / (r**3)
    return pref * r_vec_m

# Convert nm→m
z_m = z_nm*1e-9
res_pos_m = residue_positions_nm*1e-9

# Membrane field (normal); apply dielectric focusing
E_mem = E_membrane(z_m, phi0, lambda_D)
E_mem_focused = F_eps * E_mem

# Residue field (vector sum at cysteine; project onto normal and focus)
E_res_vec = np.zeros((z_m.size, 3))
for pos, q_e in zip(res_pos_m, residue_charges_e):
    r_vec = np.broadcast_to(-pos, E_res_vec.shape)  # displacement from charge to site
    E_res_vec += E_point_charge_vec(r_vec, q_e*e, eps_r_bulk, lambda_D)

E_res_normal   = E_res_vec[:, 2]
E_res_focused  = F_eps * E_res_normal

E_total_normal = E_mem_focused + E_res_focused

# ---- export + plot ----
df = pd.DataFrame({
    "z_nm": z_nm,
    "lambda_D_nm": np.full_like(z_nm, lambda_D*1e9),
    "E_mem_V_per_m": E_mem,
    "E_mem_focused": E_mem_focused,
    "E_res_normal_V_per_m": E_res_normal,
    "E_res_focused": E_res_focused,
    "E_total_normal_V_per_m": E_total_normal,
})
df.to_csv("local_field_membrane_plus_residues.csv", index=False)

plt.figure(figsize=(7,5), dpi=200)
plt.semilogy(z_nm, np.abs(E_mem_focused), label="Membrane (focused)")
plt.semilogy(z_nm, np.abs(E_res_focused), label="Pos. residues (focused normal)")
plt.semilogy(z_nm, np.abs(E_total_normal), label="Total normal field")
plt.xlabel("Distance from membrane (nm)")
plt.ylabel("|E| (V m$^{-1}$)")
plt.title(f"Local normal electric field (φ₀={phi0*1e3:.0f} mV, I={ionic_strength_M*1e3:.0f} mM, ε_shell={eps_r_shell}, Fε≈{F_eps:.1f})")
plt.grid(which='both', ls='--', lw=0.4, alpha=0.6)
plt.legend(); plt.tight_layout()
plt.savefig("local_field_vs_distance.png", dpi=600, bbox_inches="tight", transparent=True)
plt.show()
