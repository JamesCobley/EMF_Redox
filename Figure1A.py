import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- User parameters ---
k2 = 10.0        # second-order rate constant (M^-1 s^-1)
conc_min_nM = 1  # min [H2O2] in nM
conc_max_uM = 1  # max [H2O2] in µM
n_points = 300
expr_copies = 100000
hela_volume_pL = 3.0

# --- Computation ---
NA = 6.022e23
conc_range_M = np.logspace(np.log10(conc_min_nM*1e-9),
                           np.log10(conc_max_uM*1e-6), n_points)
kprime = k2 * conc_range_M
t5 = -np.log(0.95) / kprime  # seconds

df = pd.DataFrame({
    "[H2O2]_µM": conc_range_M*1e6,
    "k_prime_s^-1": kprime,
    "t5_seconds": t5,
    "t5_minutes": t5/60,
    "Copies": expr_copies,
    "HeLaVolume_pL": hela_volume_pL,
    "H2O2_molecules_present": conc_range_M * NA * (hela_volume_pL*1e-12),
})

# --- Plot ---
plt.figure(figsize=(3.5, 3.0), dpi=300)  # ~single-column width
ax = plt.gca()
ax.loglog(df["[H2O2]_µM"], df["t5_minutes"], lw=2)

# Labels & styling
ax.set_xlabel(r'[H$_2$O$_2$] (µM)', fontsize=12)
ax.set_ylabel(r'Time to 5% oxidation(min)', fontsize=12)
ax.tick_params(labelsize=10)
ax.grid(which='both', ls='--', lw=0.4, alpha=0.6)

# Framing
ax.set_xlim(1e-3, 1)
ax.set_ylim(1e2, 1e5)

# Interpretive annotation 
ax.annotate('Bulk kinetics too slow',
            xy=(1e-2, 1.5e4), xytext=(8e-2, 4e4),
            arrowprops=dict(arrowstyle='->', lw=0.8),
            fontsize=9)

plt.tight_layout()
plt.savefig("Fig1A_bulk_kinetics.png", dpi=600, bbox_inches='tight', transparent=True)
plt.show()
