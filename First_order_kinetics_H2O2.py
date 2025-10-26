import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import log

# --- User parameters ---
k2 = 10.0  #@param {type:"number"}  # second-order rate constant (M^-1 s^-1)
conc_min_nM = 1  #@param {type:"number"}
conc_max_uM = 1  #@param {type:"number"}
n_points = 300   #@param {type:"integer"}
expr_copies = 100000  #@param {type:"integer"}
hela_volume_pL = 3.0  #@param {type:"number"}

# --- Computation ---
NA = 6.022e23
conc_range_M = np.logspace(np.log10(conc_min_nM*1e-9), np.log10(conc_max_uM*1e-6), n_points)
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

# --- Display ---
from google.colab import data_table
data_table.enable_dataframe_formatter()

print("Time to 5% oxidation is independent of copy number under pseudo-first-order kinetics (with replenished H2O2).")
display(df.head())

# Plot
plt.figure(figsize=(6,4))
plt.loglog(df["[H2O2]_µM"], df["t5_minutes"])
plt.xlabel("[H$_2$O$_2$] (µM)")
plt.ylabel("Time to 5% oxidation (min)")
plt.tight_layout()
plt.savefig("pseudo-first-order-kinetics.png", dpi=300)
plt.show()
