#Local Electromagnetic Fields Enable Fast Redox Sensing
James N. Cobley, University of Dundee, Scotland, UK

Contact: j_cobley@yahoo.com

ORCID: 0000-0001-5505-7424

License: MIT

📘 Overview

This repository contains the full code and computational framework supporting the manuscript:

Cobley, J.N. (2025)
Local Electromagnetic Fields Enable Fast Redox Sensing by Physically Accelerating Cysteine Oxidation.

The work demonstrates that local electromagnetic fields (EMFs)—ubiquitous in proteins, membranes, and nanodomains—can lawfully modulate the Eyring barrier and orient reactants, accelerating cysteine oxidation by orders of magnitude without altering the underlying chemistry.

By embedding a field term into the Eyring equation, this framework shows how biologically plausible EMFs (10⁸–10⁹ V m⁻¹) reconcile the discrepancy between bulk rate constants and cellular redox sensing timescales.

🧠 Conceptual Summary

Traditional rate constants assume a field-free, isotropic solvent, but inside cells, local EMFs deform the potential energy surface and orient molecules in space.
The field-aware Eyring equation used here is:

𝑘
(
𝐸
)
=
𝜅
0
𝑘
𝐵
𝑇
ℎ
exp
⁡
[
−
(
Δ
𝐺
0
‡
−
Δ
𝜇
𝐸
−
1
2
Δ
𝛼
𝐸
2
)
𝑅
𝑇
]
k(E)=κ
0
	​

h
k
B
	​

T
	​

exp[−
RT
(ΔG
0
‡
	​

−ΔμE−
2
1
	​

ΔαE
2
)
	​

]

Where:

Δμ – dipole moment change (10–20 D)

Δα – polarizability change (Å³)

E – local electric field (10⁷–10⁹ V m⁻¹)

Fε – dielectric focusing factor (1–5)

The result is an exponential field-driven enhancement of reaction rates:

𝑘
(
𝐸
)
𝑘
0
=
𝐺
𝑜
𝑟
𝑖
𝑒
𝑛
𝑡
 
exp
⁡
(
Δ
𝜇
𝐸
+
1
2
Δ
𝛼
𝐸
2
𝑅
𝑇
)
k
0
	​

k(E)
	​

=G
orient
	​

exp(
RT
ΔμE+
2
1
	​

ΔαE
2
	​

)

Local EMFs thus act as intrinsic electrostatic catalysts—accelerating thiolate oxidation, aligning reactive geometries, and stabilizing transition states.

🧩 Repository Contents
File	Description
First_order_kinetics_H2O2.py	Computes the expected oxidation times under pseudo–first-order kinetics for given [H₂O₂] and rate constants.
Reverse_kinetic_calculator.py	Determines the [H₂O₂] or rate enhancement required to achieve a given fractional oxidation (e.g. 5%) within a target time.
Local_EMFs.py	Models local electromagnetic field magnitudes near membranes and charged residues under physiological conditions (Debye–Hückel screening).
Local_EMF_rates.py	Embeds the EMF into the Eyring equation to compute field-dependent rate constants across parameter sweeps.
EMF_accelerated_rates.py	Translates field-enhanced rate constants into observable redox sensing timescales.
t5_sweep_EMF_enhancements.csv	Tabulated results of field–rate sweeps showing time to 5% oxidation for various enhancement factors.
field_sweep_results.csv	Data of rate constant fold-accelerations across Δμ, Δα, and dielectric focusing parameters.
t5_summary_selected_doses.csv	Summary table comparing kinetic enhancements under biologically relevant [H₂O₂].
Requirements.txt	Dependencies for Python runtime.
LICENSE	MIT license for reuse and modification.
🧮 Running the Code

All scripts are self-contained and compatible with a standard Google Colab or local Python 3.9+ environment.

Clone the repository
git clone https://github.com/JamesCobley/EMF_Redox.git
cd EMF_Redox

Install dependencies
pip install -r Requirements.txt

Example: Compute field-accelerated rates
python Local_EMF_rates.py


This will output the fold acceleration k(E)/k₀ across a range of field strengths and dipole changes, as plotted in the manuscript.

🧭 Reproducibility

All numerical values, field magnitudes, and kinetic constants are directly traceable to:

Empirical literature on cysteine oxidation and EMFs (see manuscript references)

Standard physical constants and dielectric models

Fully open, reproducible Python scripts hosted here under MIT license

Each script contains detailed inline comments describing equations and assumptions.

🔬 Key Findings

Local EMFs of 10⁷–10⁸ V m⁻¹ accelerate cysteine oxidation by 10–10³×.

Dielectric focusing (ε ≈ 10) further amplifies rate enhancement ∝ E².

Falsifiable predictions: Vibrational Stark shifts of 1–10 cm⁻¹ (MV cm⁻¹)⁻¹ in thiolate–peroxide complexes.

Broader implication: Rate constants are not immutable numbers—they are field-conditioned parameters reflecting local geometry, charge, and time.

🧾 Citation

If you use this repository, please cite:

Cobley, J.N.
Local Electromagnetic Fields Enable Fast Redox Sensing by Physically Accelerating Cysteine Oxidation.
(2025)

DOI and Zenodo citation will be provided upon release 
