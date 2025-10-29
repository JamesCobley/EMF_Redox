# Local Electromagnetic Fields Enable Fast Redox Sensing  
### *James N. Cobley, University of Dundee, Scotland, UK*  
**Contact:** [j_cobley@yahoo.com](mailto:j_cobley@yahoo.com)  
**ORCID:** [0000-0001-5505-7424](https://orcid.org/0000-0001-5505-7424)  
**License:** MIT  

---

## 📘 Overview

This repository contains the full code and computational framework supporting the preprint:  

> **Cobley, J.N. (2025)**  
> *Local Electromagnetic Fields Enable Fast Redox Sensing by Physically Accelerating Cysteine Oxidation.*
> available at online[https://arxiv.org/abs/2510.24649]

The work demonstrates that **local electromagnetic fields (EMFs)**—ubiquitous in proteins, membranes, and nanodomains—can **lawfully modulate the Eyring barrier** and **orient reactants**, accelerating cysteine oxidation by orders of magnitude without altering the underlying chemistry.  

By embedding a field term into the Eyring equation, this framework shows how biologically plausible EMFs (10⁸–10⁹ V m⁻¹) reconcile the discrepancy between **bulk rate constants** and **cellular redox sensing timescales**.

---

## 🧠 Conceptual Summary

Traditional rate constants assume a **field-free, isotropic solvent**, but inside cells, **local EMFs** deform the potential energy surface and orient molecules in space.  
The field-aware Eyring equation used here is:

\[
k(E) = \kappa_0 \frac{k_B T}{h} 
\exp\left[-\frac{(\Delta G_0^{\ddagger} - \Delta \mu E - \frac{1}{2}\Delta \alpha E^2)}{RT}\right]
\]

Where:  
- **Δμ** – dipole moment change (10–20 D)  
- **Δα** – polarizability change (Å³)  
- **E** – local electric field (10⁷–10⁹ V m⁻¹)  
- **Fε** – dielectric focusing factor (1–5)  

The result is an exponential field-driven enhancement of reaction rates:  
\[
\frac{k(E)}{k_0} = G_{orient}\,\exp\left(\frac{\Delta \mu E + \frac{1}{2}\Delta \alpha E^2}{RT}\right)
\]

Local EMFs thus act as *intrinsic electrostatic catalysts*—accelerating thiolate oxidation, aligning reactive geometries, and stabilizing transition states.

---

## 🧩 Repository Contents

| File | Description |
|------|--------------|
| **`First_order_kinetics_H2O2.py`** | Computes the expected oxidation times under pseudo–first-order kinetics for given [H₂O₂] and rate constants. |
| **`Reverse_kinetic_calculator.py`** | Determines the [H₂O₂] or rate enhancement required to achieve a given fractional oxidation (e.g. 5%) within a target time. |
| **`Local_EMFs.py`** | Models local electromagnetic field magnitudes near membranes and charged residues under physiological conditions (Debye–Hückel screening). |
| **`Local_EMF_rates.py`** | Embeds the EMF into the Eyring equation to compute field-dependent rate constants across parameter sweeps. |
| **`EMF_accelerated_rates.py`** | Translates field-enhanced rate constants into observable redox sensing timescales. |
| **`t5_sweep_EMF_enhancements.csv`** | Tabulated results of field–rate sweeps showing time to 5% oxidation for various enhancement factors. |
| **`field_sweep_results.csv`** | Data of rate constant fold-accelerations across Δμ, Δα, and dielectric focusing parameters. |
| **`t5_summary_selected_doses.csv`** | Summary table comparing kinetic enhancements under biologically relevant [H₂O₂]. |
| **`Requirements.txt`** | Dependencies for Python runtime. |
| **`LICENSE`** | MIT license for reuse and modification. |

---

## 🧮 Running the Code

All scripts are self-contained and compatible with a **standard Google Colab** or local **Python 3.9+** environment.

### Clone the repository
```bash
git clone https://github.com/JamesCobley/EMF_Redox.git
cd EMF_Redox

pip install -r Requirements.txt
python Local_EMF_rates.py


