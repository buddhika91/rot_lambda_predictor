# ROT Λ Predictor

**A boundary-mode counting and horizon-thermodynamic mapping tool  
for computing the cosmological constant Λ within the Recursive Observation Theory (ROT) framework.**

---

## Overview

This repository provides a **single-file reference implementation** that connects:

- **ROT boundary Planck-saturation closure**
- **Cosmological horizon thermodynamics**
- **General Relativity ΛCDM normalization**

to compute the cosmological constant **Λ** in both:

- **time units:** `1/s²`  
- **length units:** `1/m²`

The script is designed to be:

- Minimal  
- Transparent  
- Numerically stable  
- Suitable for research notes or reproducibility supplements  

---

## Physical Idea (High-Level)

The predictor combines three ingredients:

### 1. Boundary mode-counting closure (ROT)

ROT introduces a boundary invariant:

\[
I = \frac{\kappa}{\ell_\kappa^2}
\]

Numerical RG-style relaxation selects the **Planck-saturated fixed point**:

\[
I \rightarrow \frac{1}{\ell_P^2}
\]

This removes any free ultraviolet cutoff.

---

### 2. Horizon thermodynamic mapping

Using standard cosmological horizon quantities:

- \(R_H = c/H\)  
- \(T_H = \hbar H / (2\pi k_B)\)  
- \(S_H = k_B A_H / \ell_P^2\)

vacuum energy density is assigned via:

\[
\rho_\Lambda = \frac{E}{V_H}, \quad E = f\,T_H S_H
\]

where **f** encodes the quasi-local horizon energy convention.

---

### 3. Relation to ΛCDM

General Relativity defines:

\[
\Lambda = \frac{3\,\Omega_\Lambda H_0^2}{c^2}
\]

The **GR-matched horizon convention** corresponds to:

\[
f = \tfrac{1}{4}
\]

which reproduces the observed Λ scale exactly when  
\(H_0\) and \(\Omega_\Lambda\) are supplied.

---

## Features

- Planck-saturated boundary entropy (**no UV tuning**)
- Horizon thermodynamic Λ derivation
- GR-consistent normalization option
- Exact agreement with ΛCDM relation
- No external dependencies (**pure Python**)
- Single-file reproducibility

---

## Requirements

- **Python 3.9+**
- **No third-party libraries**

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/rot-lambda-predictor.git
cd rot-lambda-predictor
