#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rot_lambda_predictor_final.py
=================================

A repo-ready, single-file “final” predictor for the cosmological constant Λ
based on your ROT boundary-mode counting closure + horizon thermodynamic mapping.

What this script provides
-------------------------
1) **ROT boundary closure (Planck saturation)**:
      I := κ / ℓ_κ^2  →  1 / ℓ_P^2
   so the UV scale is selected (not chosen):
      ℓ_κ → sqrt(κ) * ℓ_P

2) **Horizon thermodynamic mapping**:
   Uses standard horizon quantities:
      R_H = c/H,  A_H = 4πR_H^2,  V_H = (4/3)πR_H^3,  T_H = ħH/(2πk_B)
   and boundary entropy:
      S_H = (k_B/ℓ_κ^2) * A_H * κ
   Under Planck saturation:
      S_H = k_B * A_H / ℓ_P^2   (independent of κ).

3) Outputs Λ in three conventions:
   A) ROT-raw:     E = 1.0 * T*S     → Λ_time = 12 H^2
   B) GR-matched:  E = 0.25 * T*S    → Λ_time = 3 H^2 (pure de Sitter)
   C) ΛCDM ref:    Λ = 3 Ω_Λ H0^2 / c^2

Important note (honest physics bookkeeping)
-------------------------------------------
- If you supply (H0, Ω_Λ), then Λ in ΛCDM is fixed by definition:
      Λ = 3 Ω_Λ H0^2 / c^2
  This script reproduces that exactly in section (C), and shows the ROT mapping
  (via the horizon-energy convention factor) aligns with it in section (B).

Usage (PowerShell)
------------------
  python rot_lambda_predictor_final.py
  python rot_lambda_predictor_final.py --H0 2.2e-18 --Omega_Lambda 0.69
  python rot_lambda_predictor_final.py --H0 2.2e-18 --Omega_Lambda 1.0 --show_details
  python rot_lambda_predictor_final.py --H0 2.2e-18 --Omega_Lambda 0.69 --energy_prefactor_raw 1.0 --energy_prefactor_gr 0.25

Outputs
-------
- Λ_time   in 1/s^2
- Λ_length in 1/m^2

This file is designed to be uploaded directly to GitHub as-is.

MIT License
-----------
Copyright (c) 2026
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass


# -----------------------------
# Physical constants (SI)
# -----------------------------
c = 299_792_458.0
G = 6.67430e-11
hbar = 1.054_571_817e-34
kB = 1.380_649e-23
pi = math.pi

ell_P = math.sqrt(hbar * G / (c**3))
I_P = 1.0 / (ell_P**2)  # Planck saturation target for I := κ/ℓ_κ^2


# -----------------------------
# Horizon quantities
# -----------------------------
def horizon_R(H: float) -> float:
    H = max(float(H), 1e-80)
    return c / H


def horizon_A(H: float) -> float:
    R = horizon_R(H)
    return 4.0 * pi * R * R


def horizon_V(H: float) -> float:
    R = horizon_R(H)
    return (4.0 / 3.0) * pi * (R**3)


def horizon_T(H: float) -> float:
    # Gibbons–Hawking temperature for cosmological horizon
    H = max(float(H), 1e-80)
    return (hbar * H) / (2.0 * pi * kB)


# -----------------------------
# ROT boundary entropy with Planck-saturation closure
# -----------------------------
def rot_entropy_boundary_planck_saturated(A_H: float) -> float:
    """
    Boundary entropy:
        S_H = (kB / ℓ_κ^2) * A_H * κ
    Under Planck saturation closure:
        κ/ℓ_κ^2 = 1/ℓ_P^2
    so
        S_H = kB * A_H / ℓ_P^2 = kB * A_H * I_P
    """
    return kB * float(A_H) * I_P


# -----------------------------
# Energy density and Λ mapping
# -----------------------------
def rho_lambda_from_horizon(H: float, energy_prefactor: float) -> float:
    """
    ROT + horizon-thermo mapping:
        E = energy_prefactor * T_H * S_H
        ρ = E / V_H

    energy_prefactor controls the convention:
      - 1.0  -> "ROT-raw" (gives Λ_time = 12 H^2 under Planck saturation)
      - 0.25 -> "GR-matched" (gives Λ_time = 3 H^2 in pure de Sitter)
    """
    H = max(float(H), 1e-80)
    A = horizon_A(H)
    V = horizon_V(H)
    T = horizon_T(H)
    S = rot_entropy_boundary_planck_saturated(A)

    E = float(energy_prefactor) * T * S
    return E / max(V, 1e-300)  # J/m^3


def Lambda_from_rho(rho_energy: float) -> tuple[float, float]:
    """
    Convert vacuum energy density (J/m^3) to Λ.

    Using your convention from earlier scripts:
        Λ_time = (8πG/c^2) * ρ
        Λ_length = Λ_time / c^2
    """
    rho_energy = float(rho_energy)
    Lambda_time = (8.0 * pi * G / (c**2)) * rho_energy
    Lambda_length = Lambda_time / (c**2)
    return Lambda_time, Lambda_length


def Lambda_lcdm_from_H0_OmegaL(H0: float, Omega_Lambda: float) -> tuple[float, float]:
    """
    ΛCDM / GR definition:
        Λ = 3 Ω_Λ H0^2 / c^2   (length units)
    Returns (Λ_time, Λ_length).
    """
    H0 = float(H0)
    Omega_Lambda = float(Omega_Lambda)
    Lambda_length = 3.0 * Omega_Lambda * (H0**2) / (c**2)
    Lambda_time = Lambda_length * (c**2)
    return Lambda_time, Lambda_length


@dataclass(frozen=True)
class PredictorConfig:
    H0: float
    Omega_Lambda: float
    energy_prefactor_raw: float
    energy_prefactor_gr: float
    show_details: bool


def pretty(x: float) -> str:
    return f"{float(x):.12e}"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="ROT Λ predictor (final): boundary Planck-saturation closure + horizon thermodynamic mapping."
    )
    ap.add_argument("--H0", type=float, default=2.2e-18, help="Hubble rate today (1/s).")
    ap.add_argument("--Omega_Lambda", type=float, default=0.69, help="Dark energy fraction Ω_Λ (dimensionless).")
    ap.add_argument("--energy_prefactor_raw", type=float, default=1.0, help="Prefactor for ROT-raw energy E = f*T*S.")
    ap.add_argument("--energy_prefactor_gr", type=float, default=0.25, help="Prefactor for GR-matched energy E = f*T*S.")
    ap.add_argument("--show_details", action="store_true", help="Print intermediate horizon quantities.")
    args = ap.parse_args()

    cfg = PredictorConfig(
        H0=args.H0,
        Omega_Lambda=args.Omega_Lambda,
        energy_prefactor_raw=args.energy_prefactor_raw,
        energy_prefactor_gr=args.energy_prefactor_gr,
        show_details=bool(args.show_details),
    )

    # Effective “Λ-sector” H in ΛCDM: H_Λ^2 = Ω_Λ H0^2
    H_L = cfg.H0 * math.sqrt(max(cfg.Omega_Lambda, 0.0))

    # A) ROT-raw mapping
    rho_raw = rho_lambda_from_horizon(H_L, energy_prefactor=cfg.energy_prefactor_raw)
    LamT_raw, LamL_raw = Lambda_from_rho(rho_raw)

    # B) GR-matched mapping
    rho_gr = rho_lambda_from_horizon(H_L, energy_prefactor=cfg.energy_prefactor_gr)
    LamT_gr, LamL_gr = Lambda_from_rho(rho_gr)

    # C) ΛCDM reference
    LamT_lcdm, LamL_lcdm = Lambda_lcdm_from_H0_OmegaL(cfg.H0, cfg.Omega_Lambda)

    print("\nROT Λ predictor (final, GitHub-ready)")
    print("--------------------------------------------------------")
    print(f"H0                 = {pretty(cfg.H0)}  1/s")
    print(f"Omega_Lambda       = {cfg.Omega_Lambda:.6f}")
    print(f"H_L = H0*sqrt(ΩΛ)   = {pretty(H_L)}  1/s")
    print(f"Planck length ℓ_P   = {pretty(ell_P)}  m")
    print("")

    print(f"A) ROT-raw mapping:        E = {cfg.energy_prefactor_raw:g} * (T*S)")
    print(f"   Lambda_time   = {pretty(LamT_raw)}  1/s^2")
    print(f"   Lambda_length = {pretty(LamL_raw)}  1/m^2")
    print("")

    print(f"B) GR-matched mapping:     E = {cfg.energy_prefactor_gr:g} * (T*S)")
    print(f"   Lambda_time   = {pretty(LamT_gr)}  1/s^2")
    print(f"   Lambda_length = {pretty(LamL_gr)}  1/m^2")
    print("")

    print("C) ΛCDM reference (GR):    Λ = 3 ΩΛ H0^2 / c^2")
    print(f"   Lambda_time   = {pretty(LamT_lcdm)}  1/s^2")
    print(f"   Lambda_length = {pretty(LamL_lcdm)}  1/m^2")
    print("")

    if cfg.show_details:
        A = horizon_A(H_L)
        V = horizon_V(H_L)
        T = horizon_T(H_L)
        S = rot_entropy_boundary_planck_saturated(A)
        R = horizon_R(H_L)

        print("Details (evaluated at H_L):")
        print(f"   R_H = c/H_L    = {pretty(R)}  m")
        print(f"   A_H            = {pretty(A)}  m^2")
        print(f"   V_H            = {pretty(V)}  m^3")
        print(f"   T_H            = {pretty(T)}  K")
        print(f"   S_H (sat)      = {pretty(S)}  J/K")
        print(f"   rho_raw        = {pretty(rho_raw)}  J/m^3")
        print(f"   rho_gr         = {pretty(rho_gr)}  J/m^3")

    # Tiny self-check: GR-matched should equal ΛCDM when prefactor_gr = 0.25
    # We do not hard-fail; we just report mismatch if user changed prefactor.
    if abs(cfg.energy_prefactor_gr - 0.25) < 1e-15:
        rel_err = abs(LamL_gr - LamL_lcdm) / max(abs(LamL_lcdm), 1e-300)
        print(f"Check: (B) vs (C) relative difference = {rel_err:.3e}")
        print("")


if __name__ == "__main__":
    main()
