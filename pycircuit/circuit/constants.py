# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""
Physical constants
"""

## ⚠ SI 2019 exact values (k and q are DEFINED since the 2019 revision).
## They were 1.38e-23 (4.7e-4 low) and 1.602e-19 (1.1e-4 low) until
## 2026-09-05, which an external reference-simulator cross-check had to carry as
## PARAMETERS on both sides of every noise test to keep the tools from
## disagreeing over a constant.  `eps0` is CODATA 2018 (it was 8.8542e-12,
## 1.4e-6 shy).  ⚠ The COMPACT MODELS keep their own permittivity
## literals (`elements_hdl`: `8.854187817e-12`) -- a model calibration
## choice, as PSP carries its own constants in `psp_scaling` -- so this
## changes `eps0` for callers that read it, not the device physics.
## ⚠ SI-2019 EXACT, AND SPICE-DERIVED TOOLS ARE NOT.  SPICE and the
## commercial simulators built on it carry CODATA-1986,
## `k_B = 1.3806226e-23`.  The ratio is 1.0000191218, so EVERY thermal-noise
## PSD compared against such a tool sits +19.12 ppm high, constant and
## expected; small-signal gains and normalised noise integrations are
## unaffected.  ⚠ Do not write a noise gate against another simulator's PSD
## tighter than 19.12 ppm -- it would fail for a reason that is not ours, and
## do not "fix" this by moving the constant backwards.
kboltzmann=1.380649e-23   # Boltzmann's constant, exact (SI 2019)
eps0 = 8.8541878128e-12   # Vacuum permittivity (CODATA 2018)
epsRSi = 11.7             # Relative permittivity of Si
epsRSiO2 = 3.9            # Relative permittivity of SiO2 
qelectron=1.602176634e-19 # Elementary charge, exact (SI 2019)
