# STFRFT — Short-Time Fractional Fourier Transform

This repository provides a Python implementation of the **Short-Time Fractional Fourier Transform (STFRFT)** based on the discrete formulation described in:

> **Zhichun Zhao and Gang Li**,  
> *Synchrosqueezing-based Short-Time Fractional Fourier Transform*,  
> *IEEE Transactions on Signal Processing*, vol. **71**, pp. 279–294, 2023.

---

## 🧠 Description

The STFRFT is a time–fractional frequency representation that generalizes the classical Short-Time Fourier Transform (STFT).  
It allows analysis of **chirp-like** and **non-stationary signals** with improved energy concentration by introducing a fractional rotation in the time–frequency plane.

This implementation provides:
- A **class-based structure** (`STFRFT`) for clean and reusable code.
- Functions for computing:
  - The **fractional kernel** `K_alpha`.
  - The **STFRFT matrix** `V_alpha`.
- Parameters for customizing:
  - Fractional order `alpha`
  - Window width `sigma`
  - Number of fractional bins `nb_q`
  - Sampling interval `delta_t`

---

## 🧩 Example Usage

```python
import numpy as np
from stfrft import STFRFT

# Example signal
fs = 1000
T = 1
t = np.linspace(0, T, int(fs*T))
s = np.exp(1j * 2 * np.pi * (100 * t + 50 * t**2))  # linear chirp

# Initialize STFRFT
st = STFRFT(s, nb_q=256, alpha=np.pi/6, delta_t=1/fs, sigma=0.5)

# Compute kernel and transform
K = st.K_alpha()
V = st.V_alpha()
