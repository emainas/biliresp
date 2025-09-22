#!/usr/bin/env python3
"""
Plot the RESP hyperbolic restraint for one charge:
    R(q) = a * sqrt( (q - q0)**2 + b**2 )
with a=0.0005, b=0.1, q0=0, and save as PNG.
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 5e-4      # restraint strength
b = 0.1       # soft-corner parameter
q0 = 0.0      # target charge (usually 0 in stage-1)

# Domain for q (in units of elementary charge e)
q_min, q_max = -1.5, 1.5
q = np.linspace(q_min, q_max, 1000)

# Hyperbolic restraint
R = a * np.sqrt((q - q0)**2 + b**2)

# Plot
plt.figure(figsize=(6, 4))
plt.plot(q, R, linewidth=2)
plt.xlabel("charge q (e)")
plt.ylabel("RESP restraint R(q) (a.u.)")
plt.title(r"Hyperbolic RESP restraint: $R(q)=a\sqrt{(q-q_0)^2+b^2}$"
          f"\n(a={a}, b={b}, q0={q0})")
plt.grid(True, alpha=0.3)

# Optional guide lines at ±b (soft corner scale)
plt.axvline(+b, linestyle="--", linewidth=1)
plt.axvline(-b, linestyle="--", linewidth=1)

plt.tight_layout()
out = "resp_hyperbolic_restraint.png"
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"saved: {out}")

