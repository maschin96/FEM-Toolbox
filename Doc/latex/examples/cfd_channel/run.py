"""Run the reduced experimental CFD reference, without a RANS claim."""

import numpy as np

from fem_toolbox.cfd import solve_channel

result = solve_channel(
    faces_m=np.linspace(0, 1, 65),
    viscosity_Pa_s=2,
    pressure_gradient_Pa_m=-4,
)
exact = result.centers_m * (1 - result.centers_m)
error = float(np.max(np.abs(result.velocity_m_s - exact)))
if error >= 0.00025 or abs(result.power_balance_W_m2) >= 1e-12:
    raise RuntimeError("channel reference gate failed")
print(f"velocity error: {error:.6g} m/s")
print(f"flow per span: {result.flow_rate_m2_s:.9g} m²/s")
print(f"momentum balance: {result.momentum_balance_Pa:.6g} Pa")
print(f"power balance: {result.power_balance_W_m2:.6g} W/m²")
