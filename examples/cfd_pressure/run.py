"""Independent smooth Stokes reference on four uniform MAC grids."""

import numpy as np

from fem_toolbox.cfd import solve_stokes

pi = np.pi
velocity_errors = []
pressure_errors = []
for n in (8, 16, 32, 64):
    yu, xu = np.meshgrid((np.arange(n) + 0.5) / n, np.arange(1, n) / n, indexing="ij")
    yv, xv = np.meshgrid(np.arange(1, n) / n, (np.arange(n) + 0.5) / n, indexing="ij")
    fx = -4 * pi**3 * (2 * np.cos(2 * pi * xu) - 1) * np.sin(
        2 * pi * yu
    ) + 2 * pi * np.cos(2 * pi * xu) * np.sin(2 * pi * yu)
    fy = -4 * pi**3 * np.sin(2 * pi * xv) * (
        1 - 2 * np.cos(2 * pi * yv)
    ) + 2 * pi * np.sin(2 * pi * xv) * np.cos(2 * pi * yv)
    result = solve_stokes(
        cells=(n, n), size_m=(1, 1), viscosity_Pa_s=1, force_x_N_m3=fx, force_y_N_m3=fy
    )
    u = pi * (1 - np.cos(2 * pi * xu)) * np.sin(2 * pi * yu)
    v = -pi * np.sin(2 * pi * xv) * (1 - np.cos(2 * pi * yv))
    yp, xp = np.meshgrid(
        (np.arange(n) + 0.5) / n, (np.arange(n) + 0.5) / n, indexing="ij"
    )
    pressure = np.sin(2 * pi * xp) * np.sin(2 * pi * yp)
    pressure -= pressure[0, 0]
    velocity_errors.append(
        float(
            np.sqrt(
                np.mean((result.u_m_s[:, 1:-1] - u) ** 2)
                + np.mean((result.v_m_s[1:-1, :] - v) ** 2)
            )
        )
    )
    pressure_errors.append(
        float(np.sqrt(np.mean((result.pressure_Pa - pressure) ** 2)))
    )
    print(
        f"n={n}: velocity RMS={velocity_errors[-1]:.6g} m/s, pressure RMS={pressure_errors[-1]:.6g} Pa, max div={np.max(np.abs(result.divergence_s_inv)):.6g} 1/s"
    )
u_orders = np.log2(np.array(velocity_errors[:-1]) / velocity_errors[1:])
p_orders = np.log2(np.array(pressure_errors[:-1]) / pressure_errors[1:])
if (
    not np.all(u_orders > 1.8)
    or not np.all(p_orders > 1.5)
    or velocity_errors[-1] >= 0.01
    or pressure_errors[-1] >= 0.01
):
    raise RuntimeError("manufactured Stokes reference gate failed")
print(f"velocity orders: {u_orders}")
print(f"pressure orders: {p_orders}")
