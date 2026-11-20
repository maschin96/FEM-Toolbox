"""Independent channel solutions and conservation before full RANS coupling."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError


def solve(**kwargs):
    from fem_toolbox.cfd import solve_channel

    return solve_channel(**kwargs)


def test_poiseuille_conservation():
    r = solve(
        faces_m=np.linspace(0, 1, 65), viscosity_Pa_s=2, pressure_gradient_Pa_m=-4
    )
    y = r.centers_m
    exact = y * (1 - y)
    assert np.max(np.abs(r.velocity_m_s - exact)) < 0.00025
    assert r.shear_Pa[[0, -1]] == pytest.approx([2, -2], rel=0, abs=1e-12)
    assert abs(r.momentum_balance_Pa) < 1e-12
    assert abs(r.power_balance_W_m2) < 1e-12
    assert r.flow_rate_m2_s == pytest.approx(1 / 6, rel=0, abs=0.0002)
    with pytest.raises(ValueError):
        r.velocity_m_s.setflags(write=True)


def test_couette_on_stretched_mesh():
    faces = np.array([0, 0.02, 0.1, 0.4, 0.7, 1])
    r = solve(
        faces_m=faces,
        viscosity_Pa_s=3,
        pressure_gradient_Pa_m=0,
        wall_velocities_m_s=(2, 5),
    )
    np.testing.assert_allclose(r.velocity_m_s, 2 + 3 * r.centers_m, rtol=0, atol=1e-12)
    np.testing.assert_allclose(r.shear_Pa, 9, rtol=0, atol=1e-12)
    assert r.dissipation_W_m2 == pytest.approx(27, rel=0, abs=1e-12)
    assert abs(r.power_balance_W_m2) < 1e-12


def test_reverse_pressure_gradient():
    r = solve(faces_m=[0, 0.5, 1], viscosity_Pa_s=1, pressure_gradient_Pa_m=2)
    assert np.all(r.velocity_m_s < 0)
    assert r.dissipation_W_m2 > 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"faces_m": [0, 0, 1]},
        {"faces_m": [0, 1, float("nan")]},
        {"faces_m": [[0, 1]]},
        {"faces_m": [0]},
        {"viscosity_Pa_s": 0},
        {"viscosity_Pa_s": -1},
        {"viscosity_Pa_s": True},
        {"pressure_gradient_Pa_m": float("inf")},
        {"wall_velocities_m_s": (0,)},
        {"wall_velocities_m_s": (0, float("nan"))},
    ],
)
def test_invalid_data(kwargs):
    options = dict(faces_m=[0, 0.5, 1], viscosity_Pa_s=1, pressure_gradient_Pa_m=-1)
    options.update(kwargs)
    with pytest.raises(ValidationError):
        solve(**options)


def test_four_meshes_second_order():
    errors = []
    for n in [8, 16, 32, 64]:
        r = solve(
            faces_m=np.linspace(0, 1, n + 1),
            viscosity_Pa_s=1,
            pressure_gradient_Pa_m=-2,
        )
        errors.append(np.max(np.abs(r.velocity_m_s - r.centers_m * (1 - r.centers_m))))
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    assert np.all((orders > 1.99) & (orders < 2.01))
