"""Four independent refinements; analytic sine derivatives, no flow solve."""

import numpy as np

from fem_toolbox.cfd import QuadFVGeometry


def test_smooth_gradient_four_meshes(quad_mesh, record_property):
    errors = []
    for n in (16, 32, 64, 128):
        mesh = quad_mesh(n, n)
        fv = QuadFVGeometry(mesh, width_m=1)
        x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
        values = np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])
        boundary = np.sin(np.pi * b[:, 0]) * np.sin(np.pi * b[:, 1])
        exact = np.pi * np.column_stack(
            (
                np.cos(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1]),
                np.sin(np.pi * x[:, 0]) * np.cos(np.pi * x[:, 1]),
            )
        )
        error = np.sqrt(np.mean((fv.gradient(values, boundary) - exact) ** 2))
        errors.append(float(error))
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    record_property("rms_gradient_errors", errors)
    record_property("observed_orders", orders.tolist())
    assert errors[-1] < 0.02, errors
    assert np.all(orders > 1.4), (errors, orders)
