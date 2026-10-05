"""Sparse affine field maps and conservative Oseen/SIMPLE block assembly."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import (  # type: ignore[import-untyped]
    block_diag,
    coo_matrix,
    csr_matrix,
    diags,
    eye,
    kron,
)

from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.quad_fv import QuadFVGeometry, _error


@dataclass
class Map:
    """Private affine field map; caller supplies already validated arrays."""

    matrix: csr_matrix
    constant: NDArray[np.float64]

    def apply(self, value: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.asarray(self.matrix @ value + self.constant)


def sparse(
    rows: list[int], cols: list[int], data: list[float], shape: tuple[int, int]
) -> csr_matrix:
    return coo_matrix((data, (rows, cols)), shape=shape).tocsr()


def affine_product(matrix: csr_matrix, mapping: Map) -> Map:
    return Map((matrix @ mapping.matrix).tocsr(), np.asarray(matrix @ mapping.constant))


class QuadSystem:
    """Internal maps preserve common face states, full stress and BC constraints."""

    def __init__(
        self, domain: FluidDomain, mu: NDArray[np.float64], mub: NDArray[np.float64]
    ) -> None:
        self.domain = domain
        self.fv = QuadFVGeometry.from_domain(domain)
        self.mesh = domain.mesh
        self.n, self.nf = len(self.mesh.quad4), len(self.mesh.faces)
        self.kind = np.full(self.nf, "internal", dtype="U20")
        self.ub = np.zeros((self.nf, 2))
        self.pressure_rows: list[int] = []
        self.walls: list[int] = []
        self.reference = 0.0
        for name, bc in domain.boundaries.items():
            ids = domain.mesh.boundary_sets[name]
            faces = domain.mesh.boundary_face_indices[ids]
            if bc["type"] == "periodic":
                continue
            self.kind[faces] = bc["type"]
            if bc["type"] == "inlet":
                self.ub[faces] = bc["velocity_m_s"]
            elif bc["type"] == "pressure_outlet":
                self.pressure_rows.extend(int(i) for i in ids)
                self.reference = float(bc["pressure_Pa"])
            elif bc["type"] == "wall":
                self.walls.extend(int(i) for i in ids)
        r, c, d = [], [], []
        for f in self.fv.active_faces:
            r.append(int(self.mesh.face_owners[f]))
            c.append(int(f))
            d.append(1.0)
            q = int(self.fv.neighbours[f])
            if q >= 0:
                r.append(q)
                c.append(int(f))
                d.append(-1.0)
        self.balance = sparse(r, c, d, (self.n, self.nf))
        self.vector_balance = kron(self.balance, eye(2), format="csr")
        self.poly_scale = np.sqrt(self.mesh.cell_areas_m2)
        self.gu, self.poly_u = self._polynomial_reconstruction(velocity=True)
        self.gp, self.poly_p = self._polynomial_reconstruction(velocity=False)
        self.uf = self._polynomial_faces(self.poly_u, velocity=True)
        self.pf = self._polynomial_faces(self.poly_p, velocity=False)
        self.guf, self.guf_average = self._face_gradient(
            self.gu, self.uf, velocity=True
        )
        self.guf = self._polynomial_face_gradient(self.guf)
        self.gpf, self.gpf_average = self._face_gradient(
            self.gp, self.pf, velocity=False
        )
        self.mu_cell = mu
        self.mu_face = self.fv.face_values(mu, mub)
        if np.any(self.mu_face <= 0):
            _error("viscosity_Pa_s", "nonpositive reconstructed face viscosity")
        r, c, d = [], [], []
        for f, (nx, ny) in enumerate(self.mesh.face_normals):
            for component, entries in enumerate(
                ([(0, 2 * nx), (1, ny), (2, ny)], [(1, nx), (2, nx), (3, 2 * ny)])
            ):
                for axis, value in entries:
                    r.append(2 * f + component)
                    c.append(4 * f + axis)
                    d.append(float(self.mu_face[f] * value))
        self.viscous_traction = affine_product(
            sparse(r, c, d, (2 * self.nf, 4 * self.nf)), self.guf
        )
        area = diags(np.repeat(self.fv.face_areas_m2, 2))
        self.viscous_force = affine_product(area, self.viscous_traction)
        self.diffusion = affine_product(-self.vector_balance, self.viscous_force)
        r, c, d = [], [], []
        for f, normal in enumerate(self.mesh.face_normals):
            for k in range(2):
                r.append(2 * f + k)
                c.append(f)
                d.append(float(normal[k] * self.fv.face_areas_m2[f]))
        self.pressure_force = affine_product(
            sparse(r, c, d, (2 * self.nf, self.nf)), self.pf
        )
        self.pressure_momentum = affine_product(
            self.vector_balance, self.pressure_force
        )
        self.normal = sparse(
            [f for f in range(self.nf) for _ in range(2)],
            [2 * f + k for f in range(self.nf) for k in range(2)],
            self.mesh.face_normals.ravel().tolist(),
            (self.nf, 2 * self.nf),
        )
        self.velocity_flux = affine_product(
            diags(domain.fluid.density_kg_m3 * self.fv.face_areas_m2) @ self.normal,
            self.uf,
        )
        self.mass_velocity = affine_product(self.balance, self.velocity_flux)

    def _polynomial_reconstruction(self, *, velocity: bool) -> tuple[Map, Map]:
        """Consistent cubic cell and face reconstruction with BC constraints.

        Local P3 fit anchored at the cell value, using two interior or three
        boundary face-neighbour rings. No fallback on deficient stencils.
        Periodic images use translated coordinates and the same cell values.
        """
        components = 2 if velocity else 1
        adjacent: list[list[tuple[int, NDArray[np.float64]]]] = [
            [] for _ in range(self.n)
        ]
        outer: list[list[int]] = [[] for _ in range(self.n)]
        for f in self.fv.active_faces:
            p, q = int(self.mesh.face_owners[f]), int(self.fv.neighbours[f])
            if q >= 0:
                adjacent[p].append((q, self.fv._direction[f]))
                adjacent[q].append((p, -self.fv._direction[f]))
            else:
                outer[p].append(int(f))
        r: list[int] = []
        c: list[int] = []
        d: list[float] = []
        constant = np.zeros(9 * components * self.n)
        for cell in range(self.n):
            scale = self.poly_scale[cell]
            locations = [(cell, np.zeros(2))]
            frontier = [(cell, np.zeros(2))]
            seen = {(cell, 0.0, 0.0)}
            for _ in range(3 if outer[cell] else 2):
                following = []
                for p, offset in frontier:
                    for q, delta in adjacent[p]:
                        location = offset + delta
                        key = (q, *np.round(location / scale, 12))
                        if key not in seen:
                            seen.add(key)
                            locations.append((q, location))
                            following.append((q, location))
                frontier = following
            rows: list[NDArray[np.float64]] = []
            rhs_entries: list[list[tuple[int, float]]] = []
            rhs_constants: list[float] = []

            def sample_basis(
                offset: NDArray[np.float64], scale: float = scale
            ) -> NDArray[np.float64]:
                x, y = offset / scale
                return np.array(
                    [
                        x,
                        y,
                        x * x / 2,
                        x * y,
                        y * y / 2,
                        x**3 / 6,
                        x * x * y / 2,
                        x * y * y / 2,
                        y**3 / 6,
                    ]
                )

            def sample(
                offset: NDArray[np.float64],
                target: int | None,
                value: NDArray[np.float64],
                scale: float = scale,
                cell: int = cell,
                rows: list[NDArray[np.float64]] = rows,
                rhs_entries: list[list[tuple[int, float]]] = rhs_entries,
                rhs_constants: list[float] = rhs_constants,
            ) -> None:
                weight = 1 / max(float(np.linalg.norm(offset / scale)), 0.25)
                basis = sample_basis(offset) * weight
                for k in range(components):
                    row = np.zeros(9 * components)
                    row[9 * k : 9 * k + 9] = basis
                    entries = [(components * cell + k, -weight)]
                    if target is not None:
                        entries.append((components * target + k, weight))
                    rows.append(row)
                    rhs_entries.append(entries)
                    rhs_constants.append(float(value[k] * weight))

            for q, offset in locations:
                if np.linalg.norm(offset) > 1e-12 * scale:
                    sample(offset, q, np.zeros(components))
                for f in outer[q]:
                    location = (
                        offset
                        + self.mesh.face_centres_m[f]
                        - self.mesh.cell_centres_m[q]
                    )
                    kind = self.kind[f]
                    if (velocity and kind in ("wall", "inlet")) or (
                        not velocity and kind == "pressure_outlet"
                    ):
                        sample(location, None, self.ub[f] if velocity else np.zeros(1))
                    elif velocity:
                        normal = self.mesh.face_normals[f]
                        tangent = np.array([-normal[1], normal[0]])
                        x, y = location / scale
                        derivative = np.array(
                            [
                                normal[0],
                                normal[1],
                                normal[0] * x,
                                normal[0] * y + normal[1] * x,
                                normal[1] * y,
                                normal[0] * x * x / 2,
                                normal[0] * x * y + normal[1] * x * x / 2,
                                normal[0] * y * y / 2 + normal[1] * x * y,
                                normal[1] * y * y / 2,
                            ]
                        )
                        if kind == "pressure_outlet":
                            for k in range(2):
                                row = np.zeros(18)
                                row[9 * k : 9 * k + 9] = derivative
                                rows.append(row)
                                rhs_entries.append([])
                                rhs_constants.append(0.0)
                        elif kind == "symmetry":
                            rows.append(
                                np.outer(normal, sample_basis(location)).ravel()
                            )
                            rhs_entries.append(
                                [(2 * cell + k, float(-normal[k])) for k in range(2)]
                            )
                            rhs_constants.append(0.0)
                            rows.append(np.outer(tangent, derivative).ravel())
                            rhs_entries.append([])
                            rhs_constants.append(0.0)
            design = np.array(rows)
            if (
                np.linalg.matrix_rank(design) < 9 * components
                or np.linalg.cond(design) > 1e12
            ):
                _error(
                    f"boundary.cell[{cell}]", "unresolved third-degree boundary stencil"
                )
            fit = np.linalg.pinv(design)
            for obs, entries in enumerate(rhs_entries):
                for k in range(9 * components):
                    output = 9 * components * cell + k
                    coefficient = float(fit[k, obs])
                    constant[output] += coefficient * rhs_constants[obs]
                    for col, weight in entries:
                        r.append(output)
                        c.append(col)
                        d.append(coefficient * weight)
        polynomial = Map(
            sparse(r, c, d, (9 * components * self.n, components * self.n)), constant
        )
        r, c, d = [], [], []
        for cell in range(self.n):
            for k in range(components):
                for axis in range(2):
                    gradient_row = 2 * components * cell + 2 * k + axis
                    r.append(gradient_row)
                    c.append(9 * components * cell + 9 * k + axis)
                    d.append(1 / self.poly_scale[cell])
        select = sparse(r, c, d, (2 * components * self.n, 9 * components * self.n))
        return Map(
            (select @ polynomial.matrix).tocsr(),
            np.asarray(select @ polynomial.constant),
        ), polynomial

    def _polynomial_faces(self, polynomial: Map, *, velocity: bool) -> Map:
        components = 2 if velocity else 1
        r: list[int] = []
        c: list[int] = []
        d: list[float] = []
        gr: list[int] = []
        gc: list[int] = []
        gd: list[float] = []
        constant = self.ub.ravel().copy() if velocity else np.zeros(self.nf)
        for f in range(self.nf):
            q = int(self.fv.neighbours[f])
            if q < 0 and (
                (velocity and self.kind[f] in ("wall", "inlet"))
                or (not velocity and self.kind[f] == "pressure_outlet")
            ):
                continue
            p = int(self.mesh.face_owners[f])
            offset = self.mesh.face_centres_m[f] - self.mesh.cell_centres_m[p]
            fraction = self.fv._fraction[f]
            sides = (
                [
                    (p, 1 - fraction, offset),
                    (q, fraction, offset - self.fv._direction[f]),
                ]
                if q >= 0
                else [(p, 1.0, offset)]
            )
            projection = np.eye(components)
            if velocity and self.kind[f] == "symmetry":
                projection -= np.outer(
                    self.mesh.face_normals[f], self.mesh.face_normals[f]
                )
            for k in range(components):
                constant[components * f + k] = 0
                for cell, weight, delta in sides:
                    x, y = delta / self.poly_scale[cell]
                    basis = [
                        x,
                        y,
                        x * x / 2,
                        x * y,
                        y * y / 2,
                        x**3 / 6,
                        x * x * y / 2,
                        x * y * y / 2,
                        y**3 / 6,
                    ]
                    for j in range(components):
                        coefficient = float(weight * projection[k, j])
                        r.append(components * f + k)
                        c.append(components * cell + j)
                        d.append(coefficient)
                        for v in range(9):
                            gr.append(components * f + k)
                            gc.append(9 * components * cell + 9 * j + v)
                            gd.append(coefficient * basis[v])
        base = sparse(r, c, d, (components * self.nf, components * self.n))
        extrapolate = sparse(
            gr, gc, gd, (components * self.nf, 9 * components * self.n)
        )
        return Map(
            (base + extrapolate @ polynomial.matrix).tocsr(),
            constant + np.asarray(extrapolate @ polynomial.constant),
        )

    def _polynomial_face_gradient(self, mapping: Map) -> Map:
        r: list[int] = []
        c: list[int] = []
        d: list[float] = []
        keep = np.ones(4 * self.nf)
        for f in range(self.nf):
            if self.fv.neighbours[f] >= 0:
                continue
            p = int(self.mesh.face_owners[f])
            scale = self.poly_scale[p]
            x, y = (self.mesh.face_centres_m[f] - self.mesh.cell_centres_m[p]) / scale
            evaluation = (
                np.array(
                    [
                        [1, 0, x, y, 0, x * x / 2, x * y, y * y / 2, 0],
                        [0, 1, 0, x, y, 0, x * x / 2, x * y, y * y / 2],
                    ]
                )
                / scale
            )
            normal = self.mesh.face_normals[f]
            nn = np.outer(normal, normal)
            tt = np.eye(2) - nn
            transform: NDArray[np.float64] = np.eye(4)
            if self.kind[f] == "wall":
                transform = np.asarray(np.kron(np.eye(2), nn), dtype=np.float64)
            elif self.kind[f] == "pressure_outlet":
                transform = np.asarray(np.kron(np.eye(2), tt), dtype=np.float64)
            elif self.kind[f] == "symmetry":
                transform = np.kron(nn, nn) + np.kron(tt, tt)
            values = transform @ np.kron(np.eye(2), evaluation)
            for k in range(4):
                keep[4 * f + k] = 0
                for j in range(18):
                    r.append(4 * f + k)
                    c.append(18 * p + j)
                    d.append(float(values[k, j]))
        self.physical_gradient_mask = keep == 0
        self.physical_gradient = affine_product(
            sparse(r, c, d, (4 * self.nf, 18 * self.n)), self.poly_u
        )
        return Map(
            (diags(keep) @ mapping.matrix + self.physical_gradient.matrix).tocsr(),
            keep * mapping.constant + self.physical_gradient.constant,
        )

    def _face_gradient(
        self, gradient: Map, face: Map, *, velocity: bool
    ) -> tuple[Map, Map]:
        components = 2 if velocity else 1
        r: list[int] = []
        c: list[int] = []
        d: list[float] = []
        dr: list[int] = []
        dc: list[int] = []
        dd: list[float] = []
        br: list[int] = []
        bc: list[int] = []
        bd: list[float] = []
        tensors = []
        tr: list[int] = []
        tc: list[int] = []
        td: list[float] = []
        for f in range(self.nf):
            p, q = int(self.mesh.face_owners[f]), int(self.fv.neighbours[f])
            normal = self.mesh.face_normals[f]
            distance = self.fv._distance[f]
            fraction = self.fv._fraction[f]
            sides = [(p, 1 - fraction), (q, fraction)] if q >= 0 else [(p, 1.0)]
            for k in range(components):
                for axis in range(2):
                    row = 2 * components * f + 2 * k + axis
                    for cell, weight in sides:
                        if velocity:
                            offset = (
                                self.mesh.face_centres_m[f]
                                - self.mesh.cell_centres_m[p]
                            )
                            if cell == q:
                                offset = offset - self.fv._direction[f]
                            x, y = offset / self.poly_scale[cell]
                            evaluation = (
                                [1, 0, x, y, 0, x * x / 2, x * y, y * y / 2, 0]
                                if axis == 0
                                else [0, 1, 0, x, y, 0, x * x / 2, x * y, y * y / 2]
                            )
                            for j, value in enumerate(evaluation):
                                r.append(row)
                                c.append(9 * components * cell + 9 * k + j)
                                d.append(float(weight * value / self.poly_scale[cell]))
                            if q >= 0:
                                dx, dy = self.fv._direction[f] / self.poly_scale[cell]
                                cubic = [
                                    dx**3,
                                    3 * dx * dx * dy,
                                    3 * dx * dy * dy,
                                    dy**3,
                                ]
                                for j, value in enumerate(cubic):
                                    tr.append(row)
                                    tc.append(9 * components * cell + 9 * k + 5 + j)
                                    td.append(
                                        float(
                                            -weight
                                            * normal[axis]
                                            * value
                                            / (24 * distance)
                                        )
                                    )
                        else:
                            r.append(row)
                            c.append(2 * components * cell + 2 * k + axis)
                            d.append(float(weight))
                    dr.append(row)
                    dc.append(components * p + k)
                    dd.append(float(-normal[axis] / distance))
                    if q >= 0:
                        dr.append(row)
                        dc.append(components * q + k)
                        dd.append(float(normal[axis] / distance))
                    else:
                        br.append(row)
                        bc.append(components * f + k)
                        bd.append(float(normal[axis] / distance))
            tensors.append(
                np.kron(
                    np.eye(components),
                    np.eye(2) - np.outer(normal, self.fv._direction[f]) / distance,
                )
            )
        average = affine_product(
            sparse(
                r,
                c,
                d,
                (
                    2 * components * self.nf,
                    (9 if velocity else 2) * components * self.n,
                ),
            ),
            self.poly_u if velocity else gradient,
        )
        transform = block_diag([csr_matrix(t) for t in tensors], format="csr")
        boundary = sparse(br, bc, bd, (2 * components * self.nf, components * self.nf))
        direct = sparse(dr, dc, dd, (2 * components * self.nf, components * self.n))
        correction = (
            affine_product(
                sparse(tr, tc, td, (2 * components * self.nf, 9 * components * self.n)),
                self.poly_u,
            )
            if velocity
            else Map(csr_matrix(direct.shape), np.zeros(2 * components * self.nf))
        )
        return Map(
            (
                transform @ average.matrix
                + direct
                + boundary @ face.matrix
                + correction.matrix
            ).tocsr(),
            np.asarray(
                transform @ average.constant
                + boundary @ face.constant
                + correction.constant
            ),
        ), average

    def alias_flux(self, value: NDArray[np.float64]) -> NDArray[np.float64]:
        result = value.copy()
        for f, q in self.fv.periodic_partners.items():
            if f < q:
                result[q] = -result[f]
        return result

    def pressure_flux(self, diagonal: NDArray[np.float64]) -> Map:
        df, _, rho = self.fv._pressure_data(
            diagonal, self.domain.fluid.density_kg_m3, self.pressure_rows
        )
        active = (self.fv.neighbours >= 0) | (self.kind == "pressure_outlet")
        scale = diags(rho * self.fv.face_areas_m2 * df * active)
        return Map(
            (scale @ self.normal @ (self.gpf_average.matrix - self.gpf.matrix)).tocsr(),
            np.asarray(
                scale @ self.normal @ (self.gpf_average.constant - self.gpf.constant)
            ),
        )

    def transport(
        self, flux: NDArray[np.float64], u: NDArray[np.float64], scheme: str
    ) -> tuple[Map, NDArray[np.float64]]:
        if scheme == "central":
            peclet = (
                np.abs(flux)
                * self.fv._distance
                / (self.mu_face * self.fv.face_areas_m2)
            )
            if np.any(peclet >= 2):
                f = int(np.argmax(peclet))
                _error(f"face[{f}].peclet", "central scheme requires face Peclet <2")
            mapping = self.uf
            return affine_product(
                self.vector_balance @ diags(np.repeat(flux, 2)), mapping
            ), mapping.apply(u)
        r: list[int] = []
        c: list[int] = []
        d: list[float] = []
        constant = np.zeros(2 * self.nf)
        gradient = self.gu.apply(u).reshape(self.n, 2, 2)
        cell_u = u.reshape(self.n, 2)
        for f in range(self.nf):
            p, q = int(self.mesh.face_owners[f]), int(self.fv.neighbours[f])
            if q < 0:
                continue
            donor = p if flux[f] >= 0 else q
            for k in range(2):
                r.append(2 * f + k)
                c.append(2 * donor + k)
                d.append(1.0)
            if scheme == "linear-upwind":
                offset = self.mesh.face_centres_m[f] - self.mesh.cell_centres_m[p]
                if donor == q:
                    offset -= self.fv._direction[f]
                extrapolated = cell_u[donor] + gradient[donor] @ offset
                limited = np.clip(
                    extrapolated,
                    np.minimum(cell_u[p], cell_u[q]),
                    np.maximum(cell_u[p], cell_u[q]),
                )
                constant[2 * f : 2 * f + 2] = limited - cell_u[donor]
        donor_matrix = sparse(r, c, d, (2 * self.nf, 2 * self.n))
        exterior = np.repeat(self.fv.neighbours < 0, 2)
        select = diags(exterior.astype(float))
        mapping = Map(
            (donor_matrix + select @ self.uf.matrix).tocsr(),
            constant + np.asarray(select @ self.uf.constant),
        )
        cell_map = affine_product(
            self.vector_balance @ diags(np.repeat(flux, 2)), mapping
        )
        return cell_map, mapping.apply(u)
