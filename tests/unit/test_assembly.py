"""Independent small algebraic references; float tolerances are roundoff only."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.sparse.linalg import spsolve

from fem_toolbox.assembly import MatrixBlock, VectorBlock, assemble
from fem_toolbox.boundary import Dirichlet, eliminate_dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.fields import DofMap, Field


def test_fields_only_allocate_active_nodes_and_components():
    dofs = DofMap(
        (
            Field("u", ("x", "y", "z"), (4, 1), "m"),
            Field("T", ("scalar",), (1, 8, 4), "K"),
        )
    )
    assert dofs.size == 9
    assert dofs.block("T") == slice(6, 9)
    assert_array_equal(dofs.indices("u", [1, 4], ["z", "x"]), [5, 3, 2, 0])
    assert_array_equal(dofs.indices("T", [8]), [7])
    with pytest.raises(ValidationError):
        dofs.indices("T", [2])


@pytest.mark.parametrize(
    "field",
    [
        ("", ("x",), (0,), "m"),
        ("u", (), (0,), "m"),
        ("u", ("x", "x"), (0,), "m"),
        ("u", ("x",), (0, 0), "m"),
        ("u", ("x",), (-1,), "m"),
        ("u", ("x",), (1.5,), "m"),
        ("u", ("x",), (0,), ""),
    ],
)
def test_invalid_fields(field):
    with pytest.raises(ValidationError):
        Field(*field)


def test_duplicate_field_names():
    f = Field("T", ("scalar",), (0,), "K")
    with pytest.raises(ValidationError):
        DofMap((f, f))


@pytest.mark.parametrize("format", ["csr", "csc"])
def test_overlapping_assembly_and_nonzero_dirichlet(format):
    # Two unit springs: u0=2, force at node 2 = 3 -> u=[2,5,8].
    spring = [[1.0, -1.0], [-1.0, 1.0]]
    system = assemble(
        3,
        [MatrixBlock([0, 1], [0, 1], spring), MatrixBlock([1, 2], [1, 2], spring)],
        [VectorBlock([2, 2], [1.0, 2.0])],
        format=format,
    )
    assert system.matrix.format == format
    assert_allclose(
        system.matrix.toarray(), [[1, -1, 0], [-1, 2, -1], [0, -1, 1]], atol=0
    )
    constrained = eliminate_dirichlet(
        system, [Dirichlet([0], [2.0]), Dirichlet([0], [2.0])]
    )
    assert_allclose(
        constrained.matrix.toarray(), [[1, 0, 0], [0, 2, -1], [0, -1, 1]], atol=0
    )
    assert_array_equal(constrained.rhs, [2, 2, 3])
    solution = spsolve(constrained.matrix, constrained.rhs)
    assert_allclose(solution, [2, 5, 8], rtol=0, atol=1e-14)
    assert_allclose(system.residual(solution), [-3, 0, 0], rtol=0, atol=1e-14)
    assert_array_equal(system.rhs, [0, 0, 3])
    # Accessors do not expose the stored original system to mutation.
    system.matrix.data[:] = 0
    system.rhs[:] = 0
    assert_allclose(system.residual(solution), [-3, 0, 0], rtol=0, atol=1e-14)
    with pytest.raises(ValidationError, match="conflict"):
        eliminate_dirichlet(system, [Dirichlet([0, 0], [2, 3])])


def test_rectangular_coupling_blocks_and_empty_system():
    s = assemble(
        3, [MatrixBlock([0, 1], [2], [[4], [5]]), MatrixBlock([2], [0, 1], [[4, 5]])]
    )
    assert_array_equal(s.matrix.toarray(), [[0, 0, 4], [0, 0, 5], [4, 5, 0]])
    assert assemble(0, []).matrix.shape == (0, 0)
    assert_array_equal(eliminate_dirichlet(s, []).matrix.toarray(), s.matrix.toarray())
    all_fixed = eliminate_dirichlet(s, [Dirichlet([0, 1, 2], [1, 2, 3])])
    assert_array_equal(all_fixed.matrix.toarray(), np.eye(3))
    assert_array_equal(all_fixed.rhs, [1, 2, 3])


@pytest.mark.parametrize(
    "block",
    [
        MatrixBlock([3], [0], [[1]]),
        MatrixBlock([-1], [0], [[1]]),
        MatrixBlock([0.5], [0], [[1]]),
        MatrixBlock([0], [0], [[np.nan]]),
        MatrixBlock([0], [0], [1]),
        MatrixBlock([0], [0], [[1j]]),
    ],
)
def test_invalid_matrix_blocks(block):
    with pytest.raises(ValidationError):
        assemble(3, [block])


def test_invalid_vectors_constraints_and_solution():
    with pytest.raises(ValidationError):
        assemble(2, [], [VectorBlock([0], [np.inf])])
    s = assemble(2, [])
    for bc in [Dirichlet([2], [0]), Dirichlet([0], [np.nan]), Dirichlet([0, 1], [0])]:
        with pytest.raises(ValidationError):
            eliminate_dirichlet(s, [bc])
    with pytest.raises(ValidationError):
        s.residual([0])


def test_empty_support_and_lookup_errors():
    dofs = DofMap((Field("T", ("scalar",), (), "K"),))
    assert dofs.size == 0
    assert_array_equal(dofs.indices("T", []), [])
    with pytest.raises(ValidationError):
        dofs.block("missing")
    active = DofMap((Field("u", ("x",), (0,), "m"),))
    with pytest.raises(ValidationError):
        active.indices("u", [0], ["y"])
    empty = eliminate_dirichlet(assemble(0, []), [])
    assert empty.matrix.shape == (0, 0)


def test_repeated_matrix_indices_accumulate():
    system = assemble(1, [MatrixBlock([0, 0], [0, 0], [[1, 2], [3, 4]])])
    assert_array_equal(system.matrix.toarray(), [[10]])


@pytest.mark.parametrize("size", [-1, 1.5, True])
def test_invalid_system_size(size):
    with pytest.raises(ValidationError):
        assemble(size, [])
