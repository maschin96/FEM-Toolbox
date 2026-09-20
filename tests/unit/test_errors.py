"""Tests for concrete validation failures exposed by the public API."""

from fem_toolbox.errors import ShapeError, UnitError, ValidationError


def test_shape_error_preserves_machine_readable_context() -> None:
    error = ShapeError(field="points", expected="(n, 3)", actual=(4, 2))

    assert isinstance(error, ValidationError)
    assert error.field == "points"
    assert error.expected == "(n, 3)"
    assert error.actual == (4, 2)
    assert str(error) == "points has shape (4, 2); expected (n, 3)"


def test_unit_error_names_actual_and_required_si_unit() -> None:
    error = UnitError(field="length", expected="m", actual="mm")

    assert isinstance(error, ValidationError)
    assert error.field == "length"
    assert error.expected == "m"
    assert error.actual == "mm"
    assert str(error) == "length uses unit 'mm'; expected SI unit 'm'"
