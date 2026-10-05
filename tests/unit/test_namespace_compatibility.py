"""Public import behavior across the Femantis namespace transition."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def run_python(code: str) -> None:
    """Exercise clean module state with this checkout ahead of installed sources."""
    environment = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=environment,
        cwd=ROOT,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_canonical_root_import_is_lightweight_and_silent() -> None:
    run_python(
        """
import importlib.abc
import sys
import warnings

class DenyOptional(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'petsc4py', 'matplotlib'}:
            raise AssertionError('optional dependency loaded: ' + fullname)

sys.meta_path.insert(0, DenyOptional())
with warnings.catch_warnings(record=True) as observed:
    warnings.simplefilter('always')
    import femantis
assert femantis.__name__ == 'femantis'
assert not observed, observed
assert not any(name.startswith(('petsc4py', 'matplotlib', 'femantis.solvers'))
               for name in sys.modules)
"""
    )


@pytest.mark.parametrize("first", ["femantis", "fem_toolbox"])
def test_all_inventoried_modules_share_identity(first: str) -> None:
    inventory = json.loads(
        (ROOT / "Doc/brand/femantis-name-inventory.json").read_text()
    )
    modules = inventory["python_modules"]
    run_python(
        f"""
import importlib
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
importlib.import_module({first!r})
modules = {modules!r}
for entry in modules:
    a, b = entry['canonical'], entry['legacy']
    order = (a, b) if {first!r} == 'femantis' else (b, a)
    left, right = [importlib.import_module(name) for name in order]
    assert left is right, (a, b)
    assert left.__name__ == a, (a, left.__name__)
    assert left.__spec__.name == a, (a, left.__spec__.name)
    assert left.__package__ == (a if hasattr(left, '__path__') else a.rpartition('.')[0])
    if '.' in a:
        parent, _, attribute = b.rpartition('.')
        assert getattr(importlib.import_module(parent), attribute) is left
from femantis.mesh import Mesh
from fem_toolbox.mesh import Mesh as OldMesh
from femantis.errors import ValidationError
from fem_toolbox.errors import ValidationError as OldError
assert Mesh is OldMesh
assert ValidationError is OldError
assert Mesh.__module__ == 'femantis.mesh.model'
assert isinstance(ValidationError('invalid', field='input'), OldError)
"""
    )


def test_legacy_root_warns_once_at_the_caller() -> None:
    run_python(
        """
import warnings
import femantis
with warnings.catch_warnings(record=True) as observed:
    warnings.simplefilter('always', DeprecationWarning)
    import fem_toolbox
    import fem_toolbox
    from fem_toolbox.mesh import Mesh
assert fem_toolbox is femantis
assert len(observed) == 1, observed
warning = observed[0]
assert warning.category is DeprecationWarning
assert 'femantis' in str(warning.message)
assert 'v2.x' in str(warning.message)
assert warning.filename == '<string>', warning.filename
"""
    )


def test_unknown_legacy_module_raises_normal_import_error() -> None:
    run_python(
        """
import importlib
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
for name in ('femantis.does_not_exist', 'fem_toolbox.does_not_exist'):
    try:
        importlib.import_module(name)
    except ModuleNotFoundError as error:
        assert error.name == name, (name, error.name)
    else:
        raise AssertionError(name + ' unexpectedly exists')
"""
    )


@pytest.mark.parametrize("namespace", ["femantis", "fem_toolbox"])
def test_root_imports_do_not_load_numeric_or_optional_modules(namespace: str) -> None:
    run_python(
        f"""
import sys
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
__import__({namespace!r})
for name in ('numpy', 'scipy', 'meshio', 'petsc4py', 'matplotlib'):
    assert name not in sys.modules, name
assert 'femantis.model' not in sys.modules
assert 'femantis.mesh' not in sys.modules
"""
    )


@pytest.mark.parametrize("namespace", ["femantis", "fem_toolbox"])
def test_model_heat_result_and_legacy_vtu_contract(namespace: str) -> None:
    # For a unit-axis Tet4 with k=2 W/(m K), q=8 W/m³ and three nodes at
    # 300 K, the free-node equation k*V*(T1-300)=q*V/4 gives T1=301 K.
    # The total source is q*V=4/3 W; reservoirs remove that same power.
    run_python(
        f"""
import importlib
import tempfile
import warnings
from pathlib import Path
import numpy as np
import meshio
warnings.simplefilter('ignore', DeprecationWarning)
Mesh = importlib.import_module({namespace!r} + '.mesh').Mesh
Model = importlib.import_module({namespace!r} + '.model').Model
write_vtu = importlib.import_module({namespace!r} + '.io').write_vtu
read_vtu = importlib.import_module({namespace!r} + '.io').read_vtu
mesh = Mesh(
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),
    tetra4=np.array([[0,1,2,3]]),
    point_ids=np.array([10,11,12,13]),
    volume_sets={{'body': np.array([0])}},
)
document = {{
    'schema_version': 1, 'units': 'SI', 'physics': 'heat',
    'materials': {{'T': {{'thermal': {{'conductivity': 2., 'unit': 'W/(m K)'}}}}}},
    'regions': {{'T': {{'body': 'thermal'}}}},
    'boundary': {{'dirichlet': [{{'field': 'T', 'nodes': [0,2,3],
                                'value': 300., 'unit': 'K'}}]}},
    'loads': {{'source': {{'value': 8., 'unit': 'W/m^3'}}}},
}}
field = Model(mesh=mesh, document=document).solve().field_result
# 1e-12 K/W covers roundoff in this small well-conditioned SI example.
np.testing.assert_allclose(field.temperature, [300.,301.,300.,300.], rtol=0, atol=1e-12)
np.testing.assert_allclose(field.source_power, 4/3, rtol=0, atol=1e-12)
np.testing.assert_allclose(field.dirichlet_inflow, -4/3, rtol=0, atol=1e-12)
assert abs(field.power_balance) < 1e-12
with tempfile.TemporaryDirectory() as temporary:
    path = Path(temporary) / 'mesh.vtu'
    write_vtu(path, mesh)
    data = meshio.read(path)
    assert 'fem_toolbox:external_id' in data.point_data
    assert 'fem_toolbox:set:body' in data.cell_data
    np.testing.assert_array_equal(read_vtu(path).point_ids, mesh.point_ids)
"""
    )


def test_alias_import_and_reload_preserve_module_metadata() -> None:
    run_python(
        """
import importlib
import warnings
warnings.simplefilter('ignore', DeprecationWarning)
legacy = importlib.import_module('fem_toolbox.conventions')
canonical = importlib.import_module('femantis.conventions')
assert legacy is canonical
assert importlib.reload(legacy) is canonical
assert legacy.__spec__.name == 'femantis.conventions'
assert legacy.__loader__ is legacy.__spec__.loader
"""
    )
