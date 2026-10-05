"""Lazy aliases preserving canonical module identity for legacy imports."""

import sys
from collections.abc import Sequence
from importlib import import_module
from importlib.abc import Loader, MetaPathFinder
from importlib.machinery import ModuleSpec
from importlib.util import find_spec, spec_from_loader
from types import ModuleType


class _AliasLoader(Loader):
    def __init__(self, canonical: str, canonical_spec: ModuleSpec) -> None:
        self.canonical = canonical
        self.canonical_spec = canonical_spec

    def create_module(self, spec: ModuleSpec) -> ModuleType:
        return import_module(self.canonical)

    def exec_module(self, module: ModuleType) -> None:
        # create_module returns the already initialized canonical module.
        # The import machinery assigns the alias spec after create_module.
        # Restore the canonical spec without executing module code again.
        module.__spec__ = self.canonical_spec


class _LegacyFinder(MetaPathFinder):
    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None = None,
        target: ModuleType | None = None,
    ) -> ModuleSpec | None:
        if not fullname.startswith("fem_toolbox."):
            return None
        canonical = "femantis." + fullname.removeprefix("fem_toolbox.")
        canonical_spec = find_spec(canonical)
        if canonical_spec is None:
            return None
        return spec_from_loader(
            fullname,
            _AliasLoader(canonical, canonical_spec),
            origin=canonical_spec.origin,
            is_package=canonical_spec.submodule_search_locations is not None,
        )


def install_legacy_finder() -> None:
    """Handle only fem_toolbox children, leaving unrelated importers unchanged."""
    if not any(isinstance(finder, _LegacyFinder) for finder in sys.meta_path):
        sys.meta_path.insert(0, _LegacyFinder())
