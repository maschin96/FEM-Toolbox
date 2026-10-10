"""Extract declared API signatures with stdlib AST; never import optional backends."""

from __future__ import annotations

import ast
import json
from pathlib import Path
from typing import Any


def _signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> dict[str, Any]:
    return {
        "parameters": ast.unparse(node.args),
        "returns": None if node.returns is None else ast.unparse(node.returns),
        "decorators": [ast.unparse(d) for d in node.decorator_list],
    }


def module_contract(path: Path) -> dict[str, Any]:
    """Inventory explicit __all__, public functions, fields and methods.

    Implementation bodies/locations do not change signatures. Missing declared
    definitions are rejected rather than disappearing silently. Reexports are
    intentionally outside this first uncertainty-module contract.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    exports: list[str] | None = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets
        ):
            exports = list(ast.literal_eval(node.value))
    if (
        exports is None
        or len(set(exports)) != len(exports)
        or not all(isinstance(name, str) for name in exports)
    ):
        raise ValueError("module needs unique literal __all__ names")
    records = {}
    for node in tree.body:
        if (
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name in exports
        ):
            records[node.name] = {"kind": "function", **_signature(node)}
        elif isinstance(node, ast.ClassDef) and node.name in exports:
            records[node.name] = {
                "kind": "class",
                "bases": [ast.unparse(b) for b in node.bases],
                "decorators": [ast.unparse(d) for d in node.decorator_list],
                "fields": {
                    n.target.id: {
                        "annotation": ast.unparse(n.annotation),
                        "default": None if n.value is None else ast.unparse(n.value),
                    }
                    for n in node.body
                    if isinstance(n, ast.AnnAssign)
                    and isinstance(n.target, ast.Name)
                    and not n.target.id.startswith("_")
                },
                "methods": {
                    n.name: _signature(n)
                    for n in node.body
                    if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and (not n.name.startswith("_") or n.name == "__init__")
                },
            }
    if set(exports) != set(records):
        raise ValueError(f"missing definitions: {sorted(set(exports) - set(records))}")
    return {"schema": 1, "exports": sorted(exports), "symbols": records}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(module_contract(args.source), ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
