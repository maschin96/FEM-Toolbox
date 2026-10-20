"""Extract declared API signatures with stdlib AST; never import optional backends."""

from __future__ import annotations

import ast
import json
import re
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


def _definition(
    node: ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef,
) -> dict[str, Any]:
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        return {
            "kind": "function",
            **_signature(node),
            "documentation": ast.get_docstring(node),
            "direct_raise_sites": _raises(node),
        }
    return {
        "kind": "class",
        "bases": [ast.unparse(b) for b in node.bases],
        "decorators": [ast.unparse(d) for d in node.decorator_list],
        "documentation": ast.get_docstring(node),
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
            n.name: {
                **_signature(n),
                "documentation": ast.get_docstring(n),
                "direct_raise_sites": _raises(n),
            }
            for n in node.body
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
            and (
                not n.name.startswith("_")
                or n.name
                in {
                    "__init__",
                    "__call__",
                    "__getitem__",
                    "__iter__",
                    "__len__",
                    "__enter__",
                    "__exit__",
                }
            )
        },
    }


def _raises(node: ast.AST) -> list[str]:
    return sorted(
        {
            ast.unparse(n.exc.func)
            for n in ast.walk(node)
            if isinstance(n, ast.Raise) and isinstance(n.exc, ast.Call)
        }
    )


def collect_api(root: Path) -> dict[str, Any]:
    """Source inventory without importing FEM, NumPy, PETSc or Matplotlib.

    Local public definitions, annotated datatypes/constants and explicit
    reexports are included. Incidental third-party imports are excluded.
    Public helpers remain inventory entries, not implicit stability promises.
    """
    modules: dict[str, Any] = {}
    for path in sorted((root / "src/fem_toolbox").rglob("*.py")):
        parts = list(path.relative_to(root / "src").with_suffix("").parts)
        if parts[-1] == "__init__":
            parts.pop()
        module = ".".join(parts)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        definitions, reexports, aliases, imported = {}, {}, {}, {}
        exports = None
        for node in tree.body:
            if isinstance(node, ast.ImportFrom) and node.module:
                origin = node.module
                if node.level:
                    parent = (
                        module
                        if path.name == "__init__.py"
                        else module.rsplit(".", 1)[0]
                    )
                    origin = (
                        ".".join(
                            parent.split(".")[: len(parent.split(".")) - node.level + 1]
                        )
                        + "."
                        + node.module
                    )
                for item in node.names:
                    imported[item.asname or item.name] = origin + "." + item.name
                    if item.asname == item.name and origin.startswith("fem_toolbox"):
                        reexports[item.name] = origin + "." + item.name
            elif isinstance(
                node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            ) and not node.name.startswith("_"):
                definitions[node.name] = _definition(node)
            elif (
                isinstance(node, ast.AnnAssign)
                and isinstance(node.target, ast.Name)
                and not node.target.id.startswith("_")
            ):
                aliases[node.target.id] = {
                    "annotation": ast.unparse(node.annotation),
                    "value": None if node.value is None else ast.unparse(node.value),
                }
            elif isinstance(node, ast.Assign):
                for target in node.targets:
                    if not isinstance(target, ast.Name):
                        continue
                    if target.id == "__all__":
                        exports = list(ast.literal_eval(node.value))
                    elif not target.id.startswith("_") or target.id == "__version__":
                        if (
                            isinstance(node.value, ast.Attribute)
                            and isinstance(node.value.value, ast.Name)
                            and node.value.value.id in imported
                        ):
                            reexports[target.id] = (
                                imported[node.value.value.id] + "." + node.value.attr
                            )
                        else:
                            aliases[target.id] = {"value": ast.unparse(node.value)}
        if exports is not None:
            if len(exports) != len(set(exports)) or not all(
                isinstance(n, str) for n in exports
            ):
                raise ValueError(f"invalid exports in {module}")
            for name in exports:
                if name in imported and name not in definitions:
                    reexports[name] = imported[name]
                if (
                    name not in definitions
                    and name not in aliases
                    and name not in reexports
                ):
                    raise ValueError(f"unresolved reexport {module}.{name}")
        modules[module] = {
            "source": str(path.relative_to(root)),
            "documentation": ast.get_docstring(tree),
            "declared_exports": exports,
            "definitions": definitions,
            "reexports": reexports,
            "aliases": aliases,
        }

    # Verify every declared reexport, including alias chains, without imports.
    def resolve(reference: str, visited: set[str]) -> None:
        if reference in visited:
            raise ValueError(f"cyclic reexport {reference}")
        if reference in modules:
            return
        if "." not in reference:
            raise ValueError(f"missing reexport {reference}")
        module, name = reference.rsplit(".", 1)
        if module not in modules:
            raise ValueError(f"missing reexport module {reference}")
        entry = modules[module]
        if name in entry["definitions"] or name in entry["aliases"]:
            return
        if name in entry["reexports"]:
            resolve(entry["reexports"][name], visited | {reference})
            return
        raise ValueError(f"missing reexport definition {reference}")

    for entry in modules.values():
        for reference in entry["reexports"].values():
            resolve(reference, set())
    return {"schema": 2, "modules": modules}


def validate_api_map(
    root: Path, inventory: dict[str, Any], mapping: dict[str, Any]
) -> None:
    canonical = root / "Doc/29_PUBLIC_API_CONTRACT.md"
    if canonical.exists():
        for line in canonical.read_text(encoding="utf-8").splitlines():
            match = re.match(r"\| `(fem_toolbox[^`]+)` \| ([^|]+)\|", line)
            if match:
                module, names = match.groups()
                required = set(re.findall(r"`([^`]+)`", names))
                actual = set(
                    mapping["modules"].get(module, {}).get("stable_symbols", [])
                )
                if required - actual:
                    raise ValueError(
                        f"missing canonical stable API: {module}: {sorted(required - actual)}"
                    )
    if set(mapping["modules"]) != set(inventory["modules"]):
        raise ValueError("API module coverage is incomplete")
    for name, record in inventory["modules"].items():
        entry = mapping["modules"][name]
        if entry["classification"] not in {"stable", "research", "internal", "mixed"}:
            raise ValueError(f"invalid classification: {name}")
        if not entry["rationale"] or not entry["contract"] or not entry["examples"]:
            raise ValueError(f"missing API contract/example: {name}")
        if not (root / entry["contract"]).is_file():
            raise ValueError(f"missing API reference: {name}")
        public = (
            set(record["definitions"])
            | set(record["reexports"])
            | set(record["aliases"])
        )
        if set(entry["stable_symbols"]) - public:
            raise ValueError(f"unknown stable symbols: {name}")
        for example in entry["examples"]:
            if not (root / example).is_file():
                raise ValueError(f"missing API example: {name}: {example}")


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(c, c) for c in text)


def render_inventory(inventory: dict[str, Any], mapping: dict[str, Any]) -> str:
    """Render readable exact signature listings and source-backed contracts."""
    text = [
        r"\chapter{Versionsgenaues API-Inventar}",
        r"\label{api:inventory}",
        r"\begingroup\raggedright",
        r"\lstset{basicstyle=\scriptsize\ttfamily,breaklines=true,breakatwhitespace=false}",
        "Automatisch aus AST-Signaturen und Quell-Dokumentation erzeugt; keine optionalen Laufzeitimporte. Englische Quellverträge bleiben zur exakten Synchronisation erhalten. Manuelle SI-/Form-/Fehler-/Mutabilitätsregeln stehen in den vorherigen Kapiteln. Forschungs- und interne Einträge erhalten dadurch keine Stabilitätszusage.",
    ]
    for module, record in inventory["modules"].items():
        entry = mapping["modules"][module]
        text += [
            r"\section{" + tex_escape(module) + "}",
            r"\index{"
            + module
            + "@"
            + tex_escape(module.removeprefix("fem_toolbox."))
            .replace(".", r".\protect\allowbreak{}")
            .replace(r"\_", r"\_\protect\allowbreak{}")
            + "}",
            "Einordnung: "
            + tex_escape(entry["classification"])
            + ". "
            + tex_escape(entry["rationale"]),
            r"\par Fachvertrag: \nolinkurl{"
            + entry["contract"]
            + "}. Beispiele/Tests:",
        ]
        text += [r"\nolinkurl{" + example + "}." for example in entry["examples"]]
        if record["reexports"]:
            text.append(r"\subsection*{Reexportierte Einstiegspunkte}")
            for alias, reference in record["reexports"].items():
                text.append(
                    r"\noindent\nolinkurl{"
                    + alias
                    + "} $\\to$ \\nolinkurl{"
                    + reference
                    + "}\\par"
                )
        if record["aliases"]:
            text.append(r"\subsection*{Datentypen und Konstanten}")
            for alias, definition in record["aliases"].items():
                value = (
                    alias
                    + (
                        ": " + definition["annotation"]
                        if "annotation" in definition
                        else ""
                    )
                    + (
                        " = " + definition["value"]
                        if definition.get("value") is not None
                        else ""
                    )
                )
                text += [r"\begin{lstlisting}", value, r"\end{lstlisting}"]
        for name, definition in record["definitions"].items():
            status = (
                "stabiler Einstieg"
                if name in entry["stable_symbols"]
                else "Forschung/intern; keine stabile Anwenderzusage"
            )
            text += [r"\subsection{" + tex_escape(name) + "}", tex_escape(status) + "."]
            if definition["kind"] == "function":
                signature = (
                    name
                    + "("
                    + definition["parameters"]
                    + ")"
                    + (" -> " + definition["returns"] if definition["returns"] else "")
                )
                text += [r"\begin{lstlisting}", signature, r"\end{lstlisting}"]
            else:
                text.append(
                    "Basistypen: "
                    + tex_escape(", ".join(definition["bases"]) or "object")
                    + ". Dekoratoren: "
                    + tex_escape(", ".join(definition["decorators"]) or "keine")
                    + "."
                )
            if definition["documentation"]:
                text += [
                    r"\begin{quote}",
                    tex_escape(definition["documentation"])
                    .replace("\n\n", r"\par ")
                    .replace("\n", " "),
                    r"\end{quote}",
                ]
            else:
                text.append(
                    "Zusätzlich gilt der oben verknüpfte manuelle Fachvertrag; kein eigener Quell-Docstring."
                )
            for field, desc in definition.get("fields", {}).items():
                text += [
                    r"\begin{lstlisting}",
                    field
                    + ": "
                    + desc["annotation"]
                    + (" = " + desc["default"] if desc["default"] is not None else ""),
                    r"\end{lstlisting}",
                ]
            for method, desc in definition.get("methods", {}).items():
                if method == "__init__":
                    label = "Konstruktor"
                else:
                    label = method
                text += [
                    r"\subsubsection*{" + tex_escape(label) + "}",
                    r"\begin{lstlisting}",
                    method
                    + "("
                    + desc["parameters"]
                    + ")"
                    + (" -> " + desc["returns"] if desc["returns"] else ""),
                    r"\end{lstlisting}",
                ]
                if desc["documentation"]:
                    text.append(
                        tex_escape(desc["documentation"])
                        .replace("\n\n", r"\par ")
                        .replace("\n", " ")
                    )
                if desc["direct_raise_sites"]:
                    text.append(
                        "Direkte Raise-Stellen: "
                        + tex_escape(", ".join(desc["direct_raise_sites"]))
                        + "; transitive Fehler siehe Fachvertrag."
                    )
            if definition.get("direct_raise_sites"):
                text.append(
                    "Direkte Raise-Stellen: "
                    + tex_escape(", ".join(definition["direct_raise_sites"]))
                    + "; transitive Fehler siehe Fachvertrag."
                )
    text.append(r"\endgroup")
    return "\n".join(text) + "\n"


def validate_reference(root: Path) -> None:
    inventory = collect_api(root)
    source = root / "Doc/latex/api-inventory.json"
    if inventory != json.loads(source.read_text(encoding="utf-8")):
        raise ValueError("API inventory drift: regenerate and review contracts")
    mapping = json.loads((root / "Doc/latex/api-map.json").read_text(encoding="utf-8"))
    validate_api_map(root, inventory, mapping)
    generated = root / "Doc/latex/chapters/api/inventory.tex"
    if generated.read_text(encoding="utf-8") != render_inventory(inventory, mapping):
        raise ValueError("API rendered reference drift")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, nargs="?")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    if args.all:
        if args.check:
            validate_reference(args.root)
            print("API inventory, module coverage and generated chapter match source")
        else:
            inventory = collect_api(args.root)
            mapping = json.loads(
                (args.root / "Doc/latex/api-map.json").read_text(encoding="utf-8")
            )
            validate_api_map(args.root, inventory, mapping)
            (args.root / "Doc/latex/api-inventory.json").write_text(
                json.dumps(inventory, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            (args.root / "Doc/latex/chapters/api/inventory.tex").write_text(
                render_inventory(inventory, mapping), encoding="utf-8"
            )
    else:
        if args.source is None:
            parser.error("source or --all is required")
        text = (
            json.dumps(module_contract(args.source), ensure_ascii=False, indent=2)
            + "\n"
        )
        if args.output:
            args.output.write_text(text, encoding="utf-8")
        else:
            print(text, end="")
