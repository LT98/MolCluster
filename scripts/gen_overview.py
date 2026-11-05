"""Regenerate the API section of each subpackage's OVERVIEW.md from the code.

Each `src/mofsbu/<pkg>/OVERVIEW.md` has a hand-written top (purpose, rules, boundaries) and
a generated section between the two markers below: every public function and class, its
signature, the first line of its docstring, and which modules outside the package import
it. Only the generated section is rewritten.

    python scripts/gen_overview.py            # rewrite every OVERVIEW.md that has markers
    python scripts/gen_overview.py energy     # one package
    python scripts/gen_overview.py --check    # exit 1 if any generated section is stale
"""
from __future__ import annotations

import ast
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "mofsbu"
BEGIN = "<!-- BEGIN GENERATED: scripts/gen_overview.py — do not edit by hand -->"
END = "<!-- END GENERATED -->"
CALLER_DIRS = (SRC, ROOT / "scripts")


def _first_line(node: ast.AST) -> str:
    doc = ast.get_docstring(node) or ""
    para = doc.strip().split("\n\n")[0]
    return " ".join(para.split())


def _signature(fn: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    a = fn.args
    parts: list[str] = []
    pos = a.posonlyargs + a.args
    for arg in pos:
        if arg.arg not in ("self", "cls"):
            parts.append(arg.arg)
    if a.vararg:
        parts.append("*" + a.vararg.arg)
    elif a.kwonlyargs:
        parts.append("*")
    parts += [k.arg for k in a.kwonlyargs]
    if a.kwarg:
        parts.append("**" + a.kwarg.arg)
    ret = f" -> {ast.unparse(fn.returns)}" if fn.returns else ""
    return f"{fn.name}({', '.join(parts)}){ret}"


def _public(node: ast.AST) -> bool:
    return isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
        and not node.name.startswith("_")


def _callers(pkg: str) -> dict[tuple[str, str], set[str]]:
    """(module, name) -> set of importing locations outside the package."""
    out: dict[tuple[str, str], set[str]] = defaultdict(set)
    prefix = f"mofsbu.{pkg}"
    for base in CALLER_DIRS:
        for path in base.rglob("*.py"):
            rel = path.relative_to(ROOT)
            if rel.parts[:3] == ("src", "mofsbu", pkg):
                continue
            where = (".".join(rel.with_suffix("").parts[2:]) if rel.parts[0] == "src"
                     else str(rel))
            try:
                tree = ast.parse(path.read_text())
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and node.module and \
                        (node.module == prefix or node.module.startswith(prefix + ".")):
                    mod = node.module[len(prefix) + 1:] or "__init__"
                    for alias in node.names:
                        out[(mod, alias.name)].add(where)
    return out


def _reexports(pkg_dir: Path, pkg: str) -> dict[tuple[str, str], str]:
    """(module, name) -> defining module, for names a module imports from a sibling.

    `runner` importing `backend_for` from `energy.relax` is a use of the function defined in
    `energy.backends`; relax only passes it along. "__init__" stands for the package root.
    """
    prefix = f"mofsbu.{pkg}."
    out: dict[tuple[str, str], str] = {}
    for path in pkg_dir.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        here = ".".join(path.relative_to(pkg_dir).with_suffix("").parts)
        for node in ast.parse(path.read_text()).body:      # module level only
            if isinstance(node, ast.ImportFrom) and node.module \
                    and node.module.startswith(prefix):
                src = node.module[len(prefix):]
                for alias in node.names:
                    out[(here, alias.asname or alias.name)] = src
    return out


def _resolve(key: tuple[str, str], reexp: dict[tuple[str, str], str]) -> tuple[str, str]:
    seen = set()
    while key in reexp and key not in seen:
        seen.add(key)
        key = (reexp[key], key[1])
    return key


def render(pkg: str) -> str:
    pkg_dir = SRC / pkg
    callers = _callers(pkg)
    reexp = _reexports(pkg_dir, pkg)
    # a name imported through the root or a sibling is attributed to its defining module
    resolved: dict[tuple[str, str], set[str]] = defaultdict(set)
    for key, who in callers.items():
        resolved[_resolve(key, reexp)] |= who
    callers = resolved
    at_root = {_resolve(("__init__", name), reexp)
               for (m, name) in reexp if m == "__init__"}

    lines = [BEGIN, "", "## API (generated)", "",
             "`*` = re-exported from the package root. *Used by* lists importers in "
             "`src/` and `scripts/` outside this package; blank = internal or tests only.", ""]
    for path in sorted(pkg_dir.rglob("*.py")):
        if path.name == "__init__.py" or "__pycache__" in path.parts:
            continue
        mod = ".".join(path.relative_to(pkg_dir).with_suffix("").parts)
        tree = ast.parse(path.read_text())
        items = [n for n in tree.body if _public(n)]
        n_lines = len(path.read_text().splitlines())
        lines += [f"### `{mod}` ({n_lines} lines)", ""]
        summary = _first_line(tree)
        if summary:
            lines += [summary, ""]
        if not items:
            lines += ["*(no public names)*", ""]
            continue
        lines += ["| Name | Does | Used by |", "|---|---|---|"]
        for n in items:
            sig = f"class {n.name}" if isinstance(n, ast.ClassDef) else _signature(n)
            star = "*" if (mod, n.name) in at_root else ""
            used = ", ".join(f"`{c}`" for c in sorted(callers.get((mod, n.name), ())))
            does = _first_line(n).replace("|", "\\|")
            lines.append(f"| `{sig}`{star} | {does} | {used} |")
        lines.append("")
    lines.append(END)
    return "\n".join(lines)


def update(pkg: str, *, check: bool) -> bool:
    """Return True if the file is (or now is) current."""
    path = SRC / pkg / "OVERVIEW.md"
    text = path.read_text()
    if BEGIN not in text or END not in text:
        print(f"{path.relative_to(ROOT)}: no generated markers, skipped")
        return True
    head, rest = text.split(BEGIN, 1)
    tail = rest.split(END, 1)[1]
    new = head + render(pkg) + tail
    if new == text:
        return True
    if check:
        print(f"{path.relative_to(ROOT)} is stale: run scripts/gen_overview.py {pkg}")
        return False
    path.write_text(new)
    print(f"updated {path.relative_to(ROOT)}")
    return True


def main(argv: list[str]) -> int:
    check = "--check" in argv
    pkgs = [a for a in argv if not a.startswith("-")] or sorted(
        p.parent.name for p in SRC.glob("*/OVERVIEW.md"))
    ok = all([update(p, check=check) for p in pkgs])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
