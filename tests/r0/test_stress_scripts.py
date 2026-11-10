"""Ring R0: the stress scripts call the readers and builders with their current required arguments.

Static only (card A-06): every `.py` under `scripts/stress` and the `src/biobox` sources are parsed with
`ast`; nothing is imported. A signature change that makes a script call stale fails here, naming the
file and line, without anyone running the scripts.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

pytestmark = [pytest.mark.r0, pytest.mark.covers("A-06")]

ROOT = Path(__file__).resolve().parents[2]
STRESS = ROOT / "scripts" / "stress"
SRC = ROOT / "src" / "biobox"

# (module qualifier used at call sites, function) -> source file defining it
CHECKED = {
    ("tokens", "run"): SRC / "l3" / "tokens.py",
    ("postings", "run"): SRC / "l3" / "postings.py",
    ("dim", "run"): SRC / "dim.py",
    ("cluster", "run"): SRC / "l2" / "cluster.py",
    (None, "build_manifest"): SRC / "catalog.py",
    (None, "read_tapes"): SRC / "l3" / "tokens.py",
    (None, "iter_tapes"): SRC / "l3" / "tokens.py",
    (None, "libraries_in"): SRC / "l3" / "tokens.py",
}
# bare names (imported directly) and attribute calls on these qualifiers are both checked
BARE_OK = {"build_manifest", "read_tapes", "iter_tapes", "libraries_in"}
QUALIFIERS = {"build_manifest": {"catalog"}, "read_tapes": {"tokens"}, "iter_tapes": {"tokens"},
              "libraries_in": {"tokens"}}


def _required(path: Path, name: str) -> list[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            a = node.args
            pos = a.posonlyargs + a.args
            n_def = len(a.defaults)
            req = [p.arg for p in pos[:len(pos) - n_def]]
            req += [k.arg for k, d in zip(a.kwonlyargs, a.kw_defaults) if d is None]
            return req
    raise AssertionError(f"{path}: no top-level def {name}")


def _target(call: ast.Call):
    """(qualifier, name) of a call to `name(...)` or `qualifier.name(...)`; qualifier None if bare."""
    f = call.func
    if isinstance(f, ast.Name):
        return None, f.id
    if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
        return f.value.id, f.attr
    return "?", getattr(f, "attr", None)


def _unpacked(call: ast.Call) -> bool:
    return any(isinstance(a, ast.Starred) for a in call.args) or any(k.arg is None for k in call.keywords)


def problems() -> list[str]:
    required = {key: _required(path, key[1]) for key, path in CHECKED.items()}
    out = []
    for py in sorted(STRESS.rglob("*.py")):
        rel = py.relative_to(ROOT)
        tree = ast.parse(py.read_text(), filename=str(py))
        for call in (n for n in ast.walk(tree) if isinstance(n, ast.Call)):
            qual, name = _target(call)
            if _unpacked(call):
                continue
            given = len(call.args) + len(call.keywords)
            if name == "SeqStore":
                has_stages = len(call.args) >= 2 or any(k.arg and "stages" in k.arg for k in call.keywords)
                if not has_stages:
                    out.append(f"{rel}:{call.lineno}: SeqStore called without its stages")
                continue
            for (q, fn), req in required.items():
                if name != fn:
                    continue
                ok_qual = qual == q if q is not None else (qual is None or qual in QUALIFIERS[fn])
                if not ok_qual:
                    continue
                if given < len(req):
                    out.append(f"{rel}:{call.lineno}: {fn} given {given} arguments, requires {len(req)} ({', '.join(req)})")
    return out


def test_a06_stress_scripts_call_current_signatures():
    bad = problems()
    assert not bad, "stale calls in scripts/stress:\n" + "\n".join(bad)


def test_a06_checker_sees_a_stale_call():
    """The checker itself: a fixture with stale calls is flagged (guards against a silently empty walk)."""
    src = "store = SeqStore(root)\ntokens.run(root, 'r')\nok = SeqStore(root, ['S0'])\n"
    calls = [n for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Call)]
    stale = [c for c in calls if _target(c)[1] == "SeqStore" and len(c.args) < 2]
    assert [c.lineno for c in stale] == [1]
    assert len(_required(CHECKED[("tokens", "run")], "run")) == 3
    assert _required(CHECKED[(None, "read_tapes")], "read_tapes") == ["root", "stages"]
