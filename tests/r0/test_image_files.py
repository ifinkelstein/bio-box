"""Ring R0: every repository file the tests read is copied into the test image.

The image is built only on `main`, so a test that reads a top-level document the image does not ship fails there
and nowhere before it (2026-11-10: `refdb.yaml`, `RUNBOOK.md`; 2026-11-13: `GUIDE.md`, `AGENT_GUIDE.md`). This test
reads `container/Dockerfile` and the test sources: each top-level file a test names as `ROOT / "<file>"` or
`REPO / "<file>"` must be the source of a `COPY` line, and so must each top-level directory named the same way.

The Dockerfile is not itself in the image, so the test skips there; it runs in `ci` and on the pod's mirror.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

pytestmark = [pytest.mark.r0, pytest.mark.covers("C-06")]

ROOT = Path(__file__).resolve().parents[2]
DOCKERFILE = ROOT / "container" / "Dockerfile"
NAMED = re.compile(r"""\b(?:ROOT|REPO)\s*/\s*["']([A-Za-z0-9_.-]+)["']""")


def copied() -> set[str]:
    """The first path component of every source of a `COPY` line (flags and the destination left out)."""
    out: set[str] = set()
    for line in DOCKERFILE.read_text().splitlines():
        words = shlex.split(line, comments=True)
        if not words or words[0] != "COPY":
            continue
        sources = [w for w in words[1:] if not w.startswith("--")][:-1]
        out.update(Path(s).parts[0] for s in sources)
    return out


def named() -> dict[str, list[str]]:
    """{top-level name: the test files that name it} for names that exist in the repository."""
    out: dict[str, list[str]] = {}
    for path in sorted((ROOT / "tests").rglob("*.py")):
        for name in sorted(set(NAMED.findall(path.read_text()))):
            if (ROOT / name).exists():
                out.setdefault(name, []).append(str(path.relative_to(ROOT)))
    return out


def test_every_top_level_file_a_test_reads_is_copied_into_the_image():
    if not DOCKERFILE.is_file():
        pytest.skip("container/Dockerfile is not in the image: this check runs from a checkout")
    have = copied()
    assert {"tests", "src", "RUNBOOK.md"} <= have, f"the COPY lines were not read: {sorted(have)}"
    missing = {name: files for name, files in named().items() if name not in have}
    assert not missing, f"read by tests and not copied into the image by container/Dockerfile: {missing}"
