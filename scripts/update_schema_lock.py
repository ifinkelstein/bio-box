"""Rewrite tests/registry/schema_lock.json from the registry after an allowed change (LY-09).

    python scripts/update_schema_lock.py

Refuses if the change breaks LY-09 (a removed, renamed, retyped or reordered column, or an
added column that is not nullable).
"""

import json
import sys
from pathlib import Path

from biobox.registry.schemas import evolution_problems, schema_lock

LOCK = Path(__file__).resolve().parents[1] / "tests" / "registry" / "schema_lock.json"
old = json.loads(LOCK.read_text()) if LOCK.exists() else {}
problems = evolution_problems(old)
if problems:
    sys.exit("refused:\n" + "\n".join(problems))
LOCK.write_text(json.dumps(schema_lock(), indent=1) + "\n")
print(f"wrote {LOCK}")
