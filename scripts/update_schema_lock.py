"""Rewrite tests/registry/schema_lock.json from the registry after an allowed change (LY-09).

    python scripts/update_schema_lock.py
    python scripts/update_schema_lock.py --breaking DEC-17

Refuses if the change breaks LY-09 (a removed, renamed, retyped or reordered column, or an
added column that is not nullable), unless `--breaking` names the decision that allows it.
A breaking change is allowed only before the first release that holds the table: once a
release is out, its files must stay readable by the catalog (LY-09).
"""

import argparse
import json
import sys
from pathlib import Path

from biobox.registry.schemas import evolution_problems, schema_lock

LOCK = Path(__file__).resolve().parents[1] / "tests" / "registry" / "schema_lock.json"
ap = argparse.ArgumentParser()
ap.add_argument("--breaking", metavar="DECISION", help="decision record allowing the breaking change")
args = ap.parse_args()
old = json.loads(LOCK.read_text()) if LOCK.exists() else {}
problems = evolution_problems(old)
if problems and not args.breaking:
    sys.exit("refused:\n" + "\n".join(problems))
LOCK.write_text(json.dumps(schema_lock(), indent=1) + "\n")
if problems:
    print(f"breaking change allowed by {args.breaking}:\n" + "\n".join(problems))
print(f"wrote {LOCK}")
