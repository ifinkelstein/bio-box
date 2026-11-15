"""Card W-32: the tables' digests that tests/r2/test_w27_flanks_figure.py compares with (BASE_TABLES).

Run with the code before the card (the `src` of 4469fab, and only that, on PYTHONPATH) and the card's test file, in
the test image on the pod:

    PYTHONPATH=<base src> python base_digests.py <tree>/tests/r2/test_w27_flanks_figure.py <scratch dir>

It builds the test's constructed releases under the scratch directory, runs every call of the test's CALLS with that
code and prints BASE_TABLES as JSON (name -> [files, digest, README.txt sha256]). Not `python -I`: PYTHONPATH must be
read.
"""

import importlib.util
import json
import sys
from pathlib import Path


def main() -> None:
    test_file, scratch = Path(sys.argv[1]), Path(sys.argv[2])
    spec = importlib.util.spec_from_file_location("w27_for_base", test_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import biobox

    print(f"biobox from {Path(biobox.__file__).parent}", file=sys.stderr)
    scratch.mkdir(parents=True, exist_ok=False)
    print(json.dumps(mod.base_digests(scratch), indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
