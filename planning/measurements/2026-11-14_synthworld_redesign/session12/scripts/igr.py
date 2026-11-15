"""igr.py OUT SEED [path=JSON ...]: card F-02 session 12, scratch, measurement only.

Builds the micro backbone (with the overrides), sequences and replicons (S1-S5) at SEED into OUT, then runs the two
sequence-level assertions of tests/r2/test_synthworld_intergenic.py on it by calling the test functions of the tree
it runs in (the cross-phylum and burst IGR identity test, and the planted-motif scan as migrated by DEC-44 (1)); prints
each one's own report line (the motif test's counts the instances whose best windows are another planted instance of
the same emission, the case the pre-DEC-44 assertion failed on) and PASS or the assertion."""
import importlib.util
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

from biobox.synthworld import backbone, intergenic, sequences

out = Path(sys.argv[1])
seed = int(sys.argv[2])
params = backbone.load_preset("micro")
sets = []
for a in sys.argv[3:]:
    k, _, v = a.partition("=")
    d = params
    keys = k.split(".")
    for x in keys[:-1]:
        d = d[x]
    d[keys[-1]] = json.loads(v)
    sets.append((k, json.loads(v)))
print(json.dumps({"out": str(out), "seed": seed, "sets": sets}), flush=True)
if out.exists():
    shutil.rmtree(out)
t0 = time.perf_counter()
backbone.build(params, seed, out)
sequences.build(out)
rw = intergenic.build(out)
print(f"built S1-S5 in {time.perf_counter() - t0:.1f} s", flush=True)
spec = importlib.util.spec_from_file_location("tsi", "tests/r2/test_synthworld_intergenic.py")
tsi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tsi)
world = (rw, out, 0.0)
tables = tsi.load_tables(world)
for name in ("test_orthologous_igrs_are_near_identical_in_close_relatives",
             "test_every_planted_motif_is_recovered_at_its_truth_coordinates_by_an_independent_scan"):
    try:
        getattr(tsi, name)(world, tables)
        print(f"{name}: PASS", flush=True)
    except AssertionError:
        print(f"{name}: FAIL\n{traceback.format_exc(limit=1)[-3000:]}", flush=True)
