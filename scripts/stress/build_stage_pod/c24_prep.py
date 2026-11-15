"""Card C-24 rehearsal: the fixture inputs of RUNBOOK Appendix B (its `prep.py`, unchanged in substance): the genome
and object lists of stages S0 (lambda, M. genitalium) and S1 (T4) and the test HMM library, made with the helpers of
tests/r2/test_builders_out.py, in WORK/inputs. In the image, with the tree's src on PYTHONPATH:

    python c24_prep.py TREE WORK
"""
import shutil
import sys
from pathlib import Path

tree, work = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(tree / "tests" / "r2"))
import test_builders_out as c08  # noqa: E402  the helpers the tests use

w = c08.Work(work / "inputs")
w.prepare()
scratch_root = work / "libroot"
for step in c08.stage_steps("S0", w)[:3]:  # inputs, ingest, genecall: the proteins the HMMs are built from
    step.direct(scratch_root)
c08.write_library(scratch_root, w.lib)
shutil.rmtree(scratch_root)
print("inputs written under", w.base)
