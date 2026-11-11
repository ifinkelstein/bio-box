"""Intended path, step A: build the release's search pack (biobox.search.pack.build) in the shadow root.

The project code is called unchanged; its private stage functions are wrapped only to log
wall time and peak RSS per stage, and orfs.dark_orfs to log progress.

usage: s7_pack.py [--code15-as-11]
  --code15-as-11  workaround for the first failure (pyhmmer has no genetic code 15, which
                  pyrodigal-gv assigns to TAG-recoded phages): translate those contigs with code 11.
"""
import json
import resource
import sys
import time
import traceback
from pathlib import Path

from biobox.search import orfs, pack

SH = Path("/tmp/if2453/ssrnap/shadow")
REL = "rS1.pilot1"


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


for name in ("_contigs", "_slots", "_build_unique", "_build_igr", "_build_dark_orf", "_build_models"):
    f = getattr(pack, name)

    def w(*a, _f=f, _n=name, **k):
        t = time.time()
        print(f"{time.strftime('%H:%M:%S')} start {_n}", flush=True)
        r = _f(*a, **k)
        print(f"{time.strftime('%H:%M:%S')} done {_n} {time.time() - t:.1f} s, peak RSS {rss_gb():.2f} GB", flush=True)
        return r

    setattr(pack, name, w)

_dark = orfs.dark_orfs
N = {"contigs": 0, "orfs": 0, "code15": 0, "t": time.time()}


def dark(seq, circ, genes, code=11, *a, **k):
    if code == 15 and "--code15-as-11" in sys.argv:
        N["code15"] += 1
        code = 11
    r = _dark(seq, circ, genes, code, *a, **k)
    N["contigs"] += 1
    N["orfs"] += len(r)
    if N["contigs"] % 2000 == 0:
        print(f"{time.strftime('%H:%M:%S')}   dark_orf: {N['contigs']} contigs, {N['orfs']} ORFs, "
              f"{N['code15']} code-15 contigs, {time.time() - N['t']:.0f} s, peak RSS {rss_gb():.2f} GB", flush=True)
    return r


orfs.dark_orfs = dark
t0 = time.time()
try:
    man = pack.build(SH, REL, models={"ssrnap-1": Path("/tmp/if2453/ssrnap/hmm/query_lib.hmm")})
    print(json.dumps(man["targets"], indent=1), flush=True)
    print("pack_id", man["pack_id"], f"build {time.time() - t0:.1f} s; dark_orf {N}", flush=True)
    t1 = time.time()
    print("verify:", pack.verify(SH, REL), f"{time.time() - t1:.1f} s", flush=True)
except BaseException:
    traceback.print_exc()
    print(f"dark_orf progress at failure: {N}", flush=True)
    sys.exit(1)
finally:
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print(f"total {time.time() - t0:.1f} s wall, {(ru.ru_utime + ru.ru_stime) / 60:.1f} CPU-min, "
          f"peak RSS {rss_gb():.2f} GB", flush=True)
