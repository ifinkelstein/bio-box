"""Shared helpers for the ssRNAP walk-through (ad hoc; not project code)."""
from collections import defaultdict

import pyarrow as pa

from biobox.catalog import open_catalog

CAT = "/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1/catalog/release=rS1.pilot1.duckdb"
OUT = "/tmp/if2453/ssrnap/out"
Z = 1e8

KNOWN = {  # label -> accession in INPHARED (found by header text, s0_probe.py)
    "T7": "AY264774", "T3": "NC_003298", "SP6": "NC_004831", "K11": "NC_011043", "phiKMV": "NC_005045",
    "N4": "NC_008720", "gh-1": "NC_004665", "LKA1": "NC_009936", "phiYeO3-12": "NC_001271", "K1F": "NC_007456",
    "P60": "NC_003390", "Syn5": "NC_009531", "P-SSP7": "NC_006882", "phiA1122": "NC_004777", "Berlin": "NC_008694",
}


def catalog():
    return open_catalog(CAT, local_dir="/tmp/if2453/ssrnap/catalog_cache", threads=8, memory_limit="12GB")


def read_domtbl(path):
    """{(target, query): {'seq_score', 'seq_evalue', 'tlen', 'qlen', 'domains': [...]}}"""
    out = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        key = (f[0], f[3])
        h = out.setdefault(key, {"tlen": int(f[2]), "qlen": int(f[5]), "seq_evalue": float(f[6]),
                                 "seq_score": float(f[7]), "domains": []})
        h["domains"].append({"score": float(f[13]), "ievalue": float(f[12]), "hmm_from": int(f[15]),
                             "hmm_to": int(f[16]), "env_from": int(f[19]), "env_to": int(f[20])})
    return out


def by_query(dom):
    q = defaultdict(dict)
    for (t, qn), h in dom.items():
        q[qn][t] = h
    return q


def coverage(h):
    """Fraction of model positions covered by the union of domain hmm ranges."""
    cov = set()
    for d in h["domains"]:
        cov.update(range(d["hmm_from"], d["hmm_to"] + 1))
    return len(cov) / h["qlen"]


def uid_table(con, name, uids_hex, extra=None):
    """Register a DuckDB view `name` of protein_uid blobs (plus optional columns)."""
    cols = {"protein_uid": pa.array([bytes.fromhex(u) for u in uids_hex], pa.binary(16)).cast(pa.binary())}
    for k, v in (extra or {}).items():
        cols[k] = v
    con.register(name, pa.table(cols))


def hist(vals, edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append((lo, hi, sum(1 for v in vals if lo <= v < hi)))
    return out
