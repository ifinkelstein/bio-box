"""Protein clustering with MMseqs2 (operation O-CLUSTER).

NR90 clusters all distinct proteins of one gene-calling version; NR50 clusters
the NR90 representatives, and `mergeclusters` maps every protein to its NR50
representative. Each member is then re-aligned to its representative at each
level (`mmseqs align`, thresholds off) so `cluster_member_aln` has one row per
protein per level, representatives included as self rows.

Writes under `<root>/`:

    l2/cluster_member/release=<R>/part-0.parquet
    l2/cluster_member_aln/release=<R>/part-0.parquet

Inputs are written to FASTA sorted by protein_uid and MMseqs2 runs with
`--shuffle 0` and explicit `--threads`, so output does not depend on input order
or on the host's core count.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from blake3 import blake3

from biobox import ledger
from biobox.registry.schemas import table

# UniRef-style: member covered >= 80% by the alignment to its representative (cov-mode 1,
# coverage of the target = member), greedy by length (cluster-mode 2, CD-HIT-like).
LEVELS = {
    "nr90": {"min_seq_id": 0.9, "coverage": 0.8},
    "nr50": {"min_seq_id": 0.5, "coverage": 0.8},
}
COMMON = ["--cov-mode", "1", "--cluster-mode", "2"]


def params_hash(level: str) -> str:
    p = {"level": level, **LEVELS[level], "common": COMMON, "tool": "mmseqs2"}
    return hashlib.sha256(json.dumps(p, sort_keys=True).encode()).hexdigest()[:16]


def cluster_uid(member_uids: list[bytes], phash: str) -> bytes:
    return blake3(b"".join(sorted(member_uids)) + phash.encode()).digest(length=16)


class MMseqs:
    def __init__(self, workdir: Path, threads: int, exe: str = "mmseqs"):
        if shutil.which(exe) is None:
            raise FileNotFoundError(f"{exe} not on PATH")
        self.exe, self.threads, self.workdir = exe, threads, Path(workdir)
        self.log = open(self.workdir / "mmseqs.log", "a")

    def __call__(self, *args: str, threaded: bool = True) -> None:
        cmd = [self.exe, *args, "-v", "1"]
        if threaded:
            cmd += ["--threads", str(self.threads)]
        self.log.write(" ".join(cmd) + "\n")
        self.log.flush()
        subprocess.run(cmd, check=True, stdout=self.log, stderr=subprocess.STDOUT, cwd=self.workdir)


def _read_proteins(root: Path, gene_call_version: str) -> dict[bytes, str]:
    files = sorted((root / "l0seq/protein_seq").glob(f"gene_call_version={gene_call_version}/stage=*/protein_prefix=*/part-*.parquet"))
    if not files:
        raise FileNotFoundError(f"no protein_seq files for gene_call_version={gene_call_version}")
    out: dict[bytes, str] = {}
    for f in files:
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        out.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return out


def _tsv_pairs(path: Path) -> list[tuple[str, str]]:
    return [tuple(line.rstrip("\n").split("\t")[:2]) for line in path.read_text().splitlines() if line]


def ledger_stage(release: str) -> str:
    """rS0.1 -> S0, rS2C.1 -> S2C; any other release name is used as is."""
    return release[1:].split(".")[0] if release.startswith("rS") else release


def run(root: Path, release: str, gene_call_version: str, threads: int = 1, tmp: Path | None = None) -> dict:
    root = Path(root)
    prots = _read_proteins(root, gene_call_version)
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        w = Path(d)
        with open(w / "in.faa", "w") as fh:
            for uid in sorted(prots):
                fh.write(f">{uid.hex()}\n{prots[uid]}\n")
        mm = MMseqs(w, threads)
        mm("createdb", "in.faa", "db", "--dbtype", "1", "--shuffle", "0", "--write-lookup", "1", threaded=False)
        mm("cluster", "db", "clu90", "tmp90", "--min-seq-id", str(LEVELS["nr90"]["min_seq_id"]),
           "-c", str(LEVELS["nr90"]["coverage"]), *COMMON)
        mm("createsubdb", "clu90", "db", "rep90", threaded=False)
        mm("cluster", "rep90", "clu50r", "tmp50", "--min-seq-id", str(LEVELS["nr50"]["min_seq_id"]),
           "-c", str(LEVELS["nr50"]["coverage"]), *COMMON)
        mm("mergeclusters", "db", "clu50", "clu90", "clu50r")
        assign: dict[str, dict[str, str]] = {}
        alns: dict[str, list[list[str]]] = {}
        for level, clu in (("nr90", "clu90"), ("nr50", "clu50")):
            mm("createtsv", "db", "db", clu, f"{clu}.tsv")
            assign[level] = {m: r for r, m in _tsv_pairs(w / f"{clu}.tsv")}
            mm("align", "db", "db", clu, f"aln_{level}", "-a", "-e", "inf", "--min-seq-id", "0", "-c", "0",
               "--alignment-mode", "3")
            mm("convertalis", "db", "db", f"aln_{level}", f"aln_{level}.tsv",
               "--format-output", "query,target,qstart,qend,tstart,tend,fident,alnlen,cigar")
            alns[level] = [ln.split("\t") for ln in (w / f"aln_{level}.tsv").read_text().splitlines() if ln]

    lengths = {u.hex(): len(s) for u, s in prots.items()}
    mem_rows, aln_rows = [], []
    for level in LEVELS:
        phash = params_hash(level)
        members_of: dict[str, list[bytes]] = defaultdict(list)
        for m, r in assign[level].items():
            members_of[r].append(bytes.fromhex(m))
        if set(assign[level]) != set(lengths):
            raise RuntimeError(f"{level}: clustering does not cover every protein")
        cuid = {r: cluster_uid(ms, phash) for r, ms in members_of.items()}
        for m, r in assign[level].items():
            mem_rows.append({"member_uid": bytes.fromhex(m), "level": level, "cluster_uid": cuid[r],
                             "rep_uid": bytes.fromhex(r), "cluster_params_hash": phash})
        seen = set()
        for q, t, qs, qe, ts, te, fid, alen, cigar in alns[level]:
            # convertalis: query = representative, target = member; coordinates 1-based inclusive
            if assign[level].get(t) != q or t in seen:
                continue
            seen.add(t)
            aln_rows.append({"member_uid": bytes.fromhex(t), "level": level, "rep_uid": bytes.fromhex(q),
                             "qstart": int(ts) - 1, "qend": int(te), "tstart": int(qs) - 1, "tend": int(qe),
                             "pident": float(fid), "alnlen": int(alen), "cigar": cigar.encode()})
    mem_spec, aln_spec = table("cluster_member"), table("cluster_member_aln")
    mem_rows.sort(key=lambda r: (r["level"], r["cluster_uid"], r["member_uid"]))
    aln_rows.sort(key=lambda r: (r["level"], r["rep_uid"], r["member_uid"]))
    for name, rows, spec in (("cluster_member", mem_rows, mem_spec), ("cluster_member_aln", aln_rows, aln_spec)):
        p = root / f"l2/{name}/release={release}/part-0.parquet"
        p.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pylist(rows, schema=spec.schema), p, compression="zstd", row_group_size=122_880)
    counts = {lvl: len({r["cluster_uid"] for r in mem_rows if r["level"] == lvl}) for lvl in LEVELS}
    # ledger under the release's stage (rS0.1 -> S0): NR90 sees every protein, NR50 the NR90 representatives
    n_in = {"nr90": len(prots), "nr50": counts["nr90"]}
    ledger.write(root, ledger_stage(release), "cluster",
                 [ledger.row(lvl.encode(), "level", n_in[lvl], counts[lvl], n_in[lvl] - counts[lvl]) for lvl in LEVELS])
    return {"proteins": len(prots), "clusters": counts,
            "aligned": {lvl: sum(1 for r in aln_rows if r["level"] == lvl) for lvl in LEVELS}}
