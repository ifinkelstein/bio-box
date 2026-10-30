"""Release-time tokens for the tapes (PLAN.md §4.1 token namespace; data_layer_and_tests.md §1 L3).

Tapes carry stable keys (`protein_uid`); a release maps each slot to tokens at three levels:

    family  the canonical domain of the protein: the first library in `libraries` order with
            an included domain, and in it the domain with the highest bitscore (ties: model
            accession, then envelope start)
    nr50    the protein's NR50 cluster in this release
    nr90    the protein's NR90 cluster in this release

`token_id` is uint64: level (3 bits) | source (5 bits) | local (56 bits). PLAN.md §4.1 gives
32 bits with 24 local and widens to 64 when counts require it; NR90 clusters pass 2^24 by
S5, so the release line starts wide. Source is the library's 1-based position in
`libraries` for family tokens and 0 for clusters; local is the rank of the key within
(level, source), so ids are dense and reproducible but release-local (never a join key
outside the release).

Writes, under `<root>/l3/`:

    token_dict/release=<R>/part-<k>.parquet                          token_id -> level, source, key, n_slots
    tape_tokens/gene_call_version=<v>/release=<R>/part-<k>.parquet   per contig: one list per level
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox import layout
from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.registry.schemas import table

LEVELS = {"family": 1, "nr50": 2, "nr90": 3}  # 4 is reserved for subfamilies
CLUSTER_SOURCE = "mmseqs2"


def token_id(level: str, source: int, local: int) -> int:
    if not 0 <= source < 32 or not 0 <= local < 2**56:
        raise ValueError(f"token ({level}, {source}, {local}) out of range")
    return LEVELS[level] << 61 | source << 56 | local


def token_level(tid: int) -> str:
    return {v: k for k, v in LEVELS.items()}[tid >> 61]


def read_tapes(root: Path, gene_call_version: str, columns: list[str] | None = None) -> list[dict]:
    files = sorted((Path(root) / "l3/tapes").glob(f"gene_call_version={gene_call_version}/stage=*/part-*.parquet"))
    return [r for f in files for r in pq.read_table(f, columns=columns).to_pylist()]


def canonical_families(root: Path, libraries: list[str]) -> dict[bytes, tuple[str, str]]:
    """protein_uid -> (library, model_acc) under the precedence above."""
    best: dict[bytes, tuple] = {}
    for rank, lib in enumerate(libraries):
        for f in sorted((Path(root) / "l1/domains").glob(f"library={lib}/stage=*/part-*.parquet")):
            t = pq.read_table(f, columns=["protein_uid", "model_acc", "dom_bitscore", "env_start"]).to_pylist()
            for r in t:
                key = (rank, -r["dom_bitscore"], r["model_acc"], r["env_start"])
                if r["protein_uid"] not in best or key < best[r["protein_uid"]][0]:
                    best[r["protein_uid"]] = (key, lib, r["model_acc"])
    return {pu: (lib, acc) for pu, (_, lib, acc) in best.items()}


def libraries_in(root: Path) -> list[str]:
    return sorted({p.name.removeprefix("library=") for p in (Path(root) / "l1/domains").glob("library=*")})


def _write_parts(rows: list[dict], spec_name: str, base: Path, est_row: float, row_group: int, stats: list[str]) -> int:
    spec = table(spec_name)
    starts = layout.cut([est_row] * len(rows), layout.MAX_FILE_BYTES)
    base.mkdir(parents=True, exist_ok=True)
    for k, (a, b) in enumerate(zip(starts, starts[1:] + [len(rows)])):
        pq.write_table(pa.Table.from_pylist(rows[a:b], schema=spec.schema), base / f"part-{k}.parquet",
                       compression="zstd", row_group_size=row_group, write_statistics=stats)
    return len(starts)


def run(root: Path, release: str, gene_call_version: str = GENE_CALL_VERSION,
        libraries: list[str] | None = None) -> dict:
    root = Path(root)
    libraries = libraries_in(root) if libraries is None else libraries
    if len(libraries) > 31:
        raise ValueError("more than 31 libraries: widen the token source field")
    tapes = read_tapes(root, gene_call_version, ["contig_uid", "protein_uids"])
    fam = canonical_families(root, libraries)
    clusters: dict[str, dict[bytes, bytes]] = {"nr50": {}, "nr90": {}}
    p = root / f"l2/cluster_member/release={release}"
    for f in sorted(p.glob("part-*.parquet")) if p.exists() else []:
        for r in pq.read_table(f, columns=["member_uid", "level", "cluster_uid"]).to_pylist():
            if r["level"] in clusters:
                clusters[r["level"]][r["member_uid"]] = r["cluster_uid"]

    # keys per slot, then dense ids per (level, source)
    def slot_keys(pu: bytes | None) -> dict[str, tuple[int, str] | None]:
        out: dict[str, tuple[int, str] | None] = {}
        f = fam.get(pu) if pu is not None else None
        out["family"] = (libraries.index(f[0]) + 1, f[1]) if f else None
        for lvl in ("nr50", "nr90"):
            c = clusters[lvl].get(pu) if pu is not None else None
            out[lvl] = (0, c.hex()) if c else None
        return out

    per_contig = [(t["contig_uid"], [slot_keys(pu) for pu in t["protein_uids"]]) for t in tapes]
    uses: Counter = Counter()
    for _, slots in per_contig:
        for s in slots:
            for lvl, k in s.items():
                if k is not None:
                    uses[(lvl, *k)] += 1
    ids: dict[tuple, int] = {}
    for lvl in LEVELS:
        for src in sorted({s for (l, s, _) in uses if l == lvl}):
            for local, key in enumerate(sorted(k for (l, s, k) in uses if l == lvl and s == src)):
                ids[(lvl, src, key)] = token_id(lvl, src, local)

    source_name = {0: CLUSTER_SOURCE, **{i + 1: lib for i, lib in enumerate(libraries)}}
    dict_rows = sorted(({"token_id": tid, "level": lvl, "source": source_name[src], "key": key,
                         "n_slots": uses[(lvl, src, key)]} for (lvl, src, key), tid in ids.items()),
                       key=lambda r: r["token_id"])
    tok_rows = []
    for cu, slots in sorted(per_contig, key=lambda x: x[0]):
        row = {"contig_uid": cu}
        for lvl in LEVELS:
            row[f"{lvl}_tokens"] = [ids[(lvl, *s[lvl])] if s[lvl] is not None else None for s in slots]
        tok_rows.append(row)
    n_dict = _write_parts(dict_rows, "token_dict", root / f"l3/token_dict/release={release}", 60, 122_880,
                          ["token_id", "n_slots"])
    n_tok = _write_parts(tok_rows, "tape_tokens",
                         root / f"l3/tape_tokens/gene_call_version={gene_call_version}/release={release}",
                         100 + 30 * max(1, sum(len(s) for _, s in per_contig) // max(1, len(per_contig))),
                         8_192, ["contig_uid"])
    return {"tokens": {lvl: sum(1 for r in dict_rows if r["level"] == lvl) for lvl in LEVELS},
            "contigs": len(tok_rows), "files": n_dict + n_tok}
