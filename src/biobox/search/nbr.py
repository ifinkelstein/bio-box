"""Neighborhood fetch with explicit semantics (operation O-NBR; SEARCH_PLAN SD-08; S-09).

For each anchor slot, the slots in its window, with the anchor itself at offset 0:

    window      `10g`: slots within 10 slot offsets; `20kb` / `500bp`: slots whose gap to the
                anchor (max(starts) - min(ends), negative for overlaps) is at most that many bp
    strand      `any`, or `same`: only slots on the anchor's strand
    operon_gap  bp; when set, only slots joined to the anchor by a run of same-strand slots
                whose consecutive gaps are all <= operon_gap (`in_operon`)
    edge        linear contigs only. `truncate`: keep the anchor, cut its window at the
                contig end and set flag bit 0 when the cut removes window positions;
                `exclude:<bp>`: drop anchors within <bp> of a contig end (e.g. exclude:3kb);
                `require:<n>`: drop anchors with fewer than n slots on either side
    circular    `auto`: wrap on contigs whose tape says circular; `off`: treat every contig as linear

On a circular contig a slot is reached the shorter way round (by slot offset for gene
windows, forward on a tie; by gap for bp windows, forward on a tie), and flag bit 1 marks a
neighbor reached across the origin. `offset` is along the contig; `offset_rel` is along the
anchor's strand, positive downstream of it.

Slots are named by `feature_uid` (AGENTS §3.3, card B-16). An anchor is a protein occurrence,
identified by `(contig_uid, feature_uid)`; its ordinal is looked up in the one tape version the
fetch reads, and every row carries `anchor_feature_uid`, `feature_uid` and that `tape_version`
beside the ordinals, which are valid inside it only.

`summarize` (card W-06) counts which Pfam families and NR50/NR90 clusters sit in those windows, with DuckDB
over `fetch`'s rows joined to the store's `l1_domains`, `member_cluster` (a neighbor's clusters, by member;
card B-27), `cluster_member` (a cluster's members), `protein_seq`, `contig_dim` and `genome_input` (a contig's
genomes, card C-21c: counts by genome count every genome of a contig, counts by neighbour and locus never fan out).

Every read is of one release's files (card B-28): `occurrences` and `fetch` take `files`, a
`catalog.ReleaseFiles` (its manifest's files, else those of its `release_def` stages), so an anchor's
occurrence in a stage the release does not hold is not found.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from biobox.l1.genecall import GENE_CALL_VERSION
from biobox.l3.geometry import cyclic_offset, pair_geometry
from biobox.l3.tapes import tape_version as tape_version_of
from biobox.registry.schemas import table

FLAG_CENSORED = 1
FLAG_CROSSES_ORIGIN = 2
_BP = re.compile(r"^(\d+(?:\.\d+)?)(bp|kb|mb)$")
_UNITS = {"bp": 1, "kb": 1_000, "mb": 1_000_000}


class PolicyError(ValueError):
    pass


def parse_bp(text: str) -> int:
    m = _BP.match(str(text).strip().lower())
    if not m:
        raise PolicyError(f"{text!r}: expected a length like 500bp, 3kb")
    return int(float(m[1]) * _UNITS[m[2]])


@dataclass(frozen=True)
class Policy:
    window_genes: int | None = 10
    window_bp: int | None = None
    strand: str = "any"
    operon_gap: int | None = None
    edge: str = "truncate"  # truncate | exclude | require
    edge_bp: int = 0
    edge_n: int = 0
    circular: str = "auto"

    @classmethod
    def parse(cls, window: str = "10g", strand: str = "any", operon_gap: int | str | None = None,
              edge: str = "truncate", circular: str = "auto") -> "Policy":
        w = str(window).strip().lower()
        wg = wb = None
        if re.fullmatch(r"\d+g", w):
            wg = int(w[:-1])
        else:
            wb = parse_bp(w)
        if strand not in ("any", "same"):
            raise PolicyError(f"strand {strand!r}: expected any or same")
        if circular not in ("auto", "off"):
            raise PolicyError(f"circular {circular!r}: expected auto or off")
        og = None if operon_gap is None else (int(operon_gap) if str(operon_gap).isdigit() else parse_bp(operon_gap))
        e, ebp, en = edge, 0, 0
        if edge.startswith("exclude:"):
            e, ebp = "exclude", parse_bp(edge.split(":", 1)[1])
        elif edge.startswith("require:"):
            e, en = "require", int(edge.split(":", 1)[1])
        elif edge != "truncate":
            raise PolicyError(f"edge {edge!r}: expected truncate, exclude:<bp> or require:<n>")
        return cls(wg, wb, strand, og, e, ebp, en, circular)

    def hash(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:16]


def _operon(slots, i: int, n: int, circ: bool, gap: int, length: int) -> set[int]:
    """Slots joined to i by same-strand runs with every consecutive gap <= gap."""
    members = {i}
    strand = slots[i][2]
    for step in (1, -1):
        cur = i
        for _ in range(n - 1):
            nxt = cur + step
            if circ:
                nxt %= n
            elif not 0 <= nxt < n:
                break
            if nxt in members or slots[nxt][2] != strand:
                break
            _, d, _ = pair_geometry(slots, cur, nxt, step, length, circ)
            if d > gap:
                break
            members.add(nxt)
            cur = nxt
    return members


def _candidates(slots, i: int, n: int, circ: bool, p: Policy, length: int):
    """(j, offset, topo, delta, crosses) for every slot in i's window, before strand/operon filters."""
    out = []
    if p.window_genes is not None:
        k = p.window_genes
        for j in range(n):
            if j == i:
                continue
            d = cyclic_offset(i, j, n) if circ else j - i
            if abs(d) <= k:
                out.append((j, d, *pair_geometry(slots, i, j, d, length, circ)))
        return out
    for j in range(n):
        if j == i:
            continue
        if circ:
            fwd, back = (j - i) % n, -((i - j) % n)
            a = pair_geometry(slots, i, j, fwd, length, circ)
            b = pair_geometry(slots, i, j, back, length, circ)
            d, g = (back, b) if b[1] < a[1] else (fwd, a)
        else:
            d = j - i
            g = pair_geometry(slots, i, j, d, length, circ)
        if g[1] <= p.window_bp:
            out.append((j, d, *g))
    return out


def neighborhood(tape: dict, i: int, policy: Policy) -> list[dict] | None:
    """Rows for anchor slot i of one tape; None when the edge policy drops the anchor."""
    slots = list(zip(tape["start_bps"], tape["end_bps"], tape["strands"]))
    n, length = len(slots), tape["contig_len"]
    circ = bool(tape["is_circular"]) and policy.circular == "auto"
    s, e, strand = slots[i]
    censored = False
    if not circ:
        if policy.edge == "exclude" and min(s, length - e) < policy.edge_bp:
            return None
        if policy.edge == "require" and (i < policy.edge_n or n - 1 - i < policy.edge_n):
            return None
        if policy.window_genes is not None:
            censored = i < policy.window_genes or n - 1 - i < policy.window_genes
        else:
            censored = s - policy.window_bp < 0 or e + policy.window_bp > length
    operon = _operon(slots, i, n, circ, policy.operon_gap, length) if policy.operon_gap is not None else None
    base = FLAG_CENSORED if censored else 0
    rows = [_row(tape, i, i, 0, strand, None, None, operon, base, policy)]
    for j, d, topo, delta, crosses in _candidates(slots, i, n, circ, policy, length):
        if policy.strand == "same" and slots[j][2] != strand:
            continue
        if operon is not None and j not in operon:
            continue
        rows.append(_row(tape, i, j, d, strand, topo, delta, operon,
                         base | (FLAG_CROSSES_ORIGIN if crosses else 0), policy))
    rows.sort(key=lambda r: r["offset"])
    return rows


def _row(tape, i, j, d, anchor_strand, topo, delta, operon, flags, policy) -> dict:
    return {"anchor_protein_uid": tape["protein_uids"][i], "contig_uid": tape["contig_uid"], "anchor_ordinal": i,
            "ordinal": j, "offset": d, "offset_rel": d * anchor_strand, "protein_uid": tape["protein_uids"][j],
            "start": tape["start_bps"][j], "end": tape["end_bps"][j], "strand": tape["strands"][j],
            "delta_bp": delta, "junction_topology": topo, "same_strand": tape["strands"][j] == anchor_strand,
            "in_operon": None if operon is None else j in operon, "flags": flags, "policy_hash": policy.hash(),
            "anchor_feature_uid": tape["feature_uids"][i], "feature_uid": tape["feature_uids"][j],
            "tape_version": tape.get("tape_version")}


def resolve_tape_version(gene_call_version: str = GENE_CALL_VERSION, tape_version: str | None = None,
                         producers: list[str] | None = None) -> str:
    """The tape version a read uses: `tape_version` if given, else that of `producers` (default: the
    single gene-call producer, card B-07)."""
    if tape_version is not None:
        return tape_version
    return tape_version_of(producers if producers is not None else [gene_call_version])


def _tapes(files, tv: str, contigs: set[bytes]) -> dict[bytes, dict]:
    """The tapes of tape version tv for the contigs in the release's files, each with its `tape_version`."""
    out = {}
    for f in files.table("contig_tapes", tape_version=tv):
        t = pq.read_table(f, filters=[("contig_uid", "in", sorted(contigs))]) if contigs else None
        for r in (t.to_pylist() if t is not None else []):
            out[r["contig_uid"]] = {**r, "tape_version": tv}
    return out


def _ordinal(tapes: dict[bytes, dict], index: dict[bytes, dict], cu: bytes, fu: bytes, tv: str) -> int:
    """The ordinal of feature fu of contig cu inside tape version tv."""
    if cu not in tapes:
        raise ValueError(f"contig {cu.hex()} has no tape in tape_version {tv}")
    if cu not in index:
        index[cu] = {u: k for k, u in enumerate(tapes[cu]["feature_uids"])}
    if fu not in index[cu]:
        raise ValueError(f"feature {fu.hex()} of contig {cu.hex()} is not a slot of tape_version {tv}")
    return index[cu][fu]


def _release_files(root: Path, files):
    if files is None or not hasattr(files, "table"):
        raise TypeError("nbr reads one release's files: pass files=catalog.ReleaseFiles.of(root, release)")
    if Path(files.root).resolve() != Path(root).resolve():
        raise ValueError(f"files of release {files.release} are under {files.root}, not {root}")
    return files


def occurrences(root: Path, protein_uids: set[bytes], gene_call_version: str = GENE_CALL_VERSION,
                tape_version: str | None = None, *, files) -> list[dict]:
    """protein_occurrence rows of the proteins in the release's `files` (`catalog.ReleaseFiles`). With
    `tape_version`, each row also carries `ordinal`, the position of its feature in that tape version, and
    `tape_version` (an ordinal is valid inside one tape version only, so it never travels without them;
    AGENTS §3.3)."""
    files = _release_files(root, files)
    if not protein_uids:
        return []
    rows = [r for f in files.table("protein_occurrence", gene_call_version=gene_call_version)
            for r in pq.read_table(f, filters=[("protein_uid", "in", sorted(protein_uids))]).to_pylist()]
    if tape_version is None:
        return rows
    tapes, index = _tapes(files, tape_version, {r["contig_uid"] for r in rows}), {}
    return [{**r, "ordinal": _ordinal(tapes, index, r["contig_uid"], r["feature_uid"], tape_version),
             "tape_version": tape_version} for r in rows]


def fetch(root: Path, protein_uids: set[bytes], policy: Policy, gene_call_version: str = GENE_CALL_VERSION,
          tape_version: str | None = None, producers: list[str] | None = None, *, files) -> tuple[list[dict], dict]:
    """Neighborhoods of every occurrence of the anchor proteins in the release's `files`
    (`catalog.ReleaseFiles`, card B-28); (rows, counts).

    The tapes read are those of `tape_version` (default: the tape version of `producers`, which
    defaults to the gene-call producer). Each anchor is found on its tape by its feature_uid."""
    tv = resolve_tape_version(gene_call_version, tape_version, producers)
    occ = occurrences(root, protein_uids, gene_call_version, files=files)
    tapes, index = _tapes(files, tv, {r["contig_uid"] for r in occ}), {}
    rows, dropped = [], 0
    for r in sorted(occ, key=lambda r: (r["contig_uid"], r["feature_uid"])):
        i = _ordinal(tapes, index, r["contig_uid"], r["feature_uid"], tv)
        got = neighborhood(tapes[r["contig_uid"]], i, policy)
        if got is None:
            dropped += 1
        else:
            rows += got
    rows.sort(key=lambda r: (r["anchor_protein_uid"], r["contig_uid"], r["anchor_ordinal"], r["offset"]))
    return rows, {"anchors": len(occ), "anchors_dropped_by_edge": dropped, "rows": len(rows), "tape_version": tv}


def write(rows: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("neighborhoods").schema), path, compression="zstd")
    return path


# --- neighbor summary (card W-06) ------------------------------------------------------------------------

SUMMARY_LEVELS = ("pfam", "nr50", "nr90")
SUMMARY_TOP = 50

# the relations the aggregation reads (summarize: views over the store's files; summarize_rows: given rows)
_SUMMARY_SCHEMAS = {
    "domains": pa.schema([("protein_uid", pa.binary(16)), ("library", pa.string()), ("model_acc", pa.string()),
                          ("model_name", pa.string())]),
    "clusters": pa.schema([("member_uid", pa.binary(16)), ("level", pa.string()), ("cluster_uid", pa.binary(16)),
                           ("rep_uid", pa.binary(16))]),
    # the same rows, looked up by member (summarize: `member_cluster`, card B-27; summarize_rows: `clusters`)
    "members": pa.schema([("member_uid", pa.binary(16)), ("level", pa.string()), ("cluster_uid", pa.binary(16)),
                          ("rep_uid", pa.binary(16))]),
    "protein_len": pa.schema([("protein_uid", pa.binary(16)), ("length", pa.uint32())]),
    # a contig's row in the release's `contig_dim`: one per contig, its genome the least of its genomes (DEC-27)
    "genome_of": pa.schema([("contig_uid", pa.binary(16)), ("genome_uid", pa.binary(16))]),
    # every genome that holds a contig (card C-21c): the `genome_input` rows of the release's stages, as the catalog's
    # `genome_member`; a contig in several genomes has a row for each
    "genome_member": pa.schema([("contig_uid", pa.binary(16)), ("genome_uid", pa.binary(16))]),
}


def _check_levels(by) -> tuple[str, ...]:
    by = tuple(by)
    bad = [b for b in by if b not in SUMMARY_LEVELS]
    if bad or not by:
        raise ValueError(f"summary levels {bad or list(by)}: expected some of {', '.join(SUMMARY_LEVELS)}")
    return tuple(dict.fromkeys(by))


def _level_sql(level: str) -> tuple[str, str]:
    """(keyed, attrs) for one level. `keyed`: one row per (neighbor row, key) with the neighbor's columns;
    `attrs`: one row per key with the key's descriptive columns."""
    if level == "pfam":
        keyed = """SELECT n.*, d.library || ':' || d.model_acc AS key
                   FROM _nb n JOIN (SELECT DISTINCT protein_uid, library, model_acc FROM domains
                                    WHERE protein_uid IN (SELECT protein_uid FROM _nb)) d USING (protein_uid)"""
        attrs = """SELECT library || ':' || model_acc AS key, library, model_acc, min(model_name) AS model_name,
                          NULL::VARCHAR AS cluster_uid, NULL::VARCHAR AS rep_uid, NULL::BIGINT AS rep_length
                   FROM domains WHERE protein_uid IN (SELECT protein_uid FROM _nb) GROUP BY library, model_acc"""
        return keyed, attrs
    # a neighbor's cluster is looked up by member (`members`, sorted by level and member_uid); the clusters' rows by
    # cluster (`clusters`, sorted by level and cluster_uid)
    cl = f"""(SELECT DISTINCT member_uid, cluster_uid, rep_uid FROM members
              WHERE level = '{level}' AND member_uid IN (SELECT protein_uid FROM _nb))"""
    keyed = f"""SELECT n.*, lower(hex(c.cluster_uid)) AS key FROM _nb n JOIN {cl} c ON c.member_uid = n.protein_uid"""
    # the representative: one per cluster in a real clustering; the least rep_uid if a mapping gives several
    attrs = f"""WITH reps AS (SELECT cluster_uid, min(rep_uid) AS rep_uid FROM clusters
                              WHERE level = '{level}' AND cluster_uid IN (SELECT cluster_uid FROM {cl})
                              GROUP BY cluster_uid),
                     names AS (SELECT m.cluster_uid, d.model_name, count(DISTINCT m.member_uid) AS n
                               FROM clusters m JOIN domains d ON d.protein_uid = m.member_uid
                               WHERE m.level = '{level}' AND m.cluster_uid IN (SELECT cluster_uid FROM reps)
                                 AND d.protein_uid IN (SELECT member_uid FROM clusters WHERE level = '{level}'
                                                       AND cluster_uid IN (SELECT cluster_uid FROM reps))
                               GROUP BY m.cluster_uid, d.model_name),
                     top AS (SELECT cluster_uid, first(model_name ORDER BY n DESC, model_name) AS model_name
                             FROM names GROUP BY cluster_uid),
                     lens AS (SELECT protein_uid, max(length) AS length FROM protein_len
                              WHERE protein_uid IN (SELECT rep_uid FROM reps) GROUP BY protein_uid)
                SELECT lower(hex(r.cluster_uid)) AS key, NULL::VARCHAR AS library, NULL::VARCHAR AS model_acc,
                       coalesce(t.model_name, '') AS model_name, lower(hex(r.cluster_uid)) AS cluster_uid,
                       lower(hex(r.rep_uid)) AS rep_uid, l.length::BIGINT AS rep_length
                FROM reps r LEFT JOIN top t USING (cluster_uid) LEFT JOIN lens l ON l.protein_uid = r.rep_uid"""
    return keyed, attrs


# Counts by neighbour and locus read the neighbour rows only; the genomes of a key are counted apart, from its
# distinct contigs (`{genomes}`), so a contig in several genomes never repeats a neighbour (card C-21c, DEC-27)
_STATS = """WITH s AS (SELECT key, count(DISTINCT (contig_uid, anchor_feature_uid)) AS n_loci, count(*) AS n_neighbors,
                              count(*) FILTER (WHERE same_strand) AS n_same
                       FROM {k} GROUP BY key),
                 g AS (SELECT c.key, {genomes} AS n_genomes
                       FROM (SELECT DISTINCT key, contig_uid FROM {k}) c {join} GROUP BY c.key),
                 h AS (SELECT key, offset_rel, count(*) AS n FROM {k} GROUP BY key, offset_rel),
                 o AS (SELECT key, list(struct_pack(o := offset_rel, n := n) ORDER BY offset_rel) AS offsets,
                              first(offset_rel ORDER BY n DESC, abs(offset_rel), offset_rel DESC) AS modal_offset
                       FROM h GROUP BY key)
            SELECT s.key, s.n_loci, s.n_neighbors, g.n_genomes, s.n_same, o.offsets, o.modal_offset
            FROM s JOIN g USING (key) JOIN o USING (key) ORDER BY s.n_loci DESC, s.key"""


def _aggregate(con, by: tuple[str, ...], top_n: int | None) -> dict:
    """The summary over the relations `nbr_rows`, `domains`, `clusters`, `members`, `protein_len`, `genome_of` and
    `genome_member` of `con` (columns as in `_SUMMARY_SCHEMAS`; `nbr_rows` as the `neighborhoods` table). Every count
    is SQL.

    Genomes (card C-21c, DEC-27): the basis is `genome_uid` when every anchor locus's contig has a `genome_of` row (the
    release's `contig_dim`), else `contig_uid`, as before the card. On the genome basis a contig's genomes are its
    `genome_of` genome and every `genome_member` genome (the release's `genome_input`), so a key on a contig of two
    genomes counts both. Genome counts are distinct genomes over a key's distinct contigs, computed apart from the
    neighbour rows: `n_loci`, `n_neighbors`, the offsets and the header's locus counts never join the membership."""
    con.execute("""CREATE TEMP TABLE _loci AS
                   SELECT contig_uid, anchor_feature_uid, any_value(anchor_protein_uid) AS anchor_protein_uid,
                          bool_or("offset" = 0 AND (flags & 1) = 1) AS censored
                   FROM nbr_rows GROUP BY contig_uid, anchor_feature_uid""")
    # each relation filtered by the loci's contigs before the union, so each scan prunes to their row groups
    con.execute("""CREATE TEMP TABLE _gen AS SELECT DISTINCT contig_uid, genome_uid FROM (
                     SELECT contig_uid, genome_uid FROM genome_of WHERE contig_uid IN (SELECT contig_uid FROM _loci)
                     UNION ALL
                     SELECT contig_uid, genome_uid FROM genome_member WHERE contig_uid IN (SELECT contig_uid FROM _loci))""")
    n_contigs, with_dim = con.execute(
        "SELECT (SELECT count(DISTINCT contig_uid) FROM _loci), "
        "(SELECT count(DISTINCT contig_uid) FROM genome_of WHERE contig_uid IN (SELECT contig_uid FROM _loci))"
    ).fetchall()[0]
    basis = "genome_uid" if n_contigs and with_dim == n_contigs else "contig_uid"
    genomes, join = (("count(DISTINCT m.genome_uid)", "LEFT JOIN _gen m ON m.contig_uid = c.contig_uid")
                     if basis == "genome_uid" else ("count(DISTINCT c.contig_uid)", ""))
    con.execute("""CREATE TEMP TABLE _nb AS
                   SELECT n.contig_uid, n.anchor_feature_uid, n.feature_uid, n.protein_uid, n.offset_rel, n.same_strand
                   FROM nbr_rows n WHERE n."offset" <> 0""")
    n_loci, n_cens, n_prot = con.execute(
        "SELECT count(*), count(*) FILTER (WHERE censored), count(DISTINCT anchor_protein_uid) FROM _loci").fetchall()[0]
    n_genomes = con.execute(f"SELECT {genomes} FROM (SELECT DISTINCT contig_uid FROM _loci) c {join}").fetchall()[0][0]
    n_nb = con.execute("SELECT count(*) FROM _nb").fetchall()[0][0]
    rows, n_keys = [], {}

    def record(level: str, r: tuple, attrs: dict, unannotated: bool) -> dict:
        key, nl, nn, ng, same, offsets, modal = r
        return {"level": level, "key": key, "unannotated": unannotated, **attrs,
                "n_loci": nl, "share_of_loci": round(nl / n_loci, 4) if n_loci else None, "n_neighbors": nn,
                "n_genomes": ng, "genome_basis": basis,
                "offsets": {str(x["o"]): x["n"] for x in offsets or []}, "modal_offset": modal,
                "same_strand_share": round(same / nn, 4) if nn else None}

    empty_attrs = {"library": None, "model_acc": None, "model_name": None, "cluster_uid": None, "rep_uid": None,
                   "rep_length": None}
    for level in by:
        keyed, attrs_sql = _level_sql(level)
        con.execute(f"CREATE OR REPLACE TEMP TABLE _k AS {keyed}")
        con.execute("CREATE OR REPLACE TEMP TABLE _u AS SELECT n.*, '' AS key FROM _nb n "
                    "WHERE NOT EXISTS (SELECT 1 FROM _k k WHERE k.contig_uid = n.contig_uid "
                    "AND k.anchor_feature_uid = n.anchor_feature_uid AND k.feature_uid = n.feature_uid)")
        stats = con.execute(_STATS.format(k="_k", genomes=genomes, join=join)).fetchall()
        n_keys[level] = len(stats)
        shown = stats[:top_n] if top_n else stats
        cols = ["key", "library", "model_acc", "model_name", "cluster_uid", "rep_uid", "rep_length"]
        attrs = {r[0]: dict(zip(cols[1:], r[1:]))
                 for r in con.execute(f"SELECT {', '.join(cols)} FROM ({attrs_sql}) a "
                                      "WHERE key IN (SELECT unnest(?::VARCHAR[]))", [[r[0] for r in shown]]).fetchall()}
        rows += [record(level, r, attrs[r[0]], False) for r in shown]
        un = con.execute(_STATS.format(k="_u", genomes=genomes, join=join)).fetchall() or [("", 0, 0, 0, 0, [], None)]
        rows.append(record(level, un[0], empty_attrs, True))
    return {"header": {"levels": list(by), "top_n": top_n or None, "n_anchor_proteins": n_prot,
                       "n_anchor_loci": n_loci, "n_contigs": n_contigs, "n_anchor_genomes": n_genomes,
                       "genome_basis": basis, "n_censored_loci": n_cens, "n_neighbors": n_nb, "n_keys": n_keys},
            "rows": rows}


def _connect(threads: int | None = None):
    import duckdb

    return duckdb.connect(config={"threads": threads} if threads else {})


def summarize_rows(rows: list[dict] | pa.Table, by=("pfam", "nr50"), *, domains=None, clusters=None,
                   protein_len=None, genome_of=None, genome_member=None, top_n: int | None = SUMMARY_TOP) -> dict:
    """The neighbor summary of neighborhood rows (as `fetch` gives them) over the given annotations, each a
    list of dicts or an Arrow table with the columns of `_SUMMARY_SCHEMAS` (None: empty); `clusters` serves the
    lookups by member as well. `genome_of` is a contig's `contig_dim` genome, `genome_member` every genome that
    holds it (card C-21c; `_aggregate`)."""
    by = _check_levels(by)
    con = _connect()
    try:
        con.register("nbr_rows", rows if isinstance(rows, pa.Table)
                     else pa.Table.from_pylist(rows, schema=table("neighborhoods").schema))
        for name, data in (("domains", domains), ("clusters", clusters), ("protein_len", protein_len),
                           ("genome_of", genome_of), ("genome_member", genome_member)):
            t = data if isinstance(data, pa.Table) else pa.Table.from_pylist(data or [], schema=_SUMMARY_SCHEMAS[name])
            con.register(name, t)
        con.execute("CREATE VIEW members AS SELECT member_uid, level, cluster_uid, rep_uid FROM clusters")
        return _aggregate(con, by, top_n)
    finally:
        con.close()


def _sql_str(s: str) -> str:
    return "'" + str(s).replace("'", "''") + "'"


_NOTED: set | None = None  # the notes this command has printed (see `noted_once`); None outside one


@contextlib.contextmanager
def noted_once():
    """Within it, `slow_path_note` prints each (release, lookup, source) once (card W-20): `bio-box search hmm` opens
    nine stores and printed the same note nine times. Library calls outside it print every time. Re-entrant."""
    global _NOTED
    if _NOTED is not None:
        yield
        return
    _NOTED = set()
    try:
        yield
    finally:
        _NOTED = None


def slow_path_note(release: str | None, lookup: str, source: str) -> None:
    """One stderr line (card W-17): the release has no `lookup` files (it was built before card B-27), so the
    reader answers from `source`, which cannot prune a lookup by this key. Never an empty view in its place.
    Once per command inside `noted_once` (card W-20)."""
    if _NOTED is not None:
        if (release, lookup, source) in _NOTED:
            return
        _NOTED.add((release, lookup, source))
    print(f"note: release {release} has no {lookup} files (built before card B-27); reading {source} instead "
          f"(the slow path: every row group is scanned)", file=sys.stderr)


def _store_files(root: Path, name: str, want: dict[str, str], manifest) -> list[tuple[Path, dict]]:
    """(path, partitions) of one table's files whose partitions include `want`: from a `catalog.ReleaseFiles`
    (card B-28), or from a release manifest dict (AGENTS §3.9). Readers never list the registry location (no
    glob over stages): since card W-12 every caller passes its release's files, so None is a TypeError."""
    if hasattr(manifest, "entries"):  # catalog.ReleaseFiles
        return sorted(manifest.entries(name, **want), key=lambda x: str(x[0]))
    if manifest is None:
        raise TypeError(f"{name}: a reader takes its release's files (catalog.ReleaseFiles or a manifest), never a glob")
    found = [(root / f["path"], f["partitions"]) for f in manifest["files"] if f["table"] == name]
    return sorted(((p, parts) for p, parts in found if all(parts.get(k) == v for k, v in want.items())),
                  key=lambda x: str(x[0]))


def _store_view(con, name: str, files: list[tuple[Path, dict]], select: str, partition_cols: tuple[str, ...] = ()) -> None:
    """A view `name` over explicit files, with partition values as literal columns; typed and empty without files."""
    schema = _SUMMARY_SCHEMAS[name]
    if not files:
        from biobox.catalog import duck_type

        cols = ", ".join(f'CAST(NULL AS {duck_type(f.type)}) AS "{f.name}"' for f in schema)
        con.execute(f"CREATE VIEW {name} AS SELECT {cols} WHERE false")
        return
    parts = [f"SELECT {select}{''.join(f', {_sql_str(p[c])} AS {c}' for c in partition_cols)} "
             f"FROM read_parquet({_sql_str(str(f))}, hive_partitioning = false)" for f, p in files]
    con.execute(f"CREATE VIEW {name} AS SELECT {', '.join(f.name for f in schema)} FROM ({' UNION ALL '.join(parts)})")


def summarize(root: Path, anchors: set[bytes], policy: Policy, by=("pfam", "nr50"), *, release: str,
              gene_call_version: str = GENE_CALL_VERSION, tape_version: str | None = None,
              top_n: int | None = SUMMARY_TOP, manifest: dict | None = None, threads: int | None = None) -> dict:
    """Which families (`pfam`: every library of `l1_domains`) and clusters (`nr50`, `nr90`: `cluster_member` of
    `release`) sit in the windows of the anchors' loci, at which strand-relative offsets, in what share of the
    loci and how many genomes (when every anchor locus's contig has a `contig_dim` row of `release`: the genomes of
    those contigs, each contig's `contig_dim` genome and every genome the release's `genome_input` lists it in, card
    C-21c; else contigs). The rows are `fetch`'s; the anchor
    itself (offset 0) is no neighbor. A neighbor with no key at a level (no domain, no cluster, a non-CDS slot)
    counts in that level's `unannotated` row. `top_n` (0 or None: all) caps the keyed rows per level after the
    totals. File lists come from `manifest` (default: `manifests/release=<release>.json` if present; without
    one, from the release's `release_def` stages).

    Per key: `n_loci` (anchor loci with such a neighbor), `share_of_loci` (of all anchor loci), `n_neighbors`
    (neighbor rows: a gene in two anchors' windows counts twice; a gene with several domains of one family
    once), `n_genomes`, `offsets` ({offset_rel: rows}, positive downstream along the anchor's strand),
    `modal_offset` (the commonest; ties go to the one nearer the anchor, then to the downstream one) and
    `same_strand_share`. A key at `nr50`/`nr90` is the cluster uid, with its representative (the least rep_uid
    if a mapping gives several), the representative's length and the `model_name` most members carry ("" if
    none; ties by name). Rows by `n_loci` descending, then key; shares rounded to 4 decimals. The header holds
    the run-level facts, among them `n_censored_loci` (loci whose window a linear contig end cut).
    Without a manifest the files are those of the stages of the release's `release_def` (card B-28)."""
    from biobox import catalog

    by = _check_levels(by)
    root = Path(root)
    files = catalog.ReleaseFiles.of(root, release, manifest)
    rows, counts = fetch(root, set(anchors), policy, gene_call_version, tape_version=tape_version, files=files)
    con = _connect(threads)
    try:
        con.register("nbr_rows", pa.Table.from_pylist(rows, schema=table("neighborhoods").schema))
        views = [("domains", "l1_domains", {}, "protein_uid, model_acc, model_name", ("library",)),
                 ("clusters", "cluster_member", {"release": release}, "member_uid, level, cluster_uid, rep_uid", ()),
                 ("members", "member_cluster", {"release": release}, "member_uid, level, cluster_uid, rep_uid", ()),
                 ("protein_len", "protein_seq", {"gene_call_version": gene_call_version}, "protein_uid, length", ()),
                 ("genome_of", "contig_dim", {"release": release}, "contig_uid, genome_uid", ()),
                 # card C-21c: the release's membership, the catalog's `genome_member` (its stages' genome_input)
                 ("genome_member", "genome_input", {}, "contig_uid, genome_uid", ())]
        largest = 0
        for name, spec_name, want, select, parts in views:
            found = _store_files(root, spec_name, want, files)
            if name == "members" and not found:  # a release built before card B-27 (card W-17)
                found = _store_files(root, "cluster_member", want, files)
                if found:
                    slow_path_note(release, "member_cluster", "cluster_member")
            _store_view(con, name, found, select, parts)
            largest = max(largest, catalog.row_groups_of(p for p, _ in found))
        catalog.set_join_key_filter(con, largest)  # keys from _nb prune to the row groups that hold them (card B-27)
        out = _aggregate(con, by, top_n)
    finally:
        con.close()
    out["header"] = {"release": release, "tape_version": counts["tape_version"], "policy": asdict(policy),
                     "policy_hash": policy.hash(), "n_anchor_proteins_given": len(set(anchors)),
                     "anchors_dropped_by_edge": counts["anchors_dropped_by_edge"], **out["header"]}
    return out
