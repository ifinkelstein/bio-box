"""Ring R2: a contig's genomes on a store split into 20 stages (card C-21b; AGENTS §4.3, §4.4; design note
`planning/design/C-21_shared_contigs.md` §2a and §7 item 8).

`genome_member` is a projection of the release's stages' `genome_input` files, each sorted by `contig_uid`, with no
aggregate: DuckDB pushes no join filter through one (`catalog.label_view_sql`). A lookup of a contig's genomes then
reads, in each stage file, only the row groups whose footer `contig_uid` range holds the key: one, or two where the
key's rows straddle a row-group boundary. That is the bound tested here, for `genomes_of` with a literal key and a
parameter, and for 10 and 100 keys from a joined table (a join to `genome_member`, and `genomes_of` given the table's
keys as one list, as `tapes_of` takes them), against the row groups the files' footers say may hold the keys. The
reads per key grow with the number of stages (one file each), not with the release's rows.

The store: 20 stages of synthetic membership written by the repo's writers (`write_genome_input`; `io.write_parts`
for `contig_index`) with row groups scaled down to ROW_CAP rows (test_scale_layout's method: `io.write_table` patched
test side; the real caps are AGENTS §4.6), so each stage's `genome_input` file has hundreds of row groups. Contigs are
shared within a stage and across stages, and some genomes are listed again by the next stage with another assembly
quality. `dim.run` builds the release's `contig_dim` and `genome_dim`, once with `release_def` in stage order and once
reversed. The oracle is a brute-force pass in Python over the rows the generator made, not the code under test.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from biobox import catalog, dim, ids
from biobox import io as bio
from biobox.l0seq.inputs import GenomeInput, write_genome_input
from biobox.registry.schemas import table

pytestmark = [pytest.mark.r2, pytest.mark.covers("LY-06", "D-01")]

N_STAGES = 20
GENOMES = 1_500  # per stage
ROW_CAP = 16  # rows per row group of every table written here (scaled down; AGENTS §4.6 caps at 122,880)
MAX_RG_PER_KEY_AND_FILE = 2  # AGENTS §4.4, per stage file: the key's row group, or two where its rows straddle a cut
SEED = 2121
QUALITIES = ("complete", "chromosome", "scaffold")
STAGES = [f"M{k:02d}" for k in range(N_STAGES)]
REL, REL_REV = "rC21b.m20", "rC21b.m20rev"
GI = "/l0seq/genome_input/"


def _generate() -> tuple[dict, dict]:
    """contigs {uid: (owning stage, length, circular)} and genomes {stage: [(contig uids, assembly quality)]}."""
    rng = np.random.default_rng(SEED)
    contigs: dict[bytes, tuple[str, int, bool]] = {}
    genomes: dict[str, list[tuple[tuple[bytes, ...], str]]] = {}
    earlier: list[bytes] = []
    for k, s in enumerate(STAGES):
        gs: list[tuple[tuple[bytes, ...], str]] = []
        fresh: list[bytes] = []
        for _ in range(GENOMES):
            r = rng.random()
            if k and r < 0.02:  # a genome of the stage before, listed again with an assembly quality drawn again
                prev = genomes[STAGES[k - 1]][int(rng.integers(len(genomes[STAGES[k - 1]])))][0]
                if all(prev != g[0] for g in gs):
                    gs.append((prev, QUALITIES[int(rng.integers(3))]))
                continue
            cs = []
            for _ in range(int(rng.integers(1, 4))):
                u = rng.bytes(16)
                contigs[u] = (s, int(rng.integers(1_000, 300_000)), bool(rng.integers(0, 2)))
                cs.append(u)
                fresh.append(u)
            if earlier and r < 0.10:  # a contig an earlier stage stored (shared across stages)
                cs.append(earlier[int(rng.integers(len(earlier)))])
            elif gs and r < 0.18:  # the first contig of a genome of this stage (shared within the stage)
                other = gs[int(rng.integers(len(gs)))][0][0]
                if other not in cs:
                    cs.append(other)
            gs.append((tuple(cs), QUALITIES[int(rng.integers(3))]))
        genomes[s] = gs
        earlier += fresh
    return contigs, genomes


def _capped(real):
    def write_table(batches, spec, path, *, row_group_rows=bio.ROW_GROUP_ROWS, **kw):
        return real(batches, spec, path, row_group_rows=min(row_group_rows, ROW_CAP), **kw)
    return write_table


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    d = tmp_path_factory.mktemp("c21b_m20")
    root, out_rev = d / "root", d / "rev"
    contigs, genomes = _generate()
    spec = table("contig_index")
    mp = pytest.MonkeyPatch()
    mp.setattr(bio, "write_table", _capped(bio.write_table))
    try:
        for s in STAGES:
            own = sorted((u, hashlib.sha256(u).digest(), n) for u, (o, n, _) in contigs.items() if o == s)
            bio.write_parts(pa.table({"contig_uid": pa.array([x[0] for x in own], pa.binary(16)),
                                      "raw_digest": pa.array([x[1] for x in own], pa.binary(32)),
                                      "contig_len": pa.array([x[2] for x in own], pa.int64()),
                                      "is_default_variant": pa.array([True] * len(own))}, schema=spec.schema),
                            spec, root / spec.location / f"stage={s}", 36)
            write_genome_input(root, s, [GenomeInput(cs, 11, kind="prokaryote",
                                                     circular=tuple(c for c in cs if contigs[c][2]), assembly_quality=q)
                                         for cs, q in genomes[s]])
        res = dim.run(root, REL, STAGES, tmp=d)
        res_rev = dim.run(root, REL_REV, STAGES[::-1], tmp=d, out=out_rev)
    finally:
        mp.undo()
    man = catalog.build_manifest(root, REL, STAGES)
    con = catalog.open_catalog(catalog.generate_catalog(man, d / "c.duckdb"), threads=4)
    yield {"root": root, "rev": out_rev, "contigs": contigs, "genomes": genomes, "res": res, "res_rev": res_rev,
           "man": man, "con": con}
    con.close()


def _rows(store) -> set[tuple[bytes, bytes, str]]:
    """Every (contig, genome, stage) the generator wrote: the oracle of `genome_member`."""
    return {(c, ids.genome_uid(list(cs)), s) for s in STAGES for cs, _ in store["genomes"][s] for c in cs}


def _gi_files(store) -> list[Path]:
    return [store["root"] / f["path"] for f in store["man"]["files"] if f["table"] == "genome_input"]


def _footer_row_groups(files: list[Path], keys: list[bytes]) -> tuple[int, int]:
    """(row groups whose footer contig_uid range holds one of `keys`, the most any one key has in one file)."""
    total, most = 0, 0
    ks = sorted(keys)
    for p in files:
        md = pq.read_metadata(p)
        per_key: dict[bytes, int] = {}
        for g in range(md.num_row_groups):
            st = md.row_group(g).column(0).statistics  # contig_uid, the first column and sort key
            held = [k for k in ks if st.min <= k <= st.max]
            total += bool(held)
            for k in held:
                per_key[k] = per_key.get(k, 0) + 1
        most = max([most, *per_key.values()])
    return total, most


def _gi_scanned(prof: dict) -> tuple[int, int]:
    """(row groups the query's genome_input scans read, of how many)."""
    scans = [s for s in prof["scans"] if any(GI in f for f in s["files"])]
    assert scans, prof["scans"]
    return sum(s["row_groups_scanned"] for s in scans), sum(s["row_groups_total"] for s in scans)


def test_the_store_has_twenty_stage_files_of_many_row_groups_and_shares_contigs(store):
    files = _gi_files(store)
    assert len(files) == N_STAGES
    assert min(pq.read_metadata(p).num_row_groups for p in files) >= 100
    rows = _rows(store)
    by_contig: dict[bytes, set] = {}
    for c, g, s in rows:
        by_contig.setdefault(c, set()).add((g, s))
    stages_of = {c: {s for _, s in v} for c, v in by_contig.items()}
    assert sum(len(v) > 1 for v in stages_of.values()) > 100  # listed by several stages
    assert sum(len({g for g, s in v if s == store["contigs"][c][0]}) > 1 for c, v in by_contig.items()) > 100
    relisted = {ids.genome_uid(list(cs)) for s in STAGES for cs, _ in store["genomes"][s]}
    assert len(relisted) < sum(len(store["genomes"][s]) for s in STAGES)  # some genomes are in two stages


@pytest.mark.covers("C-01")
def test_contig_dim_and_genome_dim_equal_the_brute_force_oracle(store):
    contigs, rows = store["contigs"], _rows(store)
    qual = {(ids.genome_uid(list(cs)), s): q for s in STAGES for cs, q in store["genomes"][s]}
    of_contig: dict[bytes, set] = {}
    of_genome: dict[bytes, set] = {}
    for c, g, s in rows:
        of_contig.setdefault(c, set()).add((g, s))
        of_genome.setdefault(g, set()).add((c, s))
    want_c = []
    for c in sorted(contigs):
        qs = {qual[gs] for gs in of_contig[c]}
        want_c.append({"contig_uid": c, "genome_uid": min(g for g, _ in of_contig[c]), "contig_len": contigs[c][1],
                       "is_circular": contigs[c][2], "tier": None,
                       "assembly_quality": qs.pop() if len(qs) == 1 else None, "tax_order": None, "is_species_rep": False, "n_genomes": len({g for g, _ in of_contig[c]})})
    want_g = []
    for g in sorted(of_genome):
        cs = {c for c, _ in of_genome[g]}
        qs = {qual[(g, s)] for _, s in of_genome[g]}
        want_g.append({"genome_uid": g, "n_contigs": len(cs), "genome_len": sum(contigs[c][1] for c in cs),
                       "kind": "prokaryote", "genetic_code": 11, "assembly_quality": qs.pop() if len(qs) == 1 else None,
                       "tier": None, "tax_order": None, "is_species_rep": False})
    root = store["root"]
    for name, want in (("contig_dim", want_c), ("genome_dim", want_g)):
        files = bio.part_files(root / table(name).location / f"release={REL}")
        got = [r for p in files for r in pq.read_table(p, partitioning=None).to_pylist()]
        assert got == want, name
    assert any(r["assembly_quality"] is None for r in want_c) and any(r["n_genomes"] > 2 for r in want_c)
    assert any(r["assembly_quality"] is None for r in want_g)
    shared = sum(r["n_genomes"] > 1 for r in want_c)
    assert store["res"] == {"contigs": len(want_c), "genomes": len(want_g),
                            "circular": sum(r["is_circular"] for r in want_c), "shared_contigs": shared}


@pytest.mark.covers("D-01")
def test_the_dims_do_not_depend_on_the_order_release_def_lists_the_stages(store):
    assert store["res_rev"] == store["res"]
    for name in ("contig_dim", "genome_dim"):
        loc = table(name).location
        a = bio.part_files(store["root"] / loc / f"release={REL}")
        b = bio.part_files(store["rev"] / loc / f"release={REL_REV}")
        assert a and [p.name for p in a] == [p.name for p in b], name
        assert [p.read_bytes() for p in a] == [p.read_bytes() for p in b], name


def test_genome_member_holds_every_row_of_the_release(store):
    got = store["con"].execute("SELECT contig_uid, genome_uid, stage FROM genome_member").fetchall()
    assert len(got) == len(set(got)) and set(got) == _rows(store)


def _keys(store, n: int, seed: int) -> list[bytes]:
    """`n` distinct contigs, half of them (at least one) in several genomes."""
    rng = np.random.default_rng(seed)
    by: dict[bytes, set] = {}
    for c, g, _ in _rows(store):
        by.setdefault(c, set()).add(g)
    multi = sorted(c for c, v in by.items() if len(v) > 1)
    rest = sorted(set(store["contigs"]) - set(multi))
    pick = [multi[i] for i in rng.choice(len(multi), size=min(max(1, n // 2), len(multi)), replace=False)]
    pick += [rest[i] for i in rng.choice(len(rest), size=n - len(pick), replace=False)]
    return pick


def _want(store, keys: list[bytes]) -> list[tuple]:
    ks = set(keys)
    return sorted({(c, g) for c, g, _ in _rows(store) if c in ks})


@pytest.mark.parametrize("form", ["parameter", "literal"])
def test_genomes_of_one_key_reads_its_row_group_in_each_stage_file(store, form, record_property):
    (key,) = _keys(store, 1, 3)
    sql = ("SELECT contig_uid, genome_uid FROM genomes_of(?)" if form == "parameter" else
           f"SELECT contig_uid, genome_uid FROM genomes_of([from_hex('{key.hex()}')])")
    prof = catalog.profile_query(store["con"], sql, [[key]] if form == "parameter" else None)
    assert sorted(prof["result"]) == _want(store, [key]) and len(prof["result"]) >= 2
    scanned, total = _gi_scanned(prof)
    footer, most = _footer_row_groups(_gi_files(store), [key])
    record_property("row_groups", f"{scanned} of {total}; footer bound {footer}")
    assert most <= MAX_RG_PER_KEY_AND_FILE
    assert 1 <= scanned <= footer <= MAX_RG_PER_KEY_AND_FILE * N_STAGES, (scanned, footer, total)


@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("form", ["join genome_member", "genomes_of the table's list"])
def test_genomes_of_keys_from_a_joined_table_read_only_their_row_groups(store, n, form, record_property):
    keys = _keys(store, n, 5)
    con = store["con"]
    con.register("_keys", pa.table({"contig_uid": pa.array(keys, pa.binary())}))
    try:
        sql = ("SELECT k.contig_uid, m.genome_uid FROM _keys k JOIN genome_member m ON m.contig_uid = k.contig_uid"
               if form == "join genome_member" else
               "SELECT contig_uid, genome_uid FROM genomes_of((SELECT list(contig_uid) FROM _keys))")
        prof = catalog.profile_query(con, sql)
    finally:
        con.unregister("_keys")
    assert sorted(set(prof["result"])) == _want(store, keys)
    scanned, total = _gi_scanned(prof)
    footer, most = _footer_row_groups(_gi_files(store), keys)
    record_property("row_groups", f"{scanned} of {total}; footer bound {footer}")
    assert most <= MAX_RG_PER_KEY_AND_FILE
    assert scanned <= footer <= MAX_RG_PER_KEY_AND_FILE * N_STAGES * n, (scanned, footer, total)
    assert footer < total  # the bound excludes row groups: the test can fail
