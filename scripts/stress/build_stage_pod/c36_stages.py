"""Card C-36: the six S2C phase-2 stages C1 to C6, each genome's stage, and each stage's lock and selection.

    python c36_stages.py assign REMAINDER_TSV NOT_SERVED_TXT PROPOSAL_TSV EXTRACT_DIR OUT_DIR
    python c36_stages.py sample STAGE_DIR N OUT_DIR

assign: the rule of the proposal (`planning/measurements/2026-11-13_s2c_phase2_fetch.md` section 5, DEC-34 (2)),
written here again from its text, not copied from the script that made the proposal: the rows of REMAINDER_TSV (the
ANI selection's remainder, card C-27) less the accessions of NOT_SERVED_TXT, ordered by `gtdb_taxonomy`, then
`species_rank` (as a number), then `accession`; n = round(total bp / 10 Gbp) with total bp the sum of `genome_size`;
each species (one `gtdb_taxonomy` string) goes whole to stage 1 + floor(n * bp before its first genome / total bp).
Stage k of the proposal (S2C-k) is named Ck (card C-36 (1)).

The genomes fetched are those of the two archives extracted under EXTRACT_DIR (`s2c_phase2/kept.tsv` and
`s2c_phase2_rest/kept.tsv`, with each batch's `lock.json` beside its files). Refused (exit 2, nothing written) when
the fetched accessions are not exactly the remainder less the ones not served, when an accession is in two locks, or
when a kept genome's bp differs from GTDB's `genome_size`. The per-stage table, in the proposal's columns, is written
to OUT_DIR/stage_table.tsv and compared byte for byte with PROPOSAL_TSV: a difference is reported (`equal_to_proposal`
false and the differing lines), never adjusted.

Writes, for each stage, OUT_DIR/C<k>/selection.tsv (the remainder's columns, the stage's rows in stage order) and
OUT_DIR/C<k>/lock.json (`genomes`: each accession's entry of its batch lock, with `file` the extracted .fna.gz and
`archived_as` its member name; `locks`: the sha256 of each batch lock it took entries from), and OUT_DIR/stages.tsv
(accession, stage, position in stage order, species, archive member). Prints a JSON report.

sample: the first N rows of STAGE_DIR/selection.tsv (stage order) and their lock entries, written to OUT_DIR as
selection.tsv and lock.json, for a rehearsal.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

GBP = 10_000_000_000
PREFIXES = ("s2c_phase2", "s2c_phase2_rest")  # the members' top directories in the two archives (sources.yaml)
COLS = ["stage", "genomes", "bp", "species", "genomes_in_species_over_20", "phyla", "first_species", "last_species"]


class Refused(Exception):
    pass


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return list(rd.fieldnames or []), list(rd)


def write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def fetched(extract: Path) -> tuple[dict[str, dict], dict[str, str]]:
    """accession -> {prefix, batch, kept row, lock entry} over both archives, and lock member -> its sha256."""
    out: dict[str, dict] = {}
    lock_sha: dict[str, str] = {}
    for prefix in PREFIXES:
        _, kept = read_tsv(extract / prefix / "kept.tsv")
        locks: dict[str, dict] = {}
        for lp in sorted((extract / prefix).glob("B*/lock.json")):
            data = lp.read_bytes()
            lock_sha[f"{prefix}/{lp.parent.name}/lock.json"] = hashlib.sha256(data).hexdigest()
            locks[lp.parent.name] = json.loads(data)["genomes"]
        for r in kept:
            acc, batch = r["accession"], r["batch"]
            if acc in out:
                raise Refused(f"{acc} is kept in {out[acc]['prefix']} and {prefix}")
            entry = locks.get(batch, {}).get(acc)
            if entry is None:
                raise Refused(f"{acc}: no entry in {prefix}/{batch}/lock.json")
            for k in ("gz_sha256", "sha256"):
                if entry[k] != r[k]:
                    raise Refused(f"{acc}: kept.tsv {k} {r[k]} is not its lock's {entry[k]}")
            out[acc] = {"prefix": prefix, "batch": batch, "kept": r, "lock": entry}
        in_locks = {a for g in locks.values() for a in g}
        if in_locks != {r["accession"] for r in kept}:
            raise Refused(f"{prefix}: the locks hold {len(in_locks)} genomes and kept.tsv {len(kept)}")
    return out, lock_sha


def assign(remainder: Path, not_served: Path, proposal: Path, extract: Path, out: Path) -> dict:
    cols, rows = read_tsv(remainder)
    ns = [ln.strip() for ln in not_served.read_text().splitlines() if ln.strip()]
    by_acc = {r["accession"]: r for r in rows}
    if len(by_acc) != len(rows):
        raise Refused(f"{remainder}: an accession is listed twice")
    if not set(ns) <= set(by_acc):
        raise Refused(f"{not_served}: {sorted(set(ns) - set(by_acc))} not in {remainder}")
    rows = [r for r in rows if r["accession"] not in set(ns)]
    got, lock_sha = fetched(extract)
    want = {r["accession"] for r in rows}
    if set(got) != want:
        raise Refused(f"fetched and remainder less not served differ: {len(set(got) - want)} fetched only, "
                      f"{len(want - set(got))} not fetched, e.g. {sorted(set(got) ^ want)[:5]}")
    bad_bp = [a for a in sorted(want) if int(got[a]["kept"]["bp"]) != int(by_acc[a]["genome_size"])
              or int(got[a]["lock"]["bp"]) != int(by_acc[a]["genome_size"])]
    if bad_bp:
        raise Refused(f"{len(bad_bp)} genomes whose bp is not GTDB's genome_size, e.g. {bad_bp[:5]}")

    rows.sort(key=lambda r: (r["gtdb_taxonomy"], int(r["species_rank"]), r["accession"]))
    total = sum(int(r["genome_size"]) for r in rows)
    n = max(1, round(total / GBP))
    stage_of: dict[str, int] = {}  # species -> stage
    before = 0
    members: dict[int, list[dict]] = {}
    for r in rows:
        sp = r["gtdb_taxonomy"]
        if sp not in stage_of:
            stage_of[sp] = 1 + (n * before) // total
        members.setdefault(stage_of[sp], []).append(r)
        before += int(r["genome_size"])
    # a species is contiguous in the order, so its stage never comes back after a later one
    order = [stage_of[r["gtdb_taxonomy"]] for r in rows]
    if order != sorted(order):
        raise Refused("stages are not in non-decreasing order along the taxonomy order")

    lines = ["\t".join(COLS)]
    rep_stages = {}
    stages_rows = []
    for k in sorted(members):
        ms = members[k]
        species = list(dict.fromkeys(r["gtdb_taxonomy"] for r in ms))
        line = [f"S2C-{k}", len(ms), sum(int(r["genome_size"]) for r in ms), len(species),
                sum(1 for r in ms if r["ani_pick"] != ""), len({s.split(";")[1] for s in species}),
                species[0].split(";")[-1], species[-1].split(";")[-1]]
        lines.append("\t".join(map(str, line)))
        name = f"C{k}"
        sel = "\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in ms)
        sel_sha = write_text(out / name / "selection.tsv", sel)
        genomes, used = {}, set()
        for r in ms:
            g = got[r["accession"]]
            member = f"{g['prefix']}/{g['batch']}/{r['accession']}.fna.gz"
            e = dict(g["lock"])
            e["fetched_file"] = e["file"]
            e["file"] = str(extract / member)
            e["archived_as"] = member
            genomes[r["accession"]] = e
            used.add(f"{g['prefix']}/{g['batch']}/lock.json")
        lock = {"card": "C-36", "stage": name, "genomes": dict(sorted(genomes.items())),
                "locks": {m: lock_sha[m] for m in sorted(used)}}
        lock_sha_out = write_text(out / name / "lock.json", json.dumps(lock, indent=1, sort_keys=True) + "\n")
        for i, r in enumerate(ms):
            g = got[r["accession"]]
            stages_rows.append((r["accession"], name, i, r["gtdb_taxonomy"].split(";")[-1],
                                f"{g['prefix']}/{g['batch']}/{r['accession']}.fna.gz"))
        per_species = [sum(1 for r in ms if r["gtdb_taxonomy"] == s) for s in species]
        bp_first = sum(int(r["genome_size"]) for i, r in enumerate(ms)
                       if i == 0 or ms[i - 1]["gtdb_taxonomy"] != r["gtdb_taxonomy"])
        rep_stages[name] = {
            "genomes": len(ms), "bp": line[2], "species": len(species), "genomes_in_species_over_20": line[4],
            "phyla": line[5], "from_archive": {p: sum(1 for r in ms if got[r["accession"]]["prefix"] == p)
                                               for p in PREFIXES},
            "fasta_bytes": sum(int(got[r["accession"]]["lock"]["bytes"]) for r in ms),
            "gz_bytes": sum(int(got[r["accession"]]["lock"]["gz_bytes"]) for r in ms),
            "records": sum(int(got[r["accession"]]["lock"]["records"]) for r in ms),
            "genomes_per_species": round(len(ms) / len(species), 3),
            "species_with_10_or_more": sum(1 for c in per_species if c >= 10),
            "genomes_in_species_with_10_or_more": sum(c for c in per_species if c >= 10),
            "bp_not_first_of_species": line[2] - bp_first,
            "share_bp_not_first_of_species": round((line[2] - bp_first) / line[2], 4),
            "largest_species": max(per_species),
            "translation_tables": {t: sum(1 for r in ms if r["translation_table"] == t)
                                   for t in sorted({r["translation_table"] for r in ms})},
            "gtdb_protein_count": sum(int(r["protein_count"]) for r in ms if r["protein_count"] not in ("", "None")),
            "gtdb_protein_count_missing": sum(1 for r in ms if r["protein_count"] in ("", "None")),
            "retrieved_utc": [min(e["retrieved_utc"] for e in genomes.values()),
                              max(e["retrieved_utc"] for e in genomes.values())],
            "selection_sha256": sel_sha, "lock_sha256": lock_sha_out}
    table = "\n".join(lines) + "\n"
    table_sha = write_text(out / "stage_table.tsv", table)
    ref = proposal.read_text()
    diff = [(a, b) for a, b in zip(ref.splitlines(), table.splitlines()) if a != b]
    write_text(out / "stages.tsv", "accession\tstage\tposition\tspecies\tmember\n"
               + "".join("\t".join(map(str, r)) + "\n" for r in stages_rows))
    return {"remainder_rows": len(by_acc), "not_served": len(ns), "genomes": len(rows), "bp": total, "n": n,
            "fetched": {p: sum(1 for g in got.values() if g["prefix"] == p) for p in PREFIXES},
            "locks": lock_sha, "stage_table_sha256": table_sha,
            "proposal_sha256": hashlib.sha256(ref.encode()).hexdigest(), "equal_to_proposal": ref == table,
            "differing_lines": diff, "stages": rep_stages}


def sample(stage_dir: Path, n: int, out: Path) -> dict:
    cols, rows = read_tsv(stage_dir / "selection.tsv")
    lock = json.loads((stage_dir / "lock.json").read_text())
    rows = rows[:n]
    sel = "\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in rows)
    sel_sha = write_text(out / "selection.tsv", sel)
    genomes = {r["accession"]: lock["genomes"][r["accession"]] for r in rows}
    doc = {"card": "C-36", "stage": f"{lock['stage']} sample of {len(rows)}", "genomes": dict(sorted(genomes.items())),
           "locks": lock["locks"]}
    lock_sha = write_text(out / "lock.json", json.dumps(doc, indent=1, sort_keys=True) + "\n")
    species = list(dict.fromkeys(r["gtdb_taxonomy"] for r in rows))
    return {"genomes": len(rows), "bp": sum(int(r["genome_size"]) for r in rows), "species": len(species),
            "first": rows[0]["accession"], "last": rows[-1]["accession"],
            "first_species": species[0].split(";")[-1], "last_species": species[-1].split(";")[-1],
            "last_species_cut": sum(1 for r in rows if r["gtdb_taxonomy"] == species[-1]),
            "selection_sha256": sel_sha, "lock_sha256": lock_sha}


def main(argv: list[str]) -> int:
    try:
        if argv[:1] == ["assign"] and len(argv) == 6:
            rep = assign(*map(Path, argv[1:]))
        elif argv[:1] == ["sample"] and len(argv) == 4:
            rep = sample(Path(argv[1]), int(argv[2]), Path(argv[3]))
        else:
            print(__doc__.split("\n\n")[1], file=sys.stderr)
            return 2
    except Refused as e:
        print(f"c36_stages.py: refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
