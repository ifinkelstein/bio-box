"""Card C-47 (Step 0): restore per-stage selection.tsv and lock.json for S2C phase-2 stages C1..C6
from the committed C-36 stages.tsv, the extracted phase-2 archives under /tmp/if2453/s2/phase2, and
GTDB R226 metadata on /stor (zero network downloads).

    python c47_restore.py restore STAGES_TSV GTDB_TAR EXTRACT_DIR OUT_DIR
    python c47_restore.py sample STAGE_DIR N OUT_DIR

`restore` reads `stages.tsv` (accession, stage, position, species, member), matches every accession
against the extracted archives' `kept.tsv` and `B*/lock.json` (verifying sha256s), looks up
`ncbi_translation_table` and `gtdb_taxonomy` from `gtdb_r226_metadata.tar` (`ar53_metadata_r226.tsv.gz`
and `bac120_metadata_r226.tsv.gz`), and writes `OUT_DIR/C<k>/selection.tsv` and `OUT_DIR/C<k>/lock.json`
in the exact format `make_inputs.py prok` expects.

`sample` takes the first N genomes of `STAGE_DIR/selection.tsv` (stage order) and writes a subset
`selection.tsv` and `lock.json` into `OUT_DIR` for the fixture/sample rehearsal.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import sys
import tarfile
from pathlib import Path

PREFIXES = ("s2c_phase2", "s2c_phase2_rest")


class Refused(Exception):
    pass


def _write_text(path: Path, text: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    tmp.replace(path)
    return hashlib.sha256(text.encode()).hexdigest()


def _read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return list(rd.fieldnames or []), list(rd)


def _fetched(extract: Path) -> tuple[dict[str, dict], dict[str, str]]:
    out: dict[str, dict] = {}
    lock_sha: dict[str, str] = {}
    for prefix in PREFIXES:
        _, kept = _read_tsv(extract / prefix / "kept.tsv")
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
                    raise Refused(f"{acc}: kept.tsv {k} {r[k]} != lock {entry[k]}")
            out[acc] = {"prefix": prefix, "batch": batch, "kept": r, "lock": entry}
        in_locks = {a for g in locks.values() for a in g}
        if in_locks != {r["accession"] for r in kept}:
            raise Refused(f"{prefix}: locks hold {len(in_locks)} genomes and kept.tsv {len(kept)}")
    return out, lock_sha


def _gtdb_meta(gtdb_tar: Path, want: set[str]) -> dict[str, tuple[str, str]]:
    """accession (stripped of GB_/RS_) -> (translation_table, gtdb_taxonomy)."""
    out: dict[str, tuple[str, str]] = {}
    with tarfile.open(gtdb_tar, "r:*") as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".tsv.gz"):
                continue
            raw = tf.extractfile(m)
            if raw is None:
                continue
            with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                txt = io.TextIOWrapper(gz, encoding="utf-8", errors="replace", newline="")
                rd = csv.DictReader(txt, delimiter="\t")
                for row in rd:
                    acc = row["accession"]
                    if acc.startswith(("GB_", "RS_")):
                        acc = acc[3:]
                    if acc not in want:
                        continue
                    tt = (row.get("ncbi_translation_table") or "").strip()
                    if tt.lower() in ("", "none", "na"):
                        tt = ""
                    tax = (row.get("gtdb_taxonomy") or "").strip()
                    out[acc] = (tt, tax)
    missing = want - set(out)
    if missing:
        raise Refused(f"{len(missing)} accessions missing from {gtdb_tar}: {sorted(missing)[:5]}")
    return out


def restore(stages_tsv: Path, gtdb_tar: Path, extract: Path, out: Path) -> dict:
    _, rows = _read_tsv(stages_tsv)
    got, lock_sha = _fetched(extract)
    want = {r["accession"] for r in rows}
    if set(got) != want:
        raise Refused(
            f"fetched ({len(got)}) != stages.tsv ({len(want)}): "
            f"extra={sorted(set(got) - want)[:5]}, missing={sorted(want - set(got))[:5]}"
        )
    meta = _gtdb_meta(gtdb_tar, want)
    by_stage: dict[str, list[dict]] = {}
    for r in rows:
        by_stage.setdefault(r["stage"], []).append(r)
    rep = {}
    for st in sorted(by_stage):
        ms = sorted(by_stage[st], key=lambda x: int(x["position"]))
        sel_lines = ["accession\ttranslation_table\tgtdb_taxonomy\tphylum\tspecies\tmember\n"]
        genomes: dict[str, dict] = {}
        used: set[str] = set()
        for r in ms:
            acc = r["accession"]
            g = got[acc]
            member = f"{g['prefix']}/{g['batch']}/{acc}.fna.gz"
            if member != r["member"]:
                raise Refused(f"{acc}: member {member} != stages.tsv {r['member']}")
            tt, tax = meta[acc]
            parts = tax.split(";")
            phylum = parts[1] if len(parts) > 1 else ""
            sel_lines.append(f"{acc}\t{tt}\t{tax}\t{phylum}\t{r['species']}\t{member}\n")
            e = dict(g["lock"])
            e["fetched_file"] = e["file"]
            e["file"] = str(extract / member)
            e["archived_as"] = member
            genomes[acc] = e
            used.add(f"{g['prefix']}/{g['batch']}/lock.json")
        sel_sha = _write_text(out / st / "selection.tsv", "".join(sel_lines))
        lock = {
            "card": "C-47",
            "stage": st,
            "genomes": dict(sorted(genomes.items())),
            "locks": {m: lock_sha[m] for m in sorted(used)},
        }
        lock_sha_out = _write_text(out / st / "lock.json", json.dumps(lock, indent=1, sort_keys=True) + "\n")
        rep[st] = {
            "genomes": len(ms),
            "selection_sha256": sel_sha,
            "lock_sha256": lock_sha_out,
        }
    return {"total_genomes": len(rows), "stages": rep}


def sample(stage_dir: Path, n: int, out: Path) -> dict:
    cols, rows = _read_tsv(stage_dir / "selection.tsv")
    lock = json.loads((stage_dir / "lock.json").read_text())
    rows = rows[:n]
    sel = "\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in rows)
    sel_sha = _write_text(out / "selection.tsv", sel)
    genomes = {r["accession"]: lock["genomes"][r["accession"]] for r in rows}
    doc = {
        "card": "C-47",
        "stage": f"{lock['stage']} sample of {len(rows)}",
        "genomes": dict(sorted(genomes.items())),
        "locks": lock["locks"],
    }
    lock_sha = _write_text(out / "lock.json", json.dumps(doc, indent=1, sort_keys=True) + "\n")
    return {"genomes": len(rows), "selection_sha256": sel_sha, "lock_sha256": lock_sha}


def main(argv: list[str]) -> int:
    try:
        if argv[:1] == ["restore"] and len(argv) == 5:
            rep = restore(*map(Path, argv[1:]))
        elif argv[:1] == ["sample"] and len(argv) == 4:
            rep = sample(Path(argv[1]), int(argv[2]), Path(argv[3]))
        else:
            print("usage: c47_restore.py restore STAGES_TSV GTDB_TAR EXTRACT_DIR OUT_DIR | sample STAGE_DIR N OUT_DIR",
                  file=sys.stderr)
            return 2
    except Refused as e:
        print(f"c47_restore.py: refused: {e}", file=sys.stderr)
        return 2
    print(json.dumps(rep, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
