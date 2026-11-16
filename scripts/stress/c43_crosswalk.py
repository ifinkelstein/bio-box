"""Card C-43: a crosswalk of metagenome assemblies by sequencing run.

    python -I c43_crosswalk.py spire-urls STUDIES_HTML OUT_TSV
    python -I c43_crosswalk.py labels --insdc-assemblies T --ncbi-summary T [--mgnify-samples DIR] OUT_TSV
    python -I c43_crosswalk.py build --spire-dir D --read-run T [T ...] --insdc-assemblies T --ncbi-summary T
        --mgnify-assemblies DIR [--mgnify-samples DIR] --habitat-map T [--insdc-include-types L] --out-dir D
    python -I c43_crosswalk.py lookup-urls UNITS_TSV OUT_TSV

Resources, in the order of preference: `spire` (one unit per SPIRE sample), `mgnify` (one per MGnify
assembly), `insdc` (one per INSDC assembly of a metagenome taxon, of the assembly types taken, that is
not MGnify's own, from ENA's assembly table with or without a WGS set, plus any NCBI Assembly row of
genbank/metagenomes that ENA's table lacks). A unit is
linked to runs by the runs its table names, else by its sample through the read-run table (primary
or secondary sample accession), else it is unlinked. Every input path is an argument; fetched files
are read as data only.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import json
import os
import re
import sys
import urllib.parse
from dataclasses import dataclass

SPIRE_STUDY_URL = "https://spire.embl.de/spire/api/study/{}?format=tsv"
_STUDY_ROW = re.compile(r'\{\s*study:\s*"((?:[^"\\]|\\.)*)",\s*samples:\s*(\d+)\s*\}')


def parse_studies(html: str) -> list[tuple[str, int]]:
    """The (study, sample count) rows of the `studiesData` list embedded in SPIRE's studies page."""
    start = html.index("var studiesData")
    end = html.index("];", start)
    return [(s, int(n)) for s, n in _STUDY_ROW.findall(html[start:end])]


def study_url(study: str) -> str:
    return SPIRE_STUDY_URL.format(urllib.parse.quote(study, safe=""))


def parse_micron(cell: str) -> list[str]:
    """SPIRE's microntology cell, a Python-style list of strings, read as data (never evaluated as code)."""
    cell = (cell or "").strip()
    if not cell:
        return []
    v = ast.literal_eval(cell)
    if not isinstance(v, list) or not all(isinstance(x, str) for x in v):
        raise ValueError(f"not a list of strings: {cell[:80]!r}")
    return v


# SPIRE's habitat top levels (the other microntology categories describe a sample, they are not habitats).
# A sample carries several tags (SPIRE's 'multiple tag' logic pulls 'aquatic' into sewage and mangroves), so
# one label per sample is the first top level present in this order, and level 2 is that top level's
# second path component as SPIRE writes it.
HABITAT_ORDER = ("host-associated", "anthropogenic", "terrestrial", "aquatic", "air")
UNLABELLED = "unlabelled"


def spire_habitat(terms: list[str]) -> tuple[str, str]:
    tops = {t.split(":", 1)[0] for t in terms}
    for top in HABITAT_ORDER:
        if top in tops:
            seconds = sorted(":".join(t.split(":")[:2]) for t in terms
                             if t.startswith(top + ":"))
            return top, (seconds[0] if seconds else top)
    return UNLABELLED, UNLABELLED


def mapped_label(value: str | None, table: dict[str, str]) -> str:
    """The label of `value` by the longest key that equals it or is its prefix at a ':' boundary."""
    if not value:
        return UNLABELLED
    best = None
    for k, lab in table.items():
        if value == k or value.startswith(k + ":"):
            if best is None or len(k) > len(best[0]):
                best = (k, lab)
    return best[1] if best else UNLABELLED


class RunIndex:
    """Runs of the read-run table by sample (primary or secondary accession), and their read bases."""

    def __init__(self, rows):
        self.by_sample: dict[str, set[str]] = {}
        self.primary: dict[str, str] = {}
        self._bases: dict[str, int | None] = {}
        self._name: dict[str, str] = {}
        for row in rows:
            run, sample, secondary, bases = row[:4]
            self._bases[run] = bases
            if len(row) > 4 and row[4]:
                self._name[run] = row[4]
            if sample:
                self.by_sample.setdefault(sample, set()).add(run)
                if secondary:
                    for s in secondary.split(";"):
                        if s:
                            self.primary[s] = sample

    def norm(self, sample: str) -> str:
        return self.primary.get(sample, sample)

    def bases(self, run: str) -> int | None:
        return self._bases.get(run)

    def has(self, run: str) -> bool:
        return run in self._bases

    def organism(self, runs) -> str:
        """The INSDC organism name (e.g. 'human gut metagenome') of the first of `runs` that has one."""
        return next((self._name[r] for r in sorted(runs) if r in self._name), "")

    def link(self, samples, direct_runs) -> tuple[set[str], set[str], str]:
        norm = {self.norm(s) for s in samples if s}
        if direct_runs:
            return set(direct_runs), norm, "run"
        runs = set()
        for s in norm:
            runs |= self.by_sample.get(s, set())
        return runs, norm, ("sample" if runs else "none")


_RUN = re.compile(r"^[EDS]RR\d+$")
_INSDC_SAMPLE = re.compile(r"^(SAM(EA|N|D)\d+|[EDS]RS\d+)$")


def spire_link(sample_id: str, idx: RunIndex) -> tuple[set[str], set[str], str]:
    """A SPIRE sample id that is itself a run accession is that run; otherwise a sample accession."""
    if _RUN.match(sample_id):
        return {sample_id}, set(), "run"
    return idx.link([sample_id], [])


def insdc_samples(ids) -> list[str]:
    return sorted({s for s in ids if s and _INSDC_SAMPLE.match(s)})


def wgs_included(assembly_type: str, include: set[str]) -> bool:
    return (assembly_type or "").strip().lower() in include


@dataclass(frozen=True)
class Unit:
    resource: str
    uid: str
    runs: frozenset
    samples: frozenset
    linked_by: str
    asm_bases: int | None
    habitat: str
    habitat2: str = ""
    note: str = ""
    gz_bytes: int | None = None   # size of the provider's compressed contig files, where a table states it


def _ordered(units, order):
    rank = {r: i for i, r in enumerate(order)}
    return sorted(units, key=lambda u: (rank[u.resource], u.uid))


def choose(units, order) -> dict[str, str]:
    """One assembly per run: resources in `order`, units of a resource by uid; a unit is 'take' when none of its
    runs is taken yet and none of its samples was assembled by a taken unit, 'covered' when all of its runs are
    taken or its sample was, 'partial' when some but not all of its runs are taken (not taken: it would repeat
    reads), 'take_unlinked' when it has no run to check and its sample is new."""
    runs_taken: set[str] = set()
    samples_taken: set[str] = set()
    out = {}
    for u in _ordered(units, order):
        hit = u.runs & runs_taken
        if u.runs and hit == u.runs:
            out[u.uid] = "covered"
        elif hit:
            out[u.uid] = "partial"
        elif u.samples & samples_taken:
            out[u.uid] = "covered"
        else:
            out[u.uid] = "take" if u.runs else "take_unlinked"
            runs_taken |= u.runs
            samples_taken |= u.samples
    return out


def run_resources(units, order) -> dict[str, tuple[tuple[str, ...], str]]:
    """For each run: the resources holding it, in `order`, and the habitat label of the first of them."""
    held: dict[str, list] = {}
    for u in _ordered(units, order):
        for r in u.runs:
            held.setdefault(r, []).append(u)
    out = {}
    for r, us in held.items():
        res = tuple(dict.fromkeys(u.resource for u in us))
        out[r] = (res, us[0].habitat)
    return out


def combo_table(per_run, idx: RunIndex) -> dict:
    t: dict = {}
    for run, (res, hab) in per_run.items():
        cell = t.setdefault(res, {}).setdefault(hab, {"runs": 0, "read_bases": 0, "runs_without_bases": 0})
        cell["runs"] += 1
        b = idx.bases(run)
        if b is None:
            cell["runs_without_bases"] += 1
        else:
            cell["read_bases"] += b
    return t


def _spire_urls(a) -> int:
    with open(a.studies_html, encoding="utf-8") as fh:
        rows = parse_studies(fh.read())
    names = set()
    with open(a.out_tsv, "w") as out:
        for study, _ in rows:
            name = urllib.parse.quote(study, safe="") + ".tsv"
            if name in names:
                raise SystemExit(f"duplicate study file name {name}")
            names.add(name)
            out.write(f"{study_url(study)}\t{name}\n")
    print(f"{len(rows)} studies, {sum(n for _, n in rows)} samples listed")
    return 0


# ---------------------------------------------------------------------------------------------- readers

ORDER = ("spire", "mgnify", "insdc")


def _tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        head = next(r)
        while head and head[0].startswith("##"):
            head = next(r)
        head = [h.lstrip("#") for h in head]
        for row in r:
            if row and not row[0].startswith("#"):
                yield dict(zip(head, row))


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _unver(acc: str) -> str:
    return acc.split(".", 1)[0] if acc else ""


def read_runs(paths) -> RunIndex:
    def rows():
        for path in paths:
            for r in _tsv(path):
                if r.get("library_strategy", "").upper() == "AMPLICON":   # reads of one gene assemble no contigs
                    continue
                yield (r["run_accession"], r["sample_accession"], r["secondary_sample_accession"],
                       _int(r["base_count"]), r.get("scientific_name", ""))
    return RunIndex(rows())


def read_spire(study_dir) -> list[tuple[str, str, list[str]]]:
    """(study, sample, microntology terms) for every sample of every study table; probe files are skipped."""
    out = []
    for name in sorted(os.listdir(study_dir)):
        if not name.endswith(".tsv") or name.startswith("probe_"):
            continue
        study = urllib.parse.unquote(name[:-4])
        for r in _tsv(os.path.join(study_dir, name)):
            out.append((study, r["sample_id"], parse_micron(r.get("microntology", ""))))
    return out


def _pages(d):
    for name in sorted(os.listdir(d)):
        if name.startswith("page_") and name.endswith(".json"):
            with open(os.path.join(d, name)) as fh:
                yield from json.load(fh)["data"]


def _ids(rel):
    data = (rel or {}).get("data")
    if data is None:
        return []
    if isinstance(data, dict):
        return [data["id"]]
    return [x["id"] for x in data]


def read_habitat_map(path) -> dict[str, dict[str, str]]:
    m: dict[str, dict[str, str]] = {"mgnify_lineage": {}, "organism_name": {}}
    for r in _tsv(path):
        m[r["kind"]][r["key"]] = r["label"]
    return m


# ---------------------------------------------------------------------------------------------- build

def build_units(a, idx: RunIndex, hmap):
    units: list[Unit] = []
    # SPIRE: one unit per sample, linked to runs by the sample accession
    seen = set()
    for study, sample, terms in read_spire(a.spire_dir):
        if sample in seen:
            continue
        seen.add(sample)
        runs, samples, how = spire_link(sample, idx)
        h1, h2 = spire_habitat(terms)
        units.append(Unit("spire", sample, frozenset(runs), frozenset(samples), how, None, h1, h2,
                          note=f"study={study}"))
    # MGnify: assemblies with runs and samples inline; bases from the INSDC record of the same assembly
    wgs_rows = list(_tsv(a.insdc_assemblies))
    wgs_by_asm = {_unver(r["assembly_accession"]): r for r in wgs_rows if r["assembly_accession"]}
    wgs_by_set = {_wgs_key(r["wgs_set"]): r for r in wgs_rows if _wgs_key(r.get("wgs_set", ""))}
    biome = {}
    if a.mgnify_samples:
        for s in _pages(a.mgnify_samples):
            b = _ids(s["relationships"].get("biome"))
            biome[s["id"]] = b[0] if b else None
            bs = (s.get("attributes") or {}).get("biosample")
            if bs:
                biome.setdefault(bs, b[0] if b else None)
    emg = {}
    if a.mgnify_ena_analyses:
        for r in _tsv(a.mgnify_ena_analyses):
            if r["assembly_type"].strip().lower() != "primary metagenome":
                continue
            sizes = [_int(x) for x in r["submitted_bytes"].split(";") if x]
            emg[r["analysis_accession"]] = {
                "gz": sum(x for x in sizes if x is not None) if sizes else None,
                "runs": [x for x in r["run_accession"].split(";") if x],
                "samples": [x for x in r["sample_accession"].split(";") if x]}
    def mg_label(lineage, runs):
        """MGnify's biome of the sample; where MGnify gives none, the INSDC organism of a run, both mapped
        through the stated table."""
        if lineage:
            return mapped_label(lineage, hmap["mgnify_lineage"]), lineage
        org = idx.organism(runs)
        if org:
            return mapped_label(org, hmap["organism_name"]), "insdc:" + org
        return UNLABELLED, ""

    mg_asm, mg_wgs = set(), set()
    for d in _pages(a.mgnify_assemblies):
        at = d["attributes"]
        if at.get("is-private"):
            continue
        legacy, wgsacc = at.get("legacy-accession") or "", at.get("wgs-accession") or ""
        w = wgs_by_asm.get(_unver(legacy)) if legacy.startswith("GC") else None
        if w is None and _wgs_key(wgsacc):
            w = wgs_by_set.get(_wgs_key(wgsacc))
        if w is not None:
            mg_wgs.add(_insdc_uid(w))
        if legacy.startswith("GC"):
            mg_asm.add(_unver(legacy))
        samples_raw = _ids(d["relationships"].get("samples"))
        runs, samples, how = idx.link(samples_raw, _ids(d["relationships"].get("runs")))
        lineage = next((biome[s] for s in samples_raw if biome.get(s)), None)
        e = emg.pop(d["id"], None)
        units.append(Unit("mgnify", d["id"], frozenset(runs), frozenset(samples), how,
                          _int(w["base_count"]) if w else None, *mg_label(lineage, runs),
                          note=f"insdc={_insdc_uid(w) if w else ''}", gz_bytes=e["gz"] if e else None))
    # MGnify's primary assemblies that ENA holds and MGnify's API does not list
    for erz, e in sorted(emg.items()):
        runs, samples, how = idx.link(e["samples"], e["runs"])
        lineage = next((biome[s] for s in list(e["samples"]) + sorted(samples) if biome.get(s)), None)
        units.append(Unit("mgnify", erz, frozenset(runs), frozenset(samples), how, None, *mg_label(lineage, runs),
                          note="source=ena_analysis_only", gz_bytes=e["gz"]))
    # INSDC metagenome assemblies that are not MGnify's own; the NCBI Assembly table marks its subset
    ncbi = {_unver(r["assembly_accession"]): r for r in _tsv(a.ncbi_summary)}
    include = {t.strip().lower() for t in a.insdc_include_types.split(",")}
    in_wgs = set()
    excluded = {}
    for r in wgs_rows:
        asm = _unver(r["assembly_accession"])
        in_wgs.add(asm)
        if _insdc_uid(r) in mg_wgs or (asm and asm in mg_asm):
            k = ("mgnify's own", r.get("assembly_type", ""))
        elif not wgs_included(r.get("assembly_type", ""), include):
            k = ("type not taken", r.get("assembly_type", ""))
        else:
            k = None
        if k is not None:
            e = excluded.setdefault(k, [0, 0, 0, 0])
            e[0] += 1
            e[1] += _int(r["base_count"]) or 0
            e[2] += asm in ncbi
            e[3] += bool(r["run_accession"])
            continue
        direct = [x for x in r["run_accession"].split(";") if x]
        runs, samples, how = idx.link([r["sample_accession"], r["secondary_sample_accession"]], direct)
        units.append(Unit("insdc", _insdc_uid(r), frozenset(runs), frozenset(samples), how, _int(r["base_count"]),
                          mapped_label(r["scientific_name"], hmap["organism_name"]), r["scientific_name"],
                          note=f"ncbi_assembly={'yes' if asm in ncbi else 'no'};wgs_set={r.get('wgs_set', '')};"
                               f"center={r.get('center_name', '')}"))
    # NCBI Assembly rows that ENA's table lacks (kept, so nothing is dropped; their type is not stated)
    for asm, r in sorted(ncbi.items()):
        if asm in in_wgs or asm in mg_asm:
            continue
        runs, samples, how = idx.link([r["biosample"]], [])
        units.append(Unit("insdc", asm, frozenset(runs), frozenset(samples), how, _int(r["genome_size"]),
                          mapped_label(r["organism_name"], hmap["organism_name"]), r["organism_name"],
                          note="ncbi_assembly=yes;not_in_ena_table"))
    return units, excluded


_WGS = re.compile(r"^([A-Z]{4,6}\d{2})")


def _wgs_key(v: str) -> str:
    m = _WGS.match(v or "")
    return m.group(1) if m else ""


def _insdc_uid(r) -> str:
    return _unver(r.get("assembly_accession", "")) or _wgs_key(r.get("wgs_set", ""))


def _write_tsv(path, header, rows):
    with open(path + ".part", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    os.replace(path + ".part", path)


def _build(a) -> int:
    os.makedirs(a.out_dir, exist_ok=True)
    idx = read_runs(a.read_run)
    hmap = read_habitat_map(a.habitat_map)
    units, excluded = build_units(a, idx, hmap)
    _write_tsv(os.path.join(a.out_dir, "insdc_not_units.tsv"),
               ["why", "assembly_type", "assemblies", "bases", "in_ncbi_assembly", "with_run_accession"],
               [[w, t, *v] for (w, t), v in sorted(excluded.items())])
    dec = choose(units, ORDER)
    per_run = run_resources(units, ORDER)

    def rb(u):
        known = [idx.bases(r) for r in u.runs]
        return sum(b for b in known if b is not None), sum(1 for b in known if b is None)

    # per unit (large: stays on the pod)
    rows = []
    for u in _ordered(units, ORDER):
        b, nb = rb(u)
        rows.append([u.resource, u.uid, dec[u.uid], u.linked_by, len(u.runs), ";".join(sorted(u.runs)),
                     ";".join(sorted(u.samples)), "" if u.asm_bases is None else u.asm_bases, b, nb,
                     u.habitat, u.habitat2, u.note, "" if u.gz_bytes is None else u.gz_bytes])
    _write_tsv(os.path.join(a.out_dir, "units.tsv"),
               ["resource", "uid", "decision", "linked_by", "n_runs", "runs", "samples", "asm_bases",
                "read_bases", "runs_without_bases", "habitat", "habitat_detail", "note", "gz_bytes"], rows)

    # resources: units, how linked, runs, bases
    summ = {}
    for u in units:
        s = summ.setdefault(u.resource, {"units": 0, "run": 0, "sample": 0, "none": 0, "asm_bases": 0,
                                         "units_with_asm_bases": 0, "runs": set(), "gz": 0, "with_gz": 0})
        s["units"] += 1
        s[u.linked_by] += 1
        s["runs"] |= u.runs
        if u.asm_bases is not None:
            s["asm_bases"] += u.asm_bases
            s["units_with_asm_bases"] += 1
        if u.gz_bytes is not None:
            s["gz"] += u.gz_bytes
            s["with_gz"] += 1
    out = []
    for res in ORDER:
        s = summ.get(res)
        if not s:
            continue
        known = [idx.bases(r) for r in s["runs"]]
        out.append([res, s["units"], s["run"], s["sample"], s["none"], len(s["runs"]),
                    sum(b for b in known if b is not None), sum(1 for b in known if b is None),
                    s["units_with_asm_bases"], s["asm_bases"], s["with_gz"], s["gz"]])
    _write_tsv(os.path.join(a.out_dir, "resources.tsv"),
               ["resource", "units", "linked_by_run", "linked_by_sample", "unlinked", "distinct_runs",
                "read_bases", "runs_without_bases", "units_with_asm_bases", "asm_bases", "units_with_gz_bytes",
                "gz_bytes"], out)

    # runs by the combination of resources holding them, and habitat
    t = combo_table(per_run, idx)
    out = []
    for res in sorted(t, key=lambda c: ([ORDER.index(x) for x in c])):
        for hab in sorted(t[res]):
            c = t[res][hab]
            out.append(["+".join(res), hab, c["runs"], c["read_bases"], c["runs_without_bases"]])
    _write_tsv(os.path.join(a.out_dir, "runs_by_combination_and_habitat.tsv"),
               ["resources", "habitat", "runs", "read_bases", "runs_without_bases"], out)

    # what to take from each resource, by habitat and decision
    agg = {}
    for u in units:
        k = (u.resource, u.habitat, dec[u.uid])
        g = agg.setdefault(k, [0, 0, 0, 0, 0, 0])
        if u.gz_bytes is not None:
            g[5] += u.gz_bytes
        b, _ = rb(u)
        g[0] += 1
        g[1] += len(u.runs)
        g[2] += b
        if u.asm_bases is not None:
            g[3] += 1
            g[4] += u.asm_bases
    out = [[r, h, d, *v] for (r, h, d), v in sorted(agg.items(), key=lambda kv: (ORDER.index(kv[0][0]), kv[0][1:]))]
    _write_tsv(os.path.join(a.out_dir, "decisions_by_habitat.tsv"),
               ["resource", "habitat", "decision", "units", "runs", "read_bases", "units_with_asm_bases",
                "asm_bases", "gz_bytes"], out)

    # SPIRE level 2 (for card C-42)
    lvl2 = {}
    for u in units:
        if u.resource == "spire":
            g = lvl2.setdefault(u.habitat2, [0, 0, 0])
            b, _ = rb(u)
            g[0] += 1
            g[1] += len(u.runs)
            g[2] += b
    _write_tsv(os.path.join(a.out_dir, "spire_habitat_level2.tsv"),
               ["habitat_detail", "samples", "runs", "read_bases"], [[k, *v] for k, v in sorted(lvl2.items())])

    # labels the map left unlabelled, with their units (to show what the map misses)
    miss = {}
    for u in units:
        if u.resource != "spire" and u.habitat == UNLABELLED:
            k = (u.resource, u.habitat2)
            miss[k] = miss.get(k, 0) + 1
    _write_tsv(os.path.join(a.out_dir, "unlabelled_sources.tsv"), ["resource", "label_given", "units"],
               [[r, l, n] for (r, l), n in sorted(miss.items(), key=lambda kv: -kv[1])])

    # take lists (one per resource): the units to import, in the order of preference
    for res in ORDER:
        _write_tsv(os.path.join(a.out_dir, f"take_{res}.tsv"),
                   ["uid", "decision", "runs", "asm_bases", "habitat", "habitat_detail"],
                   [[u.uid, dec[u.uid], ";".join(sorted(u.runs)), "" if u.asm_bases is None else u.asm_bases,
                     u.habitat, u.habitat2]
                    for u in _ordered(units, ORDER) if u.resource == res and dec[u.uid].startswith("take")])
    print(json.dumps({"units": len(units), "runs": len(per_run)}))
    return 0


ENA_BY_SAMPLE = ("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query={}&fields=run_accession,"
                 "sample_accession,secondary_sample_accession,study_accession,base_count,library_source,"
                 "library_strategy,tax_id,scientific_name&limit=0&format=tsv")


def lookup_query(ids) -> str:
    """ENA query naming each sample by its field: BioSample accessions are primary, ERS/SRS/DRS secondary."""
    return " OR ".join(f'{"sample_accession" if s.startswith("SAM") else "secondary_sample_accession"}="{s}"'
                       for s in ids)


def lookup_url(ids) -> str:
    return ENA_BY_SAMPLE.format(urllib.parse.quote(lookup_query(ids), safe=""))


def _lookup_urls(a) -> int:
    """URLs of ENA's runs for the INSDC samples of unlinked units, in batches (the ena_portal_tables entry)."""
    ids = []
    for r in _tsv(a.units_tsv):
        if r["linked_by"] == "none":
            ids.extend(r["samples"].split(";"))
    ids = insdc_samples(ids)
    with open(a.out_tsv, "w") as out:
        for i in range(0, len(ids), a.batch):
            out.write(lookup_url(ids[i:i + a.batch]) + f"\tbatch_{i // a.batch:04d}.tsv\n")
    print(f"{len(ids)} samples in {(len(ids) + a.batch - 1) // a.batch} batches")
    return 0


def _labels(a) -> int:
    """Distinct labels the other resources give, with unit counts, to write the habitat map from."""
    cnt = {}
    for r in _tsv(a.insdc_assemblies):
        k = ("organism_name", r["scientific_name"])
        cnt[k] = cnt.get(k, 0) + 1
    for r in _tsv(a.ncbi_summary):
        k = ("organism_name", r["organism_name"])
        cnt[k] = cnt.get(k, 0) + 1
    if a.mgnify_samples:
        for s in _pages(a.mgnify_samples):
            b = _ids(s["relationships"].get("biome"))
            k = ("mgnify_lineage", b[0] if b else "")
            cnt[k] = cnt.get(k, 0) + 1
    _write_tsv(a.out_tsv, ["kind", "label_given", "units"],
               [[k, l, n] for (k, l), n in sorted(cnt.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1]))])
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("spire-urls")
    p.add_argument("studies_html")
    p.add_argument("out_tsv")
    p.set_defaults(fn=_spire_urls)
    p = sub.add_parser("labels")
    p.add_argument("--insdc-assemblies", required=True)
    p.add_argument("--ncbi-summary", required=True)
    p.add_argument("--mgnify-samples", default=None)
    p.add_argument("out_tsv")
    p.set_defaults(fn=_labels)
    p = sub.add_parser("build")
    for opt in ("--spire-dir", "--insdc-assemblies", "--ncbi-summary", "--mgnify-assemblies", "--habitat-map", "--out-dir"):
        p.add_argument(opt, required=True)
    p.add_argument("--read-run", required=True, nargs="+")
    p.add_argument("--mgnify-samples", default=None)
    p.add_argument("--mgnify-ena-analyses", default=None)
    p.add_argument("--insdc-include-types", default="primary metagenome",
                   help="comma-separated ENA assembly types taken as WGS units; an empty item takes a blank type")
    p.set_defaults(fn=_build)
    p = sub.add_parser("lookup-urls")
    p.add_argument("units_tsv")
    p.add_argument("out_tsv")
    p.add_argument("--batch", type=int, default=100)
    p.set_defaults(fn=_lookup_urls)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
