"""synthworld null worlds (planning/tests/synthworld_spec.md §4.1; card F-02 session 9; MR-18, B-NULL). A null world
is a world of its own: the same stage code as the signal world writes every stage directory, the null's change is
applied at the stage where the planted association is made, and S11 (`emit`) writes it with its own world id, truth,
digest and ``manifest.json``, so it is emitted, built and scored exactly as a signal world.

**N0-full** (`build`, ``null`` "N0-full"; §4.1 primary null, B-NULL). The signal preset and seed with zero planted
groups (`null_preset`): ``modules.count`` 0 and ``profile_only.count`` 0 (S2: the groups are not declared to ZOMBI2,
so their families never exist), ``decoys.partner_modules`` empty (S2: a ``Recipients`` decoy needs a module partner),
``fusions.nogd`` 0 (S7: NOGD displaces a module member). Everything else runs unchanged.

- *Kept exactly* (same seed): every other preset value (rates, extents, plasmids, vehicles, hubs, phages, prophages,
  CRISPR, fusions and splits, call errors, assembly, contamination, duplicates, labels); the species tree with its
  dead lineages, the clonal burst and the taxonomy (ZOMBI2 draws the tree from its own seed before any gene); the
  per-lineage codon usage, clock and intergenic composition (S4, S5 lineage tables come from the tree and the seed);
  the phage set (S6 builds it from S1-S5 sections that the null does not touch, `mobile.phage_preset`); everything
  drawn per genome from the tree: the code-4 clade and the code-25 genome, the genome of each call error, each
  genome's assembly class, N50 target and completeness target, the contaminated and duplicated genomes and their
  donors, accessions, BioProjects, labels and biomes; the background tight-junction rule (§2.4: a hash of the seed
  and the two families, so a pair of root families gets the same class in both worlds).
- *Kept in distribution*: gene content and order. ZOMBI2's gene-level draws differ once the groups change (§4.1
  "verify": the R2 test measures it), so genome sizes, family sizes, phyletic patterns, strands, operons and
  composition are new draws of the same process; the R2 test states each comparison.
- *Destroyed*: every planted association: no planted or ``profile_only`` group, no decoy, no NOGD; ``truth_module``
  has no planted or ``profile_only`` rows.

**N0-strict** (`build`, "N0-strict"). N0-full and also no vehicles (``vehicles.cargo`` empty, S2) and no CRISPR-Cas
loci (``crispr.arrays`` 0, S6: an array is planted with its Cas operon as one locus, so both go; prophages are placed
first from the same stream and are kept). Plasmids stay (§5.3). Every pair of families is null.

**N1** (`build_n1`; §4.1, MR-18). A profile-preserving relabel of a *signal* world. For each member A of each
planted (adjacent) module, in module and member order, a replacement A' is drawn uniformly (stream
``N1/relabel/<A>`` of the null seed, the world seed by default) among the free families of backbone kind
``nulls.n1.kinds`` whose phyletic profile (the set of host genomes holding a copy in ``assembly/genes.tsv``, the S8
genomes before fragmentation; a fused gene counts for both families) differs from A's in at most
floor(``max_hamming_fraction`` x genomes) genomes and whose transfer count (``backbone/events.tsv``) is within
``transfer_tolerance`` x A's; a member without such a family stays (``unmatched``). Nothing in the genomes moves: S1-S10
are the signal world's files (hard links where the file system allows), so the emitted FASTA is byte for byte the
signal world's. The swap is applied to the truth at S11: module m's member A is listed as A'. By renaming A and A'
in every table this is the spec's "swap the labels of A and A' in every genome"; the two are the same world up to
names, and keeping the gene labels keeps ``truth_gene.family_id`` the family whose sequence a gene carries.
``truth_gene.module_id`` and ``truth_fusion.in_module`` follow the new membership; the NOGD rows of the modules go
(they describe the planted unit, which no module lists any more), ``realized_loci`` of planted rows is empty.

- *Kept*: every genome base; per label, the phyletic profile within the Hamming bound, the transfer count within the
  tolerance, so the co-occurrence of two module labels changes by at most the sum of their two Hamming distances.
- *Destroyed*: the adjacency coupling between module labels (the planted units carry labels no module lists).

**N2** (`n2_seeds`). Calibration over ``nulls.n2.seeds`` N0 seeds, derived from the master seed by stream
``N2/seeds`` (one u64 per seed); each is an ordinary N0 world built with `build`.

**What the truth of a null world says** (S11 writes ``truth/truth_null.tsv`` and ``truth/truth_null_relabel.tsv`` in a
null world only, so a signal world's emit stage is unchanged): ``truth_null`` has one row per removed item (``action``
removed: the planted modules, ``profile_only`` groups, decoys and NOGD events of the signal preset, and in N0-strict
the vehicles and CRISPR-Cas loci) and one row per rule of null evaluation (``action`` excluded, null_trap, null,
evaluated or not_evaluated, ``family_ids`` the set whose within-set pairs the rule covers, ``*`` for every pair). The
decisions of session 9, for the PI to confirm:

- *Decoys* are removed from N0 (§4.1), with their partners.
- *NOGD* is removed from N0 (§4.1) and from N1's truth.
- *Cas loci*: §4.1 says N0's ``truth_module`` has no ``cas_operon`` rows, keeps CRISPR arrays, and names N0-strict as
  the variant without the Cas group. Read together: N0-full keeps every CRISPR-Cas locus as S6 makes it (an inherited
  confounder: one locus per clade, never transferred), keeps its ``cas_operon`` rows, and excludes the pairs within a
  locus's Cas families from null evaluation, as it excludes chassis pairs; N0-strict drops the loci (operon, leader
  and array) at S6 without touching the backbone, so every pair is null there.
- *Vehicle pairs*: the chassis pair of each vehicle is excluded (§3.3, §4.1); pairs involving cargo are null traps
  (``truth_traps`` mge_cargo, must_not_call: cargo always shares chassis flanks, so a count by vehicle-independent
  locus collapses them).
- *Plasmid pairs* are null in both N0 variants: a plasmid's genes share their flanks wherever it goes, so a count by
  locus collapses them as it does cargo, and §5.3 keeps plasmids in N0-strict, whose every pair is null.
- *N1* evaluates only the pairs of each module's labels that have at least one relabeled member; the planted units,
  now under unlisted labels, are excluded, and every other pair is not evaluated.
"""

from __future__ import annotations

import collections
import copy
import csv
import json
import math
import os
import shutil
from pathlib import Path
from typing import NamedTuple

import numpy as np

from biobox.synthworld import assembly, backbone, emit, fusions, intergenic, mobile, seeds, sequences, sources

NULLS = ("N0-full", "N0-strict", "N1")
#: the preset leaves each N0 variant sets (§4.1), as (section, key, value); the signal value is kept in the marker
N0_CHANGES = {
    "N0-full": (("modules", "count", 0), ("profile_only", "count", 0), ("decoys", "partner_modules", []),
                ("fusions", "nogd", 0)),
}
N0_CHANGES["N0-strict"] = N0_CHANGES["N0-full"] + (("vehicles", "cargo", []), ("crispr", "arrays", 0))
#: the stage directories S1-S10 write (the backbone also writes ``zombi2/``); N1 shares them with its signal world
SHARED = ("backbone", "zombi2", "sequences", "replicons", "mobile", "fusions", "assembly", "sources")


class RelabelRow(NamedTuple):
    module_id: str
    order_index: int
    family_id: str
    replacement_family_id: str
    status: str                  # relabeled | unmatched
    hamming: int                 # -1 when unmatched
    n_genomes: int
    replacement_n_genomes: int   # -1 when unmatched
    transfers: int
    replacement_transfers: int   # -1 when unmatched


class ScopeRow(NamedTuple):
    kind: str
    item_id: str
    family_ids: str
    action: str
    rule: str


def check_config(cfg: dict) -> None:
    n1, n2 = cfg["n1"], cfg["n2"]
    if not 0.0 <= float(n1["max_hamming_fraction"]) <= 1.0:
        raise ValueError("nulls.n1.max_hamming_fraction is a fraction of genomes")
    if float(n1["transfer_tolerance"]) < 0.0:
        raise ValueError("nulls.n1.transfer_tolerance is >= 0")
    if not n1["kinds"] or set(n1["kinds"]) - {"root", "originated"}:
        raise ValueError("nulls.n1.kinds are non-module backbone kinds: root, originated")
    if int(n2["seeds"]) < 20:
        raise ValueError("nulls.n2.seeds: §4.1 asks for at least 20 N0 seeds")


def null_preset(params: dict, null: str) -> dict:
    """The preset of an N0 world: the signal preset with the planted groups removed (module docstring). The marker
    ``null`` names the null world and keeps the signal values it changed."""
    if null == "N1":
        raise ValueError("N1 changes no preset: it relabels a signal world's truth (build_n1)")
    if null not in N0_CHANGES:
        raise ValueError(f"unknown null world {null!r}; one of {NULLS}")
    if "null" in params:
        raise ValueError(f"the preset is already a null world's ({params['null']['world']})")
    p = copy.deepcopy(params)
    removed = {}
    for section, key, value in N0_CHANGES[null]:
        removed[f"{section}.{key}"] = copy.deepcopy(p[section][key])
        p[section][key] = copy.deepcopy(value)
    p["null"] = {"world": null, "removed": removed}
    return p


def n2_seeds(seed: int, n: int) -> list[int]:
    """The N2 calibration seeds: n u64 words from stream ``N2/seeds`` (a prefix is stable as n grows)."""
    w = seeds.seed_sequence(seed, "N2/seeds").generate_state(2 * int(n), dtype=np.uint32)
    return [(int(w[2 * i]) << 32) | int(w[2 * i + 1]) for i in range(int(n))]


# --- building -------------------------------------------------------------------------------------------------

def build_stages(params: dict, seed: int, out, *, pfam_hmm=None, rfam_cm=None) -> emit.EmitWorld:
    """S1-S11 of one world into ``out``: the signal world's stage code, whatever the preset (a null marker in the
    preset makes S11 write the null's truth)."""
    out = Path(out)
    backbone.build(params, seed, out)
    sequences.build(out, pfam_hmm=pfam_hmm)
    intergenic.build(out, rfam_cm=rfam_cm)
    mobile.build(out, pfam_hmm=pfam_hmm, rfam_cm=rfam_cm)
    fusions.build(out, pfam_hmm=pfam_hmm)
    assembly.build(out)
    sources.build(out)
    marker = params.get("null")
    if marker is None:
        return emit.build(out)
    return emit.build(out, null={"world": marker["world"], "null_seed": int(seed), "relabel": [],
                                 "scope": scope_n0(out, marker)})


def build(null: str, params: dict | str, seed: int, out, *, pfam_hmm=None, rfam_cm=None) -> emit.EmitWorld:
    """An N0 world (N0-full or N0-strict) of the signal preset ``params`` at ``seed``, S1-S11, into ``out``."""
    p = backbone.load_preset(params) if isinstance(params, str) else params
    return build_stages(null_preset(p, null), seed, out, pfam_hmm=pfam_hmm, rfam_cm=rfam_cm)


def _link_tree(src: Path, dst: Path) -> None:
    def link(a, b):
        try:
            os.link(a, b)
        except OSError:
            shutil.copy2(a, b)
    shutil.copytree(src, dst, copy_function=link)


def build_n1(signal, out, *, null_seed: int | None = None) -> emit.EmitWorld:
    """The N1 world of a signal world whose S1-S10 are built: its stage files shared (hard links, else copies), the
    relabel drawn, and S11 written with it (module docstring)."""
    signal, out = Path(signal), Path(out)
    bb = json.loads((signal / "backbone" / "backbone.json").read_text(encoding="utf-8"))
    if bb["format"] != backbone.FORMAT:
        raise ValueError(f"{signal}: format {bb['format']!r}, expected {backbone.FORMAT!r}")
    if "null" in bb["preset"]:
        raise ValueError(f"{signal} is a null world ({bb['preset']['null']['world']}); N1 relabels a signal world")
    if out.exists():
        raise FileExistsError(f"{out} exists; a world is written once")
    seed = int(bb["seed"])
    ns = seed if null_seed is None else seeds.check_seed(null_seed)
    tmp = out.with_name(out.name + ".partial")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    for d in SHARED:
        _link_tree(signal / d, tmp / d)
    tmp.rename(out)
    rows = relabel(out, bb["preset"]["nulls"], ns)
    return emit.build(out, null={"world": "N1", "null_seed": ns, "relabel": rows, "scope": scope_n1(out, rows)})


# --- N1 relabel -------------------------------------------------------------------------------------------------

def match_relabel(modules, profiles: dict, transfers: dict, pool, n_genomes: int, cfg: dict, *,
                  null_seed: int) -> list[RelabelRow]:
    """The relabel of `module` members: ``modules`` is [(module_id, [family, ...])] in order; ``profiles`` maps a
    family to its frozenset of genomes, ``transfers`` to its transfer count; ``pool`` the families that may replace
    one. For each member in order, a free pool family within floor(``max_hamming_fraction`` x ``n_genomes``)
    genomes and within ``transfer_tolerance`` of its transfers, drawn uniformly from the sorted eligible ones with
    stream ``N1/relabel/<family>``; else unmatched."""
    max_h = math.floor(float(cfg["max_hamming_fraction"]) * int(n_genomes) + 1e-9)
    tol = float(cfg["transfer_tolerance"])
    free = sorted(set(pool))
    out = []
    for mod, fams in modules:
        for i, a in enumerate(fams):
            pa, ta = profiles.get(a, frozenset()), int(transfers.get(a, 0))
            ok = [x for x in free if len(pa ^ profiles.get(x, frozenset())) <= max_h
                  and abs(int(transfers.get(x, 0)) - ta) <= tol * ta + 1e-9]
            if not ok:
                out.append(RelabelRow(mod, i, a, "", "unmatched", -1, len(pa), -1, ta, -1))
                continue
            x = ok[int(seeds.generator(null_seed, f"N1/relabel/{a}").integers(len(ok)))]
            free.remove(x)
            px = profiles.get(x, frozenset())
            out.append(RelabelRow(mod, i, a, x, "relabeled", len(pa ^ px), len(pa), len(px), ta,
                                  int(transfers.get(x, 0))))
    return out


def profiles_and_transfers(world: Path) -> tuple[dict[str, frozenset], dict[str, int], int]:
    """Per host family: the set of genomes with a copy in the S8 genomes (a fused gene counts for both families),
    the transfer events of the backbone, and the number of host genomes."""
    world = Path(world)
    present = collections.defaultdict(set)
    with open(world / "assembly" / "genes.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for f in r["family_id"].split("+"):
                present[f].add(r["genome_id"])
    transfers = collections.Counter()
    with open(world / "backbone" / "events.tsv", newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"].startswith("transfer"):
                transfers[f"F{int(r['family']):06d}"] += 1
    n = len(mobile._tsv(world / "backbone" / "genomes.tsv"))
    return {f: frozenset(s) for f, s in present.items()}, dict(transfers), n


def relabel(world: Path, cfg: dict, null_seed: int) -> list[RelabelRow]:
    check_config(cfg)
    world = Path(world)
    fams = mobile._tsv(world / "backbone" / "families.tsv")
    groups = mobile._tsv(world / "backbone" / "groups.tsv")
    fid = {int(r["family"]): r["family_id"] for r in fams}
    modules = [(g["name"], [fid[int(x)] for x in g["members"].split(";")]) for g in groups if g["kind"] == "planted"]
    kinds = set(cfg["n1"]["kinds"])
    pool = [r["family_id"] for r in fams if r["kind"] in kinds]
    profiles, transfers, n = profiles_and_transfers(world)
    return match_relabel(modules, profiles, transfers, pool, n, cfg["n1"], null_seed=null_seed)


# --- what the truth says: removed items and the rules of null evaluation ----------------------------------------

def _groups(world: Path) -> tuple[list[dict], dict[int, str]]:
    fid = {int(r["family"]): r["family_id"] for r in mobile._tsv(world / "backbone" / "families.tsv")}
    return mobile._tsv(world / "backbone" / "groups.tsv"), fid


def scope_n0(world, marker: dict) -> list[ScopeRow]:
    world = Path(world)
    rm = marker["removed"]
    out = []
    for i in range(int(rm["modules.count"])):
        out.append(ScopeRow("planted_module", f"mod{i:02d}", "", "removed",
                            "§4.1: zero planted groups; the group is not declared, so its families do not exist"))
    for j in range(int(rm["profile_only.count"])):
        out.append(ScopeRow("profile_only", f"pro{j:02d}", "", "removed", "§4.1: zero planted groups"))
    for i, m in enumerate(rm["decoys.partner_modules"]):
        out.append(ScopeRow("decoy", f"dec{i}", "", "removed",
                            f"§4.1: no Recipients decoys (partner mod{int(m):02d}.0 is removed with its module)"))
    for e in range(int(rm["fusions.nogd"])):
        out.append(ScopeRow("nogd", f"NG{e:02d}", "", "removed", "§4.1: no NOGD (it displaces a module member)"))
    for v in range(len(rm.get("vehicles.cargo", []))):
        out.append(ScopeRow("vehicle", f"veh{v}", "", "removed", "§4.1 N0-strict: no vehicles"))
    for a in range(int(rm.get("crispr.arrays", 0))):
        out.append(ScopeRow("crispr_cas_locus", f"CR{a:02d}", "", "removed",
                            "§4.1 N0-strict: no Cas group; the array is planted with its operon as one locus (S6), "
                            "so both go"))
    groups, fid = _groups(world)
    for g in groups:
        mem = [fid[int(x)] for x in g["members"].split(";")]
        if g["kind"] == "vehicle":
            out.append(ScopeRow("vehicle_chassis", g["name"], f"{mem[0]};{mem[-1]}", "excluded",
                                "§3.3, §4.1 N0-full: chassis pairs are excluded from null evaluation"))
            out.append(ScopeRow("vehicle_cargo", g["name"], ";".join(mem), "null_trap",
                                "pairs with a cargo family: must_not_call (truth_traps mge_cargo); cargo always "
                                "shares chassis flanks, so a count by vehicle-independent locus collapses them"))
        elif g["kind"] == "plasmid":
            out.append(ScopeRow("plasmid", g["name"], ";".join(mem), "null",
                                "decision (session 9): a plasmid's genes share their flanks wherever it goes, so a "
                                "count by locus collapses them; §5.3 keeps plasmids in N0-strict, whose every pair "
                                "is null"))
    cas = collections.defaultdict(list)
    for r in mobile._tsv(world / "mobile" / "cas.tsv"):
        cas[r["cas_module_id"]].append((int(r["order_index"]), r["family_id"]))
    for mod in sorted(cas):
        out.append(ScopeRow("cas_operon", mod, ";".join(f for _, f in sorted(cas[mod])), "excluded",
                            "decision (session 9): N0-full keeps the CRISPR-Cas loci of S6 and their cas_operon rows; "
                            "pairs within a locus's Cas families are excluded, as chassis pairs are"))
    out.append(ScopeRow("all_pairs", "*", "*", "null", "every other pair of families, phage and prophage genes "
                                                       "included"))
    return out


def scope_n1(world, rows: list[RelabelRow]) -> list[ScopeRow]:
    out = []
    by = collections.defaultdict(list)
    for r in rows:
        by[r.module_id].append(r)
    for mod, rs in by.items():
        moved = [r for r in rs if r.status == "relabeled"]
        labels = [r.replacement_family_id if r.status == "relabeled" else r.family_id for r in rs]
        if moved:
            out.append(ScopeRow("relabeled_module", mod, ";".join(labels), "evaluated",
                                "§4.1 N1: pairs of the module's labels with at least one relabeled member; their "
                                "p-values must be uniform (MR-18)"))
            out.append(ScopeRow("planted_unit", mod, ";".join(r.family_id for r in rs), "excluded",
                                "the planted unit under labels no module lists: it keeps its adjacency"))
        else:
            out.append(ScopeRow("unmatched_module", mod, ";".join(labels), "not_evaluated",
                                "no member has a replacement: the module is the signal world's"))
        for r in rs:
            if r.status == "unmatched":
                out.append(ScopeRow("unmatched_family", mod, r.family_id, "not_evaluated",
                                    "no free family within the Hamming and transfer bounds: the label stays; a pair "
                                    "of two unmatched members keeps its adjacency"))
    for r in mobile._tsv(Path(world) / "fusions" / "events.tsv"):
        if r["kind"] == "nogd":
            out.append(ScopeRow("nogd", r["event_id"], f"{r['family_a']};{r['replacement_family']}", "removed",
                                "N1 truth_module has no NOGD rows: they describe the planted unit"))
    out.append(ScopeRow("all_pairs", "*", "*", "not_evaluated", "§4.1 N1: only the relabeled pairs are evaluated"))
    return out


__all__ = ["NULLS", "N0_CHANGES", "SHARED", "RelabelRow", "ScopeRow", "check_config", "null_preset", "n2_seeds",
           "build_stages", "build", "build_n1", "match_relabel", "profiles_and_transfers", "relabel", "scope_n0",
           "scope_n1"]
