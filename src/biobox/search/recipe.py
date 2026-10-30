"""Search recipes: parse, validate, normalize, hash (SEARCH_PLAN §4; S-04, S-12).

`load` checks a recipe on paper: every key known, every value of the right kind, a positive
control present (S-04: a recipe without one is refused), thresholds in bits (SD-04). It fills
defaults and returns the normalized recipe; `recipe_hash` is the sha256 of its canonical
JSON, so two spellings of one recipe hash alike.

`executable` lists what this build cannot yet run (kinds, targets, tools, options), so a
recipe can be valid on paper (S-12) before the stage that runs it. The runner refuses a
recipe with any such gap at submit, along with unknown model IDs and seeds or positives
that the release does not hold (`runner.resolve`).
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from pathlib import Path

import yaml

from biobox.search.nbr import Policy, PolicyError

SEED_KINDS = ("sequence", "protein_uid", "region_uid", "hmm", "msa", "protein_accession", "cm", "structure",
              "locus", "region")
TARGETS = ("unique", "dark_orf", "igr", "nr50", "nr90", "segments", "contigs", "spacers")
TOOLS = ("phmmer", "hmmsearch", "jackhmmer", "nhmmer", "cmsearch", "mmseqs", "diamond", "hhblits", "foldseek",
         "blastn", "tblastn")
STOP_RULES = ("no_new_reps", "context_lost", "drift_alarm")
OUTPUTS = ("hits", "msa_per_round", "hmm_per_round", "neighborhoods", "summary")

# what this build runs (S0, ST0)
RUN_SEED_KINDS = ("sequence", "protein_uid", "hmm", "msa")
RUN_POSITIVE_KINDS = ("sequence", "protein_uid", "region_uid")
RUN_TARGETS = {"unique": "protein", "dark_orf": "dark_orf"}
RUN_TOOLS = ("phmmer", "hmmsearch", "jackhmmer")

TOP_KEYS = {"recipe", "written_against", "seed", "split_domains", "target", "tool", "iterations", "z", "curation",
            "context", "partial", "max_output_gb", "controls", "outputs", "description", "manifest"}
_NAME = re.compile(r"^[a-z0-9][a-z0-9_]*$")


class RecipeError(ValueError):
    """A recipe refused on paper or at submit."""


def _keys(d: dict, allowed: set, where: str) -> None:
    if not isinstance(d, dict):
        raise RecipeError(f"{where}: expected a mapping")
    extra = set(d) - allowed
    if extra:
        raise RecipeError(f"{where}: unknown keys {sorted(extra)}")


def _seed(s, where: str) -> dict:
    if isinstance(s, dict) and "kind" not in s and len(s) == 1:  # shorthand: {protein_accession: MJ0613}
        (k, v), = s.items()
        s = {"kind": k, "value": v}
    _keys(s, {"kind", "value", "label", "source"}, where)
    if s.get("kind") not in SEED_KINDS:
        raise RecipeError(f"{where}: kind {s.get('kind')!r} not one of {SEED_KINDS}")
    if s.get("value") in (None, ""):
        raise RecipeError(f"{where}: value missing")
    out = {"kind": s["kind"], "value": str(s["value"])}
    if s["kind"] == "sequence":
        seq = re.sub(r"\s+", "", out["value"]).upper().rstrip("*")
        if not re.fullmatch(r"[A-Z]+", seq):
            raise RecipeError(f"{where}: sequence has characters outside A-Z")
        out["value"] = seq
    if s["kind"] in ("protein_uid", "region_uid"):
        if not re.fullmatch(r"[0-9a-f]{32}", out["value"]):
            raise RecipeError(f"{where}: {s['kind']} must be 32 lowercase hex digits")
    if "label" in s:
        out["label"] = str(s["label"])
    if "source" in s:
        src = s["source"]
        _keys(src, {"paper_id", "evidence", "database", "release", "entry"}, f"{where}.source")
        if "paper_id" in src and "evidence" not in src:
            raise RecipeError(f"{where}.source: a paper seed needs evidence (file:line)")
        if "database" in src and not {"release", "entry"} <= set(src):
            raise RecipeError(f"{where}.source: a database seed needs release and entry")
        out["source"] = dict(src)
    return out


def normalize(raw: dict) -> dict:
    """Validate on paper and fill defaults; raises RecipeError."""
    if not isinstance(raw, dict):
        raise RecipeError("recipe: expected a mapping")
    _keys(raw, TOP_KEYS, "recipe")
    r = copy.deepcopy(raw)
    for k in ("recipe", "written_against", "seed", "target", "tool"):
        if k not in r:
            raise RecipeError(f"recipe: {k} is required")
    if not _NAME.match(str(r["recipe"])):
        raise RecipeError("recipe: name must be lowercase letters, digits and underscores")
    seeds = r["seed"] if isinstance(r["seed"], list) else [r["seed"]]
    if not seeds:
        raise RecipeError("seed: at least one")
    out = {"recipe": str(r["recipe"]), "written_against": str(r["written_against"]),
           "seed": [_seed(s, f"seed[{i}]") for i, s in enumerate(seeds)],
           "split_domains": bool(r.get("split_domains", False))}
    if r["target"] not in TARGETS:
        raise RecipeError(f"target {r['target']!r} not one of {TARGETS}")
    out["target"] = r["target"]
    tool = r["tool"] if isinstance(r["tool"], dict) else {"name": r["tool"]}
    _keys(tool, {"name", "engine", "version", "options"}, "tool")
    if tool.get("name") not in TOOLS:
        raise RecipeError(f"tool {tool.get('name')!r} not one of {TOOLS}")
    if tool.get("version", "from_image") != "from_image":
        raise RecipeError("tool.version: only from_image (the pinned container version)")
    out["tool"] = {"name": tool["name"], "engine": tool.get("engine", "pyhmmer"), "version": "from_image",
                   "options": dict(tool.get("options", {}))}
    it = r.get("iterations", {})
    _keys(it, {"max", "inclusion_bits", "domain_bits", "inclusion_evalue", "stop", "alarm_fraction", "alarm_min",
               "context_retention"}, "iterations")
    if "inclusion_evalue" in it:
        raise RecipeError("iterations.inclusion_evalue: thresholds are stored as bits (SD-04); convert the paper's "
                          "E-value at its database size and give inclusion_bits")
    if "inclusion_bits" not in it:
        raise RecipeError("iterations.inclusion_bits is required (thresholds are bits, SD-04)")
    stop = list(it.get("stop", ["no_new_reps", "drift_alarm"]))
    bad = set(stop) - set(STOP_RULES)
    if bad:
        raise RecipeError(f"iterations.stop: unknown rules {sorted(bad)}")
    mx = int(it.get("max", 1))
    if mx < 1:
        raise RecipeError("iterations.max must be >= 1")
    if tool["name"] in ("phmmer", "hmmsearch") and mx != 1:
        raise RecipeError(f"{tool['name']} runs one round; use jackhmmer to iterate")
    bits = float(it["inclusion_bits"])
    out["iterations"] = {"max": mx, "inclusion_bits": bits, "domain_bits": float(it.get("domain_bits", bits / 2)),
                         "stop": sorted(stop), "alarm_fraction": float(it.get("alarm_fraction", 0.5)),
                         "alarm_min": int(it.get("alarm_min", 3)),
                         "context_retention": float(it.get("context_retention", 0.5))}
    if r.get("z", "release_default") != "release_default":
        raise RecipeError("z: only release_default (Z_prot = domZ = 1e8, Z_nt = 1e5 Mb; SD-04)")
    out["z"] = "release_default"
    cur = r.get("curation", {})
    _keys(cur, {"min_query_coverage", "required_motif", "foreign_domain_alarm", "promiscuity"}, "curation")
    out["curation"] = {"min_query_coverage": float(cur.get("min_query_coverage", 0.0)),
                       "required_motif": cur.get("required_motif"),
                       "foreign_domain_alarm": bool(cur.get("foreign_domain_alarm", True)),
                       "promiscuity": cur.get("promiscuity", "none")}
    if not 0 <= out["curation"]["min_query_coverage"] <= 1:
        raise RecipeError("curation.min_query_coverage must be in [0, 1]")
    if out["curation"]["promiscuity"] not in ("none", "down_weight", "blocklist"):
        raise RecipeError("curation.promiscuity: none, down_weight or blocklist")
    ctx = r.get("context", {})
    _keys(ctx, {"window", "strand", "edge", "operon_gap", "circular"}, "context")
    out["context"] = {"window": str(ctx.get("window", "10g")), "strand": ctx.get("strand", "any"),
                      "edge": str(ctx.get("edge", "truncate")), "operon_gap": ctx.get("operon_gap"),
                      "circular": ctx.get("circular", "auto")}
    try:
        Policy.parse(**out["context"])
    except PolicyError as e:
        raise RecipeError(f"context: {e}") from None
    out["partial"] = r.get("partial", "flag")
    if out["partial"] not in ("exclude", "include", "flag"):
        raise RecipeError("partial: exclude, include or flag (SD-20)")
    out["max_output_gb"] = float(r.get("max_output_gb", 25))
    if out["max_output_gb"] <= 0:
        raise RecipeError("max_output_gb must be > 0")
    ctl = r.get("controls")
    if not ctl or not ctl.get("positive"):
        raise RecipeError("controls.positive is required: a recipe without a positive control is refused (S-04)")
    _keys(ctl, {"positive", "negative"}, "controls")
    pos = ctl["positive"] if isinstance(ctl["positive"], list) else [ctl["positive"]]
    neg = ctl.get("negative", "shuffled_seed")
    if neg not in ("shuffled_seed", "none"):
        raise RecipeError("controls.negative: shuffled_seed or none")
    out["controls"] = {"positive": [_seed(p, f"controls.positive[{i}]") for i, p in enumerate(pos)], "negative": neg}
    outs = list(r.get("outputs", ["hits", "hmm_per_round", "msa_per_round", "neighborhoods", "summary"]))
    if set(outs) - set(OUTPUTS):
        raise RecipeError(f"outputs: unknown {sorted(set(outs) - set(OUTPUTS))}")
    out["outputs"] = sorted(outs)
    if "description" in r:
        out["description"] = str(r["description"])
    if "manifest" in r:
        out["manifest"] = str(r["manifest"])
    return out


def recipe_hash(recipe: dict) -> str:
    return hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load(path_or_text: Path | str) -> tuple[dict, str]:
    """(normalized recipe, original YAML text)."""
    text = Path(path_or_text).read_text() if isinstance(path_or_text, Path) else path_or_text
    try:
        raw = yaml.safe_load(text)
    except yaml.YAMLError as e:
        raise RecipeError(f"not YAML: {e}") from None
    return normalize(raw), text


def executable(recipe: dict) -> list[str]:
    """What this build cannot run in the recipe; empty means the runner can execute it."""
    gaps = []
    for s in recipe["seed"]:
        if s["kind"] not in RUN_SEED_KINDS:
            gaps.append(f"seed kind {s['kind']} (runs: {', '.join(RUN_SEED_KINDS)})")
    if len(recipe["seed"]) != 1:
        gaps.append("more than one seed (one profile per run at S0)")
    for p in recipe["controls"]["positive"]:
        if p["kind"] not in RUN_POSITIVE_KINDS:
            gaps.append(f"positive control kind {p['kind']} (runs: {', '.join(RUN_POSITIVE_KINDS)})")
    if recipe["target"] not in RUN_TARGETS:
        gaps.append(f"target {recipe['target']} (runs: {', '.join(RUN_TARGETS)})")
    if recipe["tool"]["name"] not in RUN_TOOLS or recipe["tool"]["engine"] != "pyhmmer":
        gaps.append(f"tool {recipe['tool']['name']}/{recipe['tool']['engine']} (runs: {', '.join(RUN_TOOLS)} on pyhmmer)")
    if recipe["split_domains"]:
        gaps.append("split_domains (needs seed envelopes as separate profiles)")
    if recipe["curation"]["required_motif"]:
        gaps.append("curation.required_motif (needs curated catalytic columns)")
    if recipe["curation"]["promiscuity"] != "none":
        gaps.append("curation.promiscuity (needs domain_promiscuity)")
    if recipe["tool"]["options"]:
        gaps.append("tool.options")
    return gaps
