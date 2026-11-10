"""Intended path, step B: recipes through biobox.search.runner against the shadow root's pack.

usage: s8_runner.py [names...]   (default: all recipes)
"""
import json
import subprocess
import sys
import time
import traceback
from pathlib import Path

from biobox.search import runner
from biobox.search.recipe import RecipeError

SH = Path("/tmp/if2453/ssrnap/shadow")
REL = "rS1.pilot1"
WORK = Path("/tmp/if2453/ssrnap/work")
T7 = "2eef240afefddcd2c8ac7660eb50e384"  # T7 gp1, found in s4_seeds.py
T7_SEQ = [l.strip() for l in open("/tmp/if2453/ssrnap/hmm/seeds.faa")][1]
COMMON = f"""written_against: {REL}
target: unique
controls: {{positive: [{{kind: protein_uid, value: {T7}, label: T7 gp1}}]}}
context: {{window: 5g, edge: truncate}}
"""
RECIPES = {
    "pfam_hmm": "recipe: ssrnap_pfam\nseed: [{hmm: PF00940}]\ntool: hmmsearch\n"
                "iterations: {max: 1, inclusion_bits: 30}\n" + COMMON,
    "seed_msa": "recipe: ssrnap_seed8\nseed: [{msa: /tmp/if2453/ssrnap/hmm/seeds.sto}]\ntool: hmmsearch\n"
                "iterations: {max: 1, inclusion_bits: 30}\n" + COMMON,
    "hmm_file": "recipe: ssrnap_file\nseed: [{hmm: /tmp/if2453/ssrnap/hmm/RNA_pol.hmm}]\ntool: hmmsearch\n"
                "iterations: {max: 1, inclusion_bits: 30}\n" + COMMON,
    "jack_t7": f"recipe: ssrnap_jack\nseed: [{{sequence: {T7_SEQ}}}]\ntool: jackhmmer\n"
               "iterations: {max: 3, inclusion_bits: 30, stop: [no_new_reps, drift_alarm, context_lost]}\n" + COMMON,
}
MEM = None if "--nomem" in sys.argv else 18  # --nomem: no RLIMIT_AS, as the CLI runs it
names = [a for a in sys.argv[1:] if not a.startswith("--")] or list(RECIPES)
out = {}
for name in names:
    text = RECIPES[name]
    Path(f"/tmp/if2453/ssrnap/out/recipe_{name}.yaml").write_text(text)
    t0 = time.time()
    try:
        row = runner.run(SH, text, REL, work_base=WORK, limits=runner.Limits(cpus=8, wall_s=1800, mem_gb=MEM))
    except RecipeError as e:
        print(f"== {name}: REFUSED after {time.time() - t0:.1f} s: {e}", flush=True)
        out[name] = {"refused": str(e)}
        continue
    except Exception:
        print(f"== {name}: EXCEPTION after {time.time() - t0:.1f} s", flush=True)
        traceback.print_exc()
        continue
    hits = runner.read_hits(SH, row["run_id"].hex())
    uids = {h["target_uid"].hex() for h in hits}
    print(f"== {name}: status {row['status']} error {row['error']} | wall {row['wall_s']:.1f} s (call {time.time() - t0:.1f} s), "
          f"worker CPU {row['cpu_min_max_pid']:.2f} min, peak RSS {(row['peak_rss_bytes'] or 0) / 1e9:.2f} GB | "
          f"target_size {row['target_size']} rounds {row['rounds']} final {row['final_round']} stop {row['stop_reason']} | "
          f"final hit proteins {len(uids)}", flush=True)
    print("   controls", row["controls_json"], "\n   warnings", row["warnings_json"], "\n   rounds", row["rounds_json"],
          flush=True)
    files = sorted(p.name for p in runner.run_dir(SH, row["run_id"]).iterdir())
    print("   run dir", runner.run_dir(SH, row["run_id"]), files, flush=True)
    out[name] = {"run_id": row["run_id"].hex(), "status": row["status"], "hits": sorted(uids)}
suffix = "_" + "_".join(names) if sys.argv[1:] else ""
json.dump(out, open(f"/tmp/if2453/ssrnap/out/runner_results{suffix}.json", "w"), indent=1)
if sys.argv[1:]:
    sys.exit(0)

# CLI on the first finished run: hits as TSV, neighborhoods as NDJSON
done = [v["run_id"] for v in out.values() if "run_id" in v]
if done:
    rid = done[0]
    for args, dest in ((["search", "hits", rid, "--root", str(SH), "--top_n", "0"], "cli_hits.tsv"),
                       (["search", "hits", rid, "--root", str(SH), "--json", "--top_n", "0"], "cli_hits.ndjson"),
                       (["nbr", "--root", str(SH), "--anchors", rid, "--window", "5g", "--json"], "cli_nbr.ndjson")):
        t0 = time.time()
        p = subprocess.run([sys.executable, "-m", "biobox.cli"] + args, capture_output=True, text=True)
        Path(f"/tmp/if2453/ssrnap/out/{dest}").write_text(p.stdout)
        print(f"CLI {' '.join(args[:2])}: exit {p.returncode}, {len(p.stdout.splitlines())} lines, "
              f"{time.time() - t0:.1f} s; stderr tail: {p.stderr.strip().splitlines()[-3:]}", flush=True)
