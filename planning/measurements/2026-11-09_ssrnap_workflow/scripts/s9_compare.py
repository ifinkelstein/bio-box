"""Compare the runner's final hit sets with the direct hmmsearch hit sets at the same threshold."""
import json

from common import OUT, by_query, read_domtbl

res = json.load(open(f"{OUT}/runner_results.json"))
res["jack_t7_nomem"] = json.load(open(f"{OUT}/runner_results_jack_t7.json"))["jack_t7"]
pfam = by_query(read_domtbl(f"{OUT}/pfam_RNA_pol.domtbl"))["RNA_pol"]
s8 = by_query(read_domtbl(f"{OUT}/seed8.domtbl"))["ssRNAP_seed8"]
direct = {"pfam_hmm": {u for u, h in pfam.items() if h["seq_score"] >= 30},
          "seed_msa": {u for u, h in s8.items() if h["seq_score"] >= 30}}
for name, r in res.items():
    if "hits" not in r:
        print(name, r)
        continue
    got = set(r["hits"])
    d = direct.get(name)
    if d is None:
        print(f"{name}: {r['status']}, {len(got)} hits; vs Pfam>=30: both {len(got & direct['pfam_hmm'])}, "
              f"runner only {len(got - direct['pfam_hmm'])}; vs seed8>=30: both {len(got & direct['seed_msa'])}, "
              f"runner only {len(got - direct['seed_msa'])}")
        union = direct["pfam_hmm"] | direct["seed_msa"]
        extra = got - union
        lens = sorted((pfam.get(u) or s8.get(u) or {}).get("tlen", -1) for u in extra)
        print(f"   vs union of both direct sets ({len(union)}): runner only {len(extra)} (lengths known for "
              f"{sum(1 for x in lens if x > 0)}: {lens[:40]}), union only {len(union - got)}")
    else:
        print(f"{name}: {r['status']}, runner {len(got)}, direct {len(d)}, both {len(got & d)}, "
              f"runner only {len(got - d)}, direct only {len(d - got)}")
