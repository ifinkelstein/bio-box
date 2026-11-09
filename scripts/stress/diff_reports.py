"""Compare two build_s0.py reports: manifest hash and per-file sha256.

    python diff_reports.py A.json B.json
"""
import json
import sys

a, b = (json.load(open(p)) for p in sys.argv[1:3])
print("hosts:", a["machine"], a["python"], "|", b["machine"], b["python"])
print("manifest_hash:", a["manifest_hash"], b["manifest_hash"], "EQUAL" if a["manifest_hash"] == b["manifest_hash"] else "DIFFERENT")
fa, fb = a["file_sha256"], b["file_sha256"]
same = [p for p in fa if fb.get(p) == fa[p]]
diff = [p for p in fa if p in fb and fb[p] != fa[p]]
print(f"files: {len(fa)} vs {len(fb)}; identical {len(same)}; different {len(diff)}; only A {len(set(fa) - set(fb))}; only B {len(set(fb) - set(fa))}")
for p in diff:
    print("  DIFF", p)
print("timings:", a["timings_s"], "|", b["timings_s"])
print("peak_rss_mb:", a["peak_rss_mb"], "|", b["peak_rss_mb"])
