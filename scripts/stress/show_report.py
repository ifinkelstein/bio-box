"""Compact view of a build_s0.py report.

    python show_report.py REPORT.json
"""
import json
import sys

if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    sys.exit(0 if sys.argv[1:2] in (["-h"], ["--help"]) else 2)

r = json.load(open(sys.argv[1]))
print(r["machine"], r["python"], r["release"], "manifest", r["manifest_hash"], "bp", r["total_bp"], "peak_rss_mb", r["peak_rss_mb"])
print("timings", r["timings_s"])
print("steps", json.dumps(r["steps"]))
genes = r["steps"]["genecall"]["genes"]
for name, t in r["tables"].items():
    rows = t.get("rows") or 0
    print(f"{name:28s} files {t['files']:3d} rows {rows:10d} bytes {t['bytes']:11d} B/row {t['bytes'] / rows if rows else 0:8.2f} "
          f"B/gene {t['bytes'] / genes:8.2f} B/Mbp {t['bytes'] / (r['total_bp'] / 1e6):12.0f}")
