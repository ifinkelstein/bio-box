"""Explain the terminase count: slots with any PF03237.22 domain vs slots whose canonical family token is it.

    python check_term.py RELEASE
"""

import json
import sys
import time
from pathlib import Path

from biobox import catalog

rel = sys.argv[1]
W = Path("/tmp/if2453/ladder")
con = catalog.open_catalog(W / f"root/catalog/release={rel}.duckdb", local_dir=W / "catalog_cache/local", threads=4)
t = time.perf_counter()
any_dom = con.execute("""WITH p AS (SELECT DISTINCT protein_uid FROM l1_domains WHERE model_acc = 'PF03237.22')
    SELECT count(*), count(DISTINCT f.protein_uid) FROM feature_protein f JOIN p USING (protein_uid)""").fetchall()[0]
canon = con.execute(f"""WITH p AS (SELECT DISTINCT protein_uid FROM l1_domains WHERE model_acc = 'PF03237.22'),
    s AS (SELECT f.feature_uid FROM feature_protein f JOIN p USING (protein_uid))
    SELECT d.key, count(*) AS n FROM s JOIN token_postings tp ON tp.feature_uid = s.feature_uid
      AND tp.release = '{rel}' AND tp.token_id >= (1::UBIGINT << 61) AND tp.token_id < (2::UBIGINT << 61)
    JOIN token_dict d ON d.token_id = tp.token_id AND d.release = '{rel}'
    GROUP BY 1 ORDER BY n DESC, 1 LIMIT 8""").fetchall()
names = dict(con.execute("SELECT DISTINCT model_acc, model_name FROM l1_domains").fetchall())
print(json.dumps({"release": rel, "slots_with_any_PF03237": any_dom[0], "proteins": any_dom[1],
                  "canonical_family_of_those_slots": [[k, names.get(k), n] for k, n in canon],
                  "wall_s": round(time.perf_counter() - t, 2)}, indent=1))
