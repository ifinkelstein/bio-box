"""Print a search run's record (`run.parquet`) as JSON, scalar columns only: r50_runrec.py RUN_DIR. Read only."""
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

row = pq.read_table(Path(sys.argv[1]) / "run.parquet").to_pylist()[0]
print(json.dumps({k: (v.hex() if isinstance(v, bytes) else v) for k, v in row.items()
                  if not isinstance(v, (list, dict))}, default=str))
