"""Pod side: run the Prodigal C binary (single mode) on the golden prokaryotes (stdlib only).

    python3 prodigal_run.py GOLDEN_DIR OUT_DIR
"""
import gzip
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

CODES = {"ecoli_mg1655": 11, "bsub_168": 11, "mgenitalium": 4, "sr1_raac1": 25, "macetivorans_c2a": 11, "hvolcanii_ds2": 11,
         "sacidocaldarius_dsm639": 11, "mtb_h37rv": 11, "btheta_vpi5482": 11, "spyogenes_sf370": 11}
golden, out = Path(sys.argv[1]), Path(sys.argv[2])
for gid, code in CODES.items():
    src = golden / gid / f"{gid}.fna.gz"
    fna = out / f"{gid}.fna"
    with gzip.open(src, "rb") as fi, open(fna, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    sha = hashlib.sha256(fna.read_bytes()).hexdigest()
    subprocess.run(["nice", "-n", "19", "prodigal", "-i", str(fna), "-f", "gff", "-o", str(out / f"{gid}.gff"), "-g", str(code), "-q"],
                   check=True)
    n = sum(1 for ln in open(out / f"{gid}.gff") if "\tCDS\t" in ln)
    print(gid, code, sha, n, flush=True)
