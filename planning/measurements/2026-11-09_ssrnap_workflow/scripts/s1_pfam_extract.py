"""Step 1a: find the single-subunit RNAP models in Pfam-A and write them as standalone .hmm files.

Also writes DNA_pol_A (the Pol I-fold relative) for later false-positive checks.
"""
import gzip
import sys
import time

import pyhmmer

PFAM = "/refdb/pfam/38.2/Pfam-A.hmm.gz"
OUT = "/tmp/if2453/ssrnap/hmm"
WANT = {"RNA_pol", "RPOL_N", "DNA_pol_A", "RNA_pol_Rpb1_1", "RNA_pol_Rpb1_2", "RNA_pol_Rpb2_6",
        "N4_RNAP_helical", "vRNAP_dom", "vRNAP_dom_2", "vRNAP_plug", "Phage_RNApol", "nvRNAP_subunit"}
PANEL = open("/tmp/if2453/ssrnap/hmm/panel.hmm", "wb")

t0 = time.time()
import os
os.makedirs(OUT, exist_ok=True)
n = 0
print("acc\tname\tM\tGA\tdesc")
with gzip.open(PFAM, "rb") as gz, pyhmmer.plan7.HMMFile(gz) as hf:
    for hmm in hf:
        n += 1
        name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
        desc = hmm.description or ""
        desc = desc.decode() if isinstance(desc, bytes) else desc
        acc = hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession
        if "RNA polymerase" in desc or "RNA_pol" in name or name in WANT:
            ga = hmm.cutoffs.gathering if hmm.cutoffs.gathering_available() else None
            print(f"{acc}\t{name}\t{hmm.M}\t{ga}\t{desc}")
        if name in WANT:
            with open(f"{OUT}/{name}.hmm", "wb") as fh:
                hmm.write(fh)
            hmm.write(PANEL)
PANEL.close()
print(f"# {n} models scanned in {time.time() - t0:.1f} s", file=sys.stderr)
