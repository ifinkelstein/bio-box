"""Read the header fields of every model of an Infernal 1.1 CM file into a TSV (card 0-08 design note).

    python -I cmheaders.py <Rfam.cm> <out.tsv>

Per model: ACC, NAME, CLEN, W, GA, TC, NC, NSEQ (seed sequences), the number of MATP nodes in the CM body
(consensus base pairs; 0 means cmsearch runs the model in HMM-only mode by default), MAXL of the filter
HMM, and DESC. Reads line by line; holds one model's fields at a time.
"""

import sys

FIELDS = ["acc", "name", "clen", "w", "ga", "tc", "nc", "nseq", "bps", "maxl", "desc"]


def main(cm, out):
    rows, head, where = [], {}, None
    with open(cm, errors="replace") as fh:
        for line in fh:
            if line.startswith("INFERNAL1/"):
                head, where = {"bps": 0}, "cmhead"
                continue
            if line.startswith("HMMER3/"):
                where = "hmmhead"
                continue
            if where is None:
                continue
            if where == "cmhead":
                parts = line.rstrip("\n").split(None, 1)
                if not parts:
                    continue
                k = parts[0]
                v = parts[1].strip() if len(parts) > 1 else ""
                if k == "CM":
                    where = "cmbody"
                elif k in ("NAME", "ACC", "CLEN", "W", "NSEQ"):
                    head[k.lower()] = v.split()[0]
                elif k == "DESC":
                    head["desc"] = v.replace("\t", " ")
                elif k in ("GA", "TC", "NC"):
                    head[k.lower()] = v.split()[0]
            elif where == "cmbody":
                s = line.lstrip()
                if s.startswith("[ MATP"):
                    head["bps"] += 1
            elif where == "hmmhead":
                parts = line.split(None, 1)
                if parts and parts[0] == "MAXL":
                    head["maxl"] = parts[1].split()[0]
                elif parts and parts[0] == "HMM":
                    rows.append(head)
                    head, where = {}, None
    rows.sort(key=lambda r: r["acc"])
    with open(out, "w") as fo:
        fo.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fo.write("\t".join(str(r.get(f, "")) for f in FIELDS) + "\n")
    print(len(rows), "models;", sum(1 for r in rows if r["bps"] == 0), "without base pairs;",
          sum(1 for r in rows if not r.get("ga")), "without GA")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
