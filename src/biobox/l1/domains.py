"""Library domain scans with pyhmmer (operation O-SCAN).

Every distinct protein of a gene-calling version is searched with every model of
one HMM library, using the library's own bit-score cutoffs (Pfam GA, NCBIfam TC)
and a fixed search space (Z = domZ = 1e8, SEARCH_PLAN SD-04, DEC-04), so stored
E-values do not change as the corpus grows. Only included domains are kept.

Writes `<root>/l1/domains/library=<name>-<version>/stage=<S>/part-0.parquet`.
"""

from __future__ import annotations

import gzip
import io
import tarfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyhmmer

from biobox import ledger
from biobox.registry.schemas import table

Z_PROT = 1e8
DOM_Z = 1e8


def iter_hmms(path: Path):
    """Read HMMs from a .hmm, .hmm.gz, or a .tgz of .HMM files (NCBIfam), in file order."""
    path = Path(path)
    name = path.name
    if name.endswith((".tgz", ".tar.gz")):
        with tarfile.open(path, "r:gz") as tar:
            for member in sorted((m for m in tar.getmembers() if m.isfile() and m.name.upper().endswith(".HMM")),
                                 key=lambda m: m.name):
                data = tar.extractfile(member).read()
                with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as hf:
                    yield from hf
        return
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rb") as fh, pyhmmer.plan7.HMMFile(fh) as hf:
        yield from hf


def load_proteins(root: Path, gene_call_version: str, stage: str) -> dict[bytes, str]:
    files = sorted((Path(root) / "l0seq/protein_seq").glob(
        f"gene_call_version={gene_call_version}/stage={stage}/protein_prefix=*/part-0.parquet"))
    out: dict[bytes, str] = {}
    for f in files:
        t = pq.read_table(f, columns=["protein_uid", "seq"])
        out.update(zip(t.column("protein_uid").to_pylist(), t.column("seq").to_pylist()))
    return out


def _decode(x) -> str:
    return x.decode() if isinstance(x, bytes) else (x or "")


def scan(
    proteins: dict[bytes, str],
    hmm_path: Path,
    cutoff: str,
    cpus: int = 1,
) -> list[dict]:
    """Search `proteins` with every model in `hmm_path`; return included domains as rows."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    uids = sorted(proteins)
    block = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=u.hex().encode(), sequence=proteins[u].replace("*", "")) for u in uids
    ).digitize(alphabet)
    kwargs = {"cpus": cpus, "Z": Z_PROT, "domZ": DOM_Z}
    if cutoff != "none":
        kwargs["bit_cutoffs"] = cutoff
    rows = []
    for hits in pyhmmer.hmmsearch(iter_hmms(hmm_path), block, **kwargs):
        q = hits.query
        acc, name, mlen = _decode(q.accession) or _decode(q.name), _decode(q.name), q.M
        for hit in hits:
            if not hit.included:
                continue
            for dom in hit.domains:
                if not dom.included:
                    continue
                a = dom.alignment
                rows.append({
                    "protein_uid": bytes.fromhex(_decode(hit.name)),
                    "model_acc": acc,
                    "model_name": name,
                    "model_len": mlen,
                    "seq_bitscore": hit.score,
                    "seq_evalue_zref": hit.evalue,
                    "dom_bitscore": dom.score,
                    "dom_ievalue_zref": dom.i_evalue,
                    "env_start": dom.env_from - 1,
                    "env_end": dom.env_to,
                    "ali_start": a.target_from - 1,
                    "ali_end": a.target_to,
                    "hmm_from": a.hmm_from - 1,
                    "hmm_to": a.hmm_to,
                    "cutoff": cutoff,
                    "z": Z_PROT,
                    "dom_z": DOM_Z,
                })
    rows.sort(key=lambda r: (r["protein_uid"], r["env_start"], r["model_acc"]))
    return rows


def run(root: Path, stage: str, gene_call_version: str, library: str, library_version: str,
        hmm_path: Path, cutoff: str, cpus: int = 1) -> dict:
    proteins = load_proteins(root, gene_call_version, stage)
    rows = scan(proteins, hmm_path, cutoff, cpus=cpus)
    out = Path(root) / f"l1/domains/library={library}-{library_version}/stage={stage}/part-0.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows, schema=table("l1_domains").schema), out,
                   compression="zstd", row_group_size=122_880)
    hit = len({r["protein_uid"] for r in rows})
    lib = f"{library}-{library_version}".encode()
    ledger.write(root, stage, f"domains-{library}-{library_version}",
                 [ledger.row(lib, "library", len(proteins), hit),
                  ledger.row(lib, "library", len(proteins), hit, n_rejected=len(proteins) - hit, reason="no_hit")])
    return {"proteins": len(proteins), "domains": len(rows),
            "proteins_with_hit": hit,
            "models_with_hit": len({r["model_acc"] for r in rows})}
