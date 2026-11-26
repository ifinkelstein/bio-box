"""Card R-15: `fetch_to_tower.sh` takes an `ftp://` URL only where the `sources.yaml` form it matches is `ftp://`.

DEC-81 (5): a source with no https route may be fetched over `ftp://`. `ftp.cngb.org` refuses connections on port 443
(card R-15's note), so batch 7 (`EEMC1`, `GPB_PHAGE`) has no https route. The rule tested here: an ftp URL is matched
only against forms written `ftp://...`, an https URL only against `https://...` forms; curl is held to the job's own
protocol for the request and any redirect; `http://` and every other scheme stay refused, as do the https cases the
R-08 tests cover (those tests are not changed).

Same stub `curl` as tests/r1/test_r08_fetch_to_tower.py: it records its arguments and serves bytes from local files
named by a JSON map; no request leaves the machine. The last tests read the tree's `sources.yaml` the way the R-08 test
does (skipped where the file is absent) and check one real URL of each form of the batch-7 job lists.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("R-15")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "fetch_to_tower.sh"

CURL = """#!{python}
import json, os, sys
args = sys.argv[1:]
with open(os.environ["STUB_LOG"], "a") as f:
    f.write(json.dumps({{"args": args}}) + "\\n")
out, url = args[args.index("-o") + 1], args[-1]
spec = json.load(open(os.environ["STUB_MAP"])).get(url)
if spec is None:
    sys.exit(78)
with open(out, "wb") as f:
    f.write(open(spec, "rb").read())
"""

SOURCES = """# fixture
version: 1
sources:
  - name: zen
    url: https://zenodo.org/api/records/111/files/{file}/content
  - name: cngb_set
    url: https://db.cngb.org/data_resources/project/CNP0000001/
    url_ftp: ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS{sample}/CNA{assembly}/{file}
    url_ftp_meta: ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/{file}
  - name: https_only
    url: https://same.example.org/pub/{path}
  - name: plain_http
    url: https://plain.example.org/landing
    url_http: http://plain.example.org/files/{file}
"""

FTP_FA = "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS0000010/CNA0000020/BGI1.fasta.gz"
FTP_META = "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/metadata_CNP0000001_assembly.tsv"
HTTPS_ZEN = "https://zenodo.org/api/records/111/files/seqs.fasta/content"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()  # noqa: S324 - the provider's checksum, not a security use


def gz(data: bytes) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as f:
        f.write(data)
    return buf.getvalue()


SERVED = {
    FTP_FA: gz(b">BGI1 contig_1\nACGTACGT\n"),
    FTP_META: b"assembly\tsample\nCNA0000020\tCNS0000010\n",
    HTTPS_ZEN: b">a\nACGT\n",
}


@pytest.fixture
def world(tmp_path: Path) -> dict:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    p = bindir / "curl"
    p.write_text(CURL.format(python=sys.executable))
    p.chmod(0o755)
    src = tmp_path / "sources.yaml"
    src.write_text(SOURCES)
    blobs = tmp_path / "blobs"
    blobs.mkdir()
    served = {}
    for i, (url, data) in enumerate(SERVED.items()):
        (blobs / f"b{i}").write_bytes(data)
        served[url] = str(blobs / f"b{i}")
    mp = tmp_path / "map.json"
    mp.write_text(json.dumps(served))
    env = dict(os.environ)
    env.update({"PATH": f"{bindir}{os.pathsep}{env['PATH']}", "STUB_LOG": str(tmp_path / "calls.jsonl"),
                "STUB_MAP": str(mp), "TMPDIR": str(tmp_path)})
    return {"tmp": tmp_path, "src": src, "dest": tmp_path / "incoming", "env": env, "map": mp}


def jobs(w: dict, rows: list[tuple], name: str = "jobs.tsv") -> Path:
    p = w["tmp"] / name
    p.write_text("".join("\t".join(r) + "\n" for r in rows))
    return p


def run(w: dict, *lists: Path, args: tuple = (), sources: Path | None = None) -> subprocess.CompletedProcess:
    cmd = ["bash", str(SCRIPT), "--dest", str(w["dest"]), "--sources", str(sources or w["src"]), "--pause", "0",
           *args, *map(str, lists)]
    return subprocess.run(cmd, env=w["env"], capture_output=True, text=True, timeout=300)


def calls(w: dict) -> list[list[str]]:
    p = Path(w["env"]["STUB_LOG"])
    return [json.loads(line)["args"] for line in p.read_text().splitlines()] if p.exists() else []


def ftp_rows(with_md5: bool = True) -> list[tuple]:
    return [
        ("batch7_cngb/SET/isolates", "BGI1.fasta.gz", FTP_FA, md5(SERVED[FTP_FA]) if with_md5 else "-", "CNGBdb-open"),
        ("batch7_cngb/SET", "metadata_CNP0000001_assembly.tsv", FTP_META, "-", "CNGBdb-open"),
    ]


def test_an_ftp_url_of_an_ftp_form_lands_and_is_recorded(world):
    w = world
    r = run(w, jobs(w, ftp_rows()))
    assert r.returncode == 0, r.stdout + r.stderr
    fa = w["dest"] / "batch7_cngb" / "SET" / "isolates" / "BGI1.fasta.gz"
    meta = w["dest"] / "batch7_cngb" / "SET" / "metadata_CNP0000001_assembly.tsv"
    assert fa.read_bytes() == SERVED[FTP_FA] and meta.read_bytes() == SERVED[FTP_META]
    rows = [x.split("\t") for x in (fa.parent / "SOURCE_URLS.tsv").read_text().splitlines()]
    assert rows[1][:5] == ["BGI1.fasta.gz", FTP_FA, str(len(SERVED[FTP_FA])), sha(SERVED[FTP_FA]), "CNGBdb-open"]
    assert (fa.parent / "SHA256SUMS").read_text().splitlines() == [f"{sha(SERVED[FTP_FA])}  BGI1.fasta.gz"]
    assert "DONE\tok=2\tskip=0\tfail=0" in r.stdout
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]


def test_each_request_is_held_to_its_own_protocol_and_carries_no_identity(world):
    w = world
    rows = ftp_rows() + [("S1", "seqs.fasta", HTTPS_ZEN, md5(SERVED[HTTPS_ZEN]), "CC BY 4.0")]
    r = run(w, jobs(w, rows))
    assert r.returncode == 0, r.stdout + r.stderr
    seen = calls(w)
    assert [c[-1] for c in seen] == [FTP_FA, FTP_META, HTTPS_ZEN]
    for c in seen:
        want = "=ftp" if c[-1].startswith("ftp://") else "=https"
        assert c[c.index("--proto") + 1] == want and c[c.index("--proto-redir") + 1] == want, c
        assert c.count("--proto") == 1 and c.count("--proto-redir") == 1, c
        assert c[0] == "-fsSL" and c[c.index("-A") + 1] == "curl", c
        for bad in ("-H", "--header", "-u", "--user", "-e", "--referer", "--user-agent", "-b", "--cookie", "-n",
                    "--netrc", "--ftp-account", "-Q", "--quote"):
            assert bad not in c, (bad, c)
        assert not re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+", " ".join(c)), c


def test_check_names_the_ftp_entry_and_fetches_nothing(world):
    w = world
    r = run(w, jobs(w, ftp_rows()), args=("--check",))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f"match cngb_set {FTP_FA}" in r.stdout and f"match cngb_set {FTP_META}" in r.stdout
    assert "check: 2 URLs, all covered" in r.stdout
    assert calls(w) == [] and not w["dest"].exists()


def test_an_ftp_file_whose_md5_differs_is_removed(world):
    w = world
    rows = ftp_rows()
    rows[0] = rows[0][:3] + (md5(b"another"),) + rows[0][4:]
    r = run(w, jobs(w, rows))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "FAIL_MD5" in r.stderr
    assert not (w["dest"] / "batch7_cngb" / "SET" / "isolates" / "BGI1.fasta.gz").exists()
    assert "DONE\tok=1\tskip=0\tfail=1" in r.stdout


CROSSED = [
    # an https URL of a path only an ftp form lists
    "https://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS0000010/CNA0000020/BGI1.fasta.gz",
    "https://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/metadata_CNP0000001_assembly.tsv",
    # an ftp URL of a path only an https form lists
    "ftp://same.example.org/pub/a/b.fa",
    "ftp://zenodo.org/api/records/111/files/seqs.fasta/content",
    "ftp://db.cngb.org/data_resources/project/CNP0000001/",
]


@pytest.mark.parametrize("url", CROSSED)
def test_a_form_of_one_protocol_never_matches_a_url_of_the_other(world, url):
    w = world
    r = run(w, jobs(w, ftp_rows() + [("S9", "x.fa", url, "-", "x")]))
    assert r.returncode == 2, r.stdout + r.stderr
    assert "not covered by any sources.yaml entry" in r.stderr and url in r.stderr
    assert calls(w) == [] and not w["dest"].exists()


FTP_UNCOVERED = [
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000002/CNS0000010/CNA0000020/BGI1.fasta.gz",  # another project
    "ftp://ftp.cngb.org/pub/CNSA/data5/CNP0000001/CNS0000010/CNA0000020/BGI1.fasta.gz",  # another data volume
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS0000010/CNA0000020/x/BGI1.fasta.gz",  # one segment more
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS0000010/BGI1.fasta.gz",  # one segment less
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/SAM0000010/CNA0000020/BGI1.fasta.gz",  # not CNS
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0000001/CNS0000010/CNA0000020/",  # no file
    "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/",
    "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a/b.tsv",
    "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv?x=1",
    "ftp://ftp.cngb.org:2121/pub/CNSA/data7/public_info/CNP0000001/a.tsv",  # another port
    "ftp://anonymous:x@ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv",  # a login in the URL
    "ftp://ftp.cngb.org@evil.example/pub/CNSA/data7/public_info/CNP0000001/a.tsv",
    "ftp://ftp.cngb.org.evil.example/pub/CNSA/data7/public_info/CNP0000001/a.tsv",
    "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/../../CNP0000002/a.tsv",
]


@pytest.mark.parametrize("url", FTP_UNCOVERED)
def test_an_ftp_url_no_ftp_form_covers_is_refused_before_any_request(world, url):
    w = world
    r = run(w, jobs(w, ftp_rows() + [("S9", "x.fa", url, "-", "x")]))
    assert r.returncode == 2, r.stdout + r.stderr
    assert "not covered by any sources.yaml entry" in r.stderr and url in r.stderr
    assert calls(w) == [] and not w["dest"].exists()


OTHER_SCHEMES = [
    ("ftps://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv", "refused URL"),
    ("sftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv", "refused URL"),
    ("tftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv", "refused URL"),
    ("FTP://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv", "refused URL"),
    ("ftp:/ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/a.tsv", "refused URL"),
    ("file:///etc/passwd", "refused URL"),
    ("http://plain.example.org/files/a.fa", "not https"),  # listed as an http:// form, and still refused
]


@pytest.mark.parametrize("url,why", OTHER_SCHEMES)
def test_every_other_scheme_stays_refused(world, url, why):
    w = world
    r = run(w, jobs(w, ftp_rows() + [("S9", "x.fa", url, "-", "x")]))
    assert r.returncode == 2, r.stdout + r.stderr
    assert why in r.stderr and url in r.stderr
    assert calls(w) == [] and not w["dest"].exists()


def test_a_sources_file_with_only_ftp_forms_is_usable(world):
    w = world
    w["src"].write_text("version: 1\nsources:\n  - name: only_ftp\n"
                        "    url_ftp: ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0000001/{file}\n")
    r = run(w, jobs(w, [ftp_rows()[1]]), args=("--check",))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f"match only_ftp {FTP_META}" in r.stdout


def test_an_ftp_form_with_a_variable_host_matches_nothing(world):
    w = world
    w["src"].write_text("version: 1\nsources:\n  - name: open\n    url: ftp://{host}/pub/{path}\n")
    r = run(w, jobs(w, [ftp_rows()[1]]))
    assert r.returncode == 2 and calls(w) == [], r.stdout + r.stderr


# --- the tree's sources.yaml and the batch-7 job lists ------------------------------------------------------------------


def real_sources() -> Path:
    """The tree's sources.yaml; the test image (built from main, container/Dockerfile) does not ship it, so skip there."""
    p = ROOT.joinpath("sources.yaml")
    if not p.is_file():
        pytest.skip("sources.yaml is not in the image: this check runs from a checkout (ci and the pod's mirror)")
    return p


REAL_BATCH7 = {  # one URL of each form, from _batch7_cngb_urls.tsv on tower (copy read 2026-11-26)
    "ftp://ftp.cngb.org/pub/CNSA/data5/public_info/CNP0004979/metadata_CNP0004979_assembly.tsv": "tower_gpb_phage",
    "ftp://ftp.cngb.org/pub/CNSA/data5/CNP0004979/CNS0943843/CNA0072848/CPB0155.fa.gz": "tower_gpb_phage",
    "ftp://ftp.cngb.org/pub/CNSA/data7/public_info/CNP0007106/metadata_CNP0007106_sample_Microbial_sample.tsv":
        "tower_eemc1",
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0007106/CNS1317800/CNA0429852/BGI80.fasta.gz": "tower_eemc1",
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0007106/CNS1322970/CNA0464303/EEMC.bin.95.fasta.gz": "tower_eemc1",
}


def test_the_trees_sources_yaml_covers_one_url_of_each_batch_7_form(world):
    w = world
    rows = [("S1", f"f{i}.bin", u, "-", "x") for i, u in enumerate(REAL_BATCH7)]
    r = subprocess.run(["bash", str(SCRIPT), "--check", "--sources", str(real_sources()), str(jobs(w, rows))],
                       env=w["env"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    for u, entry in REAL_BATCH7.items():
        assert f"match {entry} {u}" in r.stdout, (u, r.stdout)
    assert calls(w) == []


REAL_REFUSED = [
    # the https form of a batch-7 path: ftp.cngb.org refuses port 443, and no https form is listed for it
    "https://ftp.cngb.org/pub/CNSA/data7/CNP0007106/CNS1317800/CNA0429852/BGI80.fasta.gz",
    "https://ftp.cngb.org/pub/CNSA/data5/CNP0004979/CNS0943843/CNA0072848/CPB0155.fa.gz",
    # a CNGB project that is not listed (ELEPHANT1 of _batch7_skipped.tsv, controlled access)
    "ftp://ftp.cngb.org/pub/CNSA/data5/CNP0004215/CNS0000001/CNA0000001/x.fa.gz",
    "ftp://ftp.cngb.org/pub/CNSA/data5/public_info/CNP0004215/metadata_CNP0004215_assembly.tsv",
    # the batch-7 project in another data volume, and a listing of its directory
    "ftp://ftp.cngb.org/pub/CNSA/data5/CNP0007106/CNS1317800/CNA0429852/BGI80.fasta.gz",
    "ftp://ftp.cngb.org/pub/CNSA/data7/CNP0007106/",
    # an ftp URL of an https form of the tree (an NCBI genome path)
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz",
]


@pytest.mark.parametrize("url", REAL_REFUSED)
def test_the_trees_sources_yaml_refuses_what_batch_7_does_not_list(world, url):
    w = world
    rows = [("S1", "f.bin", url, "-", "x")]
    r = subprocess.run(["bash", str(SCRIPT), "--check", "--sources", str(real_sources()), str(jobs(w, rows))],
                       env=w["env"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 2, r.stdout + r.stderr
    assert "not covered" in r.stderr
