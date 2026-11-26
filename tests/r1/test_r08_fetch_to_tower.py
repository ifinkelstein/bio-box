"""Card R-08: `scripts/orchestrator/fetch_to_tower.sh`, the fetch script of `tower` with the `sources.yaml` rule.

A stub `curl` comes first on PATH, records its arguments and serves bytes from local files named by a JSON map; no
request leaves the machine, `tower` is not touched, and every run uses a temporary DEST and a temporary sources file.
The script is tested on what it must do (land, record, skip, verify) and on what it must refuse: a URL that no
`sources.yaml` entry covers, a list line with a character the shell could read, a second run on the same DEST.
The last test matches one real URL of every source of the 2026-11-25 job lists against the tree's `sources.yaml`.
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
import tarfile
from pathlib import Path

import pytest

pytestmark = pytest.mark.covers("R-08")

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "orchestrator" / "fetch_to_tower.sh"

CURL = """#!{python}
import json, os, sys
args = sys.argv[1:]
with open(os.environ["STUB_LOG"], "a") as f:
    f.write(json.dumps({{"args": args}}) + "\\n")
out, url = args[args.index("-o") + 1], args[-1]
spec = json.load(open(os.environ["STUB_MAP"])).get(url)
if spec in (None, "FAIL"):
    sys.exit(22)
with open(out, "wb") as f:
    if spec != "EMPTY":
        f.write(open(spec, "rb").read())
"""

SOURCES = """# fixture
version: 1
sources:
  - name: zen
    what: a Zenodo record
    url: https://zenodo.org/api/records/111/files/{file}/content
    note: >-
      url: https://not-a-form.example/x is folded text and is not a form
  - name: fig
    url: https://ndownloader.figshare.com/files/{file_id}
  - name: lit
    url: https://example.org/literal/one.tar.gz   # comment after the form
  - name: deep
    url: https://legacy.example.org/dl/{strain}/fasta/data/{file}
    url_stacks: https://stacks.example.org/file/druid:abc/{path}
  - name: prose
    url: not recorded when it was fetched. The page is https://elsewhere.example/page
  - name: quoted
    url: "https://quoted.example/x/{file}"
"""

ZEN = "https://zenodo.org/api/records/111/files/{}/content"
FIG = "https://ndownloader.figshare.com/files/{}"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def md5(data: bytes) -> str:
    return hashlib.md5(data).hexdigest()  # noqa: S324 - the provider's checksum, not a security use


def gz(data: bytes) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as f:
        f.write(data)
    return buf.getvalue()


def tar(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as t:
        for name, data in members.items():
            ti = tarfile.TarInfo(name)
            ti.size = len(data)
            ti.mtime = 0
            t.addfile(ti, io.BytesIO(data))
    return buf.getvalue()


# name -> (file, url, bytes served)
GOOD = {
    "fa": ("seqs.fasta", ZEN.format("seqs.fasta"), b">a\nACGT\n>b\nGGCC\n"),
    "gz": ("table.tsv.gz", ZEN.format("table.tsv.gz"), gz(b"x\ty\n1\t2\n")),
    "tar": ("bins.tar", FIG.format("424242"), tar({"a.fna": b">a\nAC\n"})),
    "txt": ("README, v1 (final).txt", ZEN.format("README,%20v1%20(final).txt"), b"read me\n"),
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
    for key, (_f, url, data) in GOOD.items():
        (blobs / key).write_bytes(data)
        served[url] = str(blobs / key)
    mp = tmp_path / "map.json"
    env = dict(os.environ)
    env.update({"PATH": f"{bindir}{os.pathsep}{env['PATH']}", "STUB_LOG": str(tmp_path / "calls.jsonl"),
                "STUB_MAP": str(mp), "TMPDIR": str(tmp_path)})
    w = {"tmp": tmp_path, "src": src, "dest": tmp_path / "incoming", "env": env, "map": mp, "served": served,
         "blobs": blobs}
    mp.write_text(json.dumps(served))
    return w


def serve(w: dict, url: str, spec: str | bytes) -> None:
    m = json.loads(w["map"].read_text())
    if isinstance(spec, bytes):
        b = w["blobs"] / f"extra{len(m)}"
        b.write_bytes(spec)
        spec = str(b)
    m[url] = spec
    w["map"].write_text(json.dumps(m))


def jobs(w: dict, rows: list[tuple], name: str = "jobs.tsv") -> Path:
    p = w["tmp"] / name
    p.write_text("".join("\t".join(r) + "\n" for r in rows))
    return p


def good_rows(stage: str = "S1", keys=("fa", "gz", "tar", "txt"), lic: str = "CC BY 4.0", with_md5: bool = True):
    rows = []
    for k in keys:
        f, url, data = GOOD[k]
        rows.append((stage, f, url, md5(data) if with_md5 else "-", lic))
    return rows


def run(w: dict, *lists: Path, args: tuple = (), env: dict | None = None) -> subprocess.CompletedProcess:
    cmd = ["bash", str(SCRIPT), "--dest", str(w["dest"]), "--sources", str(w["src"]), "--pause", "0", *args,
           *map(str, lists)]
    return subprocess.run(cmd, env=env or w["env"], capture_output=True, text=True, timeout=120)


def calls(w: dict) -> list[list[str]]:
    p = Path(w["env"]["STUB_LOG"])
    return [json.loads(line)["args"] for line in p.read_text().splitlines()] if p.exists() else []


def tree(base: Path) -> dict[str, str]:
    return {p.relative_to(base).as_posix(): sha(p.read_bytes()) for p in sorted(base.rglob("*")) if p.is_file()}


def test_a_good_list_lands_every_file_and_records_it(world):
    w = world
    r = run(w, jobs(w, good_rows(lic="CC BY 4.0") + good_rows(stage="batch5_x/Sub-2", keys=("fa",), lic="{}")))
    assert r.returncode == 0, r.stdout + r.stderr
    d = w["dest"] / "S1"
    for f, _u, data in GOOD.values():
        assert (d / f).read_bytes() == data
    assert (w["dest"] / "batch5_x" / "Sub-2" / "seqs.fasta").read_bytes() == GOOD["fa"][2]
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]
    assert sorted((d / "SHA256SUMS").read_text().splitlines()) == sorted(f"{sha(x[2])}  {x[0]}" for x in GOOD.values())
    rows = [line.split("\t") for line in (d / "SOURCE_URLS.tsv").read_text().splitlines()]
    assert rows[0] == ["filename", "url", "bytes", "sha256", "license", "date_utc"]
    by = {x[0]: x for x in rows[1:]}
    for f, u, data in GOOD.values():
        assert by[f][1:5] == [u, str(len(data)), sha(data), "CC BY 4.0"], by[f]
        assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", by[f][5])
    other = (w["dest"] / "batch5_x" / "Sub-2" / "SOURCE_URLS.tsv").read_text().splitlines()[1].split("\t")
    assert other[4] == "{}"  # a licence the list gives as {} is recorded as it stands
    assert r.stdout.count("OK\t") == 5 and "DONE\tok=5\tskip=0\tfail=0" in r.stdout
    assert not (w["dest"] / ".fetch_to_tower.lock").exists()


def test_every_request_is_a_plain_curl_with_no_identity(world):
    w = world
    rows = good_rows()
    r = run(w, jobs(w, rows))
    assert r.returncode == 0, r.stdout + r.stderr
    seen = calls(w)
    assert [c[-1] for c in seen] == [x[2] for x in rows]
    for c in seen:
        assert c[c.index("-A") + 1] == "curl"
        joined = " ".join(c)
        for bad in ("-H", "--header", "-u", "--user", "-e", "--referer", "--user-agent", "-b", "--cookie", "@"):
            assert bad not in c, (bad, joined)
        assert not re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+", joined), joined
        assert c[c.index("--proto") + 1] == "=https" and c[c.index("--proto-redir") + 1] == "=https"
        assert c[c.index("-o") + 1].startswith(".") and c[c.index("-o") + 1].endswith(".partial")


def test_a_second_run_skips_what_is_in_place(world):
    w = world
    lst = jobs(w, good_rows())
    assert run(w, lst).returncode == 0
    n, before = len(calls(w)), tree(w["dest"])
    r = run(w, lst)
    assert r.returncode == 0, r.stdout + r.stderr
    assert len(calls(w)) == n, "a file in place was fetched again"
    assert r.stdout.count("SKIP\t") == 4 and "DONE\tok=0\tskip=4\tfail=0" in r.stdout
    assert tree(w["dest"]) == before


def test_a_four_field_list_of_batch_one_is_accepted_and_the_licence_is_a_dash(world):
    w = world
    rows = [x[:4] for x in good_rows(keys=("fa",))]
    r = run(w, jobs(w, rows))
    assert r.returncode == 0, r.stdout + r.stderr
    line = (w["dest"] / "S1" / "SOURCE_URLS.tsv").read_text().splitlines()[1].split("\t")
    assert line[4] == "-"


def test_comments_and_blank_lines_are_skipped_and_several_lists_are_read(world):
    w = world
    a = w["tmp"] / "a.tsv"
    a.write_text("# stage\tfilename\turl\n\n" + "\t".join(good_rows(keys=("fa",))[0]) + "\n# end\n")
    b = jobs(w, good_rows(stage="S2", keys=("gz",)), name="b.tsv")
    r = run(w, a, b)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (w["dest"] / "S1" / "seqs.fasta").exists() and (w["dest"] / "S2" / "table.tsv.gz").exists()


UNCOVERED = [
    ZEN.format("x").replace("/111/", "/222/"),  # another record
    "https://zenodo.org/api/records/111/files/a/b/content",  # one segment became two
    "https://zenodo.org/api/records/111/files/x/content?download=1",  # a query the form does not have
    "https://zenodo.org/api/records/111/files/x/content/",
    "https://zenodo.org.evil.example/api/records/111/files/x/content",
    "https://zenodo.org@evil.example/api/records/111/files/x/content",
    "http://zenodo.org/api/records/111/files/x/content",  # not https
    "https://ndownloader.figshare.com/files/12/34",
    "https://ndownloader.figshare.com/files/",
    "https://example.org/literal/one.tar.gz.bak",
    "https://example.org/literal/two.tar.gz",
    "https://elsewhere.example/page",  # named in a prose `url:` value, which is not a form
    "https://not-a-form.example/x",  # inside folded text, not a key
    "https://legacy.example.org/dl/a/b/fasta/data/x.fa",  # {strain} is one segment
    "https://legacy.example.org/dl/a/gff/data/x.gff",  # only the fasta form is listed
    "https://stacks.example.org/file/druid:abc/",  # {path} needs one segment at least
    "https://stacks.example.org/file/druid:abc/x?y=1",
    "https://quoted.example/x/a/b",
    "https://github.com/anything/at/all",
]


@pytest.mark.parametrize("url", UNCOVERED)
def test_a_url_no_entry_covers_is_refused_before_any_request(world, url):
    w = world
    rows = good_rows() + [("S9", "evil.fa", url, "-", "x")]
    r = run(w, jobs(w, rows))
    assert r.returncode == 2, r.stdout + r.stderr
    assert ("not covered by any sources.yaml entry" in r.stderr or "not https" in r.stderr) and url in r.stderr
    assert calls(w) == [], "a request was made although a URL was refused"
    assert not w["dest"].exists()


COVERED = [
    ZEN.format("anything.fa"),
    FIG.format("7"),
    "https://example.org/literal/one.tar.gz",
    "https://legacy.example.org/dl/Strain-1/fasta/data/g.fasta",
    "https://stacks.example.org/file/druid:abc/a/b/c.tar.gz",
    "https://stacks.example.org/file/druid:abc/one.tsv",
    "https://quoted.example/x/f.txt",
]


@pytest.mark.parametrize("url", COVERED)
def test_a_url_of_an_entry_form_is_matched(world, url):
    w = world
    r = run(w, jobs(w, [("S1", "f.fa", url, "-", "x")]), args=("--check",))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f"match " in r.stdout and url in r.stdout
    assert calls(w) == [] and not w["dest"].exists()


def test_check_names_the_entry_and_fetches_nothing(world):
    w = world
    r = run(w, jobs(w, good_rows()), args=("--check",))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f"match zen {GOOD['fa'][1]}" in r.stdout and f"match fig {GOOD['tar'][1]}" in r.stdout
    assert "check: 4 URLs, all covered" in r.stdout
    assert calls(w) == [] and not w["dest"].exists()


def test_an_unusable_sources_file_stops_the_run(world):
    w = world
    w["src"].write_text("version: 1\nsources:\n  - name: none\n    url: not a url\n")
    r = run(w, jobs(w, good_rows()))
    assert r.returncode == 2 and calls(w) == []
    r = run(w, jobs(w, good_rows()), args=("--sources", str(w["tmp"] / "absent.yaml")))
    assert r.returncode == 2 and calls(w) == []


def test_a_form_with_a_variable_host_matches_nothing(world):
    w = world
    w["src"].write_text("version: 1\nsources:\n  - name: open\n    url: https://{host}/files/{file}\n")
    r = run(w, jobs(w, good_rows(keys=("fa",))))
    assert r.returncode == 2 and calls(w) == [], r.stdout + r.stderr


def test_a_md5_that_differs_removes_the_file_and_the_rest_still_land(world):
    w = world
    rows = good_rows()
    rows[1] = rows[1][:3] + (md5(b"another"),) + rows[1][4:]
    r = run(w, jobs(w, rows))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "FAIL_MD5" in r.stderr and rows[1][1] in r.stderr
    d = w["dest"] / "S1"
    assert not (d / rows[1][1]).exists()
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]
    for k in ("fa", "tar", "txt"):
        assert (d / GOOD[k][0]).read_bytes() == GOOD[k][2]
    assert not [x for x in (d / "SHA256SUMS").read_text().splitlines() if x.endswith(rows[1][1])]
    assert not [x for x in (d / "SOURCE_URLS.tsv").read_text().splitlines() if x.startswith(rows[1][1])]
    assert "DONE\tok=3\tskip=0\tfail=1" in r.stdout


@pytest.mark.parametrize("fault", ["curl", "empty", "gzip", "tar", "fasta"])
def test_a_failed_download_is_reported_removed_and_the_exit_is_one(world, fault):
    w = world
    f, url = {"curl": ("seqs.fasta", GOOD["fa"][1]), "empty": ("seqs.fasta", GOOD["fa"][1]),
              "gzip": ("table.tsv.gz", GOOD["gz"][1]), "tar": ("bins.tar", GOOD["tar"][1]),
              "fasta": ("seqs.fasta", GOOD["fa"][1])}[fault]
    serve(w, url, {"curl": "FAIL", "empty": "EMPTY", "gzip": b"not gzip", "tar": b"not a tar" * 100,
                   "fasta": b"no header\n"}[fault])
    r = run(w, jobs(w, good_rows(with_md5=False)))
    assert r.returncode == 1, r.stdout + r.stderr
    assert {"curl": "FAIL_CURL", "empty": "FAIL_EMPTY", "gzip": "FAIL_GZIP", "tar": "FAIL_TAR",
            "fasta": "FAIL_FASTA"}[fault] in r.stderr
    d = w["dest"] / "S1"
    assert not (d / f).exists()
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]
    assert sum(1 for x in (d / "SHA256SUMS").read_text().splitlines()) == 3


def test_a_file_in_place_that_is_not_in_sha256sums_is_fetched_not_trusted(world):
    w = world
    d = w["dest"] / "S1"
    d.mkdir(parents=True)
    (d / "seqs.fasta").write_bytes(b"left by a killed run, never recorded")
    (d / ".table.tsv.gz.partial").write_bytes(b"half")
    r = run(w, jobs(w, good_rows(keys=("fa", "gz"))))
    # mv -n never overwrites: the unrecorded file stays and the new copy is reported as a failed move
    assert r.returncode == 1, r.stdout + r.stderr
    assert (d / "seqs.fasta").read_bytes() == b"left by a killed run, never recorded"
    assert "FAIL_MV" in r.stderr
    assert (d / "table.tsv.gz").read_bytes() == GOOD["gz"][2]
    assert not [p for p in w["dest"].rglob("*") if p.name.endswith(".partial")]
    assert (d / "SHA256SUMS").read_text().splitlines() == [f"{sha(GOOD['gz'][2])}  table.tsv.gz"]


BAD_LINES = [
    ("S1", "a.fa", ZEN.format("a.fa")),  # three fields
    ("S1", "a.fa", ZEN.format("a.fa"), "-", "x", "extra"),
    ("../S1", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("a/../b", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("/abs", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("S 1", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("S;1", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("-S1", "a.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a/b.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "..", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a;rm.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a$(id).fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a`id`.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a'b.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "-rf.fa", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a.fa.partial", ZEN.format("a.fa"), "-", "x"),
    ("S1", "SHA256SUMS", ZEN.format("a.fa"), "-", "x"),
    ("S1", "SOURCE_URLS.tsv", ZEN.format("a.fa"), "-", "x"),
    ("S1", "a.fa", ZEN.format("a.fa") + " --output /tmp/x", "-", "x"),
    ("S1", "a.fa", ZEN.format("a.fa") + "'", "-", "x"),
    ("S1", "a.fa", ZEN.format("a.fa") + "`id`", "-", "x"),
    ("S1", "a.fa", "-o/tmp/x", "-", "x"),
    ("S1", "a.fa", "file:///etc/passwd", "-", "x"),
    ("S1", "a.fa", ZEN.format("a.fa"), "ABCDEF" * 5 + "AB", "x"),
    ("S1", "a.fa", ZEN.format("a.fa"), "c" * 31, "x"),
    ("S1", "a.fa", ZEN.format("a.fa"), "c" * 31 + ";", "x"),
    ("S1", "a.fa", ZEN.format("a.fa"), "-", "$(id)"),
    ("S1", "a.fa", ZEN.format("a.fa"), "-", "x;y"),
]


@pytest.mark.parametrize("row", BAD_LINES)
def test_a_list_line_with_anything_the_shell_could_read_is_refused_before_any_request(world, row):
    w = world
    p = w["tmp"] / "bad.tsv"
    p.write_text("\t".join(good_rows(keys=("fa",))[0]) + "\n" + "\t".join(row) + "\n")
    r = run(w, p)
    assert r.returncode == 2, r.stdout + r.stderr
    assert calls(w) == [] and not w["dest"].exists()


def test_the_same_file_twice_in_the_lists_is_refused(world):
    w = world
    rows = good_rows(keys=("fa",))
    r = run(w, jobs(w, rows + rows))
    assert r.returncode == 2 and "twice" in r.stderr and calls(w) == []


def test_a_destination_with_a_metacharacter_is_refused(world):
    w = world
    for bad in ("relative/dir", str(w["tmp"]) + "/a b", str(w["tmp"]) + "/x;y", str(w["tmp"]) + "/../y"):
        cmd = ["bash", str(SCRIPT), "--dest", bad, "--sources", str(w["src"]), str(jobs(w, good_rows()))]
        r = subprocess.run(cmd, env=w["env"], capture_output=True, text=True, timeout=60)
        assert r.returncode == 2, (bad, r.stdout + r.stderr)
    assert calls(w) == []


def test_a_stop_file_ends_the_run_cleanly_before_the_next_file(world):
    w = world
    w["dest"].mkdir()
    (w["dest"] / "STOP").write_text("")
    r = run(w, jobs(w, good_rows()))
    assert r.returncode == 0, r.stdout + r.stderr
    assert "STOP file" in r.stdout
    assert calls(w) == []
    assert not [p for p in w["dest"].rglob("*.fasta")]


def test_a_second_run_on_the_same_destination_is_refused_and_a_stale_lock_is_cleared(world):
    w = world
    lock = w["dest"] / ".fetch_to_tower.lock"
    lock.mkdir(parents=True)
    (lock / "pid").write_text(str(os.getpid()))  # this test process is alive
    r = run(w, jobs(w, good_rows()))
    assert r.returncode == 3, r.stdout + r.stderr
    assert calls(w) == [] and lock.exists()
    dead = subprocess.Popen(["true"])
    dead.wait()
    (lock / "pid").write_text(str(dead.pid))
    r = run(w, jobs(w, good_rows()))
    assert r.returncode == 0, r.stdout + r.stderr
    assert not lock.exists()


def test_the_script_calls_only_curl_and_never_reads_or_writes_outside_dest(world):
    w = world
    lst = jobs(w, good_rows())
    before = tree(w["tmp"])
    r = run(w, lst)
    assert r.returncode == 0, r.stdout + r.stderr
    after = tree(w["tmp"])
    new = sorted(set(after) - set(before))
    assert all(p.startswith("incoming/") or p == "calls.jsonl" for p in new), new
    assert not [p for p in after if p.startswith("fetch_to_tower.")], "a temporary directory was left behind"
    assert all(c[0] == "-fsSL" for c in calls(w))


def real_sources() -> Path:
    """The tree's sources.yaml; the test image (built from main, container/Dockerfile) does not ship it, so skip there."""
    p = ROOT.joinpath("sources.yaml")
    if not p.is_file():
        pytest.skip("sources.yaml is not in the image: this check runs from a checkout (ci and the pod's mirror)")
    return p


REAL_URLS = """
https://ndownloader.figshare.com/files/66974843
https://github.com/thackl/pro-tycheposons/archive/a0f299fd6932e5ea81c52cbaf6f2c81aa19be724.tar.gz
https://zenodo.org/api/records/16656823/files/Inocle_001.fasta/content
https://ndownloader.figshare.com/files/54652007
https://g-83fdd0.1beed.03c0.data.globus.org/cartuids/35846_b8b9149f2f59ff950549b5c89b8b35e2.tar
https://ndownloader.figshare.com/files/60829597
https://zenodo.org/api/records/3247846/files/mmetsp_dib_trinity2.2.0_contigIDnamemaps.tar.gz/content
https://zenodo.org/api/records/22794128/files/Donor%20strain-pKJK5%20AtraI.fasta/content
https://ndownloader.figshare.com/files/57603328
https://zenodo.org/api/records/13143815/files/metagenome_plasmidome_plasmid2.fasta/content
https://ndownloader.figshare.com/files/43971075
https://zenodo.org/api/records/14540677/files/README.txt/content
https://ndownloader.figshare.com/files/65710359
https://portal.nersc.gov/cfs/nelli/ms_protistsymb/viruses.fna.tar.gz
https://zenodo.org/api/records/11215303/files/Metadata%20for%20DPANN%20viruses%20and%20unclassified%20MGEs%20assembled%20from%20metagenomes.txt/content
https://zenodo.org/api/records/20174424/files/Nanobdellota_Zen.zip/content
https://zenodo.org/api/records/19235833/files/Metadata%20-%20Metagenomes%20from%20the%20Red%20Sea%20hydrothermal%20vents.xlsx/content
https://zenodo.org/api/records/19510112/files/Soil_activate_prophage.fasta.gz/content
https://zenodo.org/api/records/21773290/files/vOTU-abundance-profile.txt/content
https://ndownloader.figshare.com/files/51500924
https://ndownloader.figshare.com/files/54602522
https://zenodo.org/api/records/7412085/files/RVD_taxonomy_classification.xlsx/content
https://ndownloader.figshare.com/files/71721769
https://zenodo.org/api/records/23322224/files/mobilome-recovery-workflow.zip/content
https://ndownloader.figshare.com/files/61896136
https://zenodo.org/api/records/13285428/files/cFMD_humann.tar.gz/content
https://ndownloader.figshare.com/files/62219453
https://ndownloader.figshare.com/files/55385462
https://zenodo.org/api/records/15856862/files/HAGs_vHAGs_madorio-Jul25%20freshwater.fasta/content
https://ndownloader.figshare.com/files/46078842
https://zenodo.org/api/records/21563086/files/all.MAGs.tar.gz/content
https://ndownloader.figshare.com/files/41695503
https://zenodo.org/api/records/7847962/files/mOTU_table.xlsx/content
https://zenodo.org/api/records/20276352/files/RdRp_CAB_v1.faa/content
https://ndownloader.figshare.com/files/51521069
https://zenodo.org/api/records/14537680/files/MFD-LR-drep.tar.gz/content
https://zenodo.org/api/records/19066579/files/midas_core_genera_hqmag_summary.tsv/content
https://ndownloader.figshare.com/files/69455682
https://zenodo.org/api/records/23755631/files/141_SR_MAGs.tar.gz/content
https://ndownloader.figshare.com/files/40534991
https://zenodo.org/api/records/17533381/files/starship-scan.py/content
https://ndownloader.figshare.com/files/47671063
https://zenodo.org/api/records/7368133/files/RVMT_Zenodo_V4.zip/content
https://zenodo.org/api/records/5202946/files/Benjamin-Lee/viroiddb-2021-06-06.zip/content
https://zenodo.org/api/records/6859104/files/analysis.tar.gz/content
https://zenodo.org/api/records/21287531/files/sandpiper2.0.1.parsed_metadata.tsv.gz/content
https://zenodo.org/api/records/20503966/files/prokaryote_bins_5283.tar.bz2/content
https://ndownloader.figshare.com/files/68890756
https://zenodo.org/api/records/16263068/files/Archaeal-viruses.fasta/content
https://stacks.stanford.edu/file/druid:wb363nt3637/Obelisk_MS_Oct2024/preprint/Clustered_Obelisks.tar.gz
https://zenodo.org/api/records/17823213/files/globalRNAvirus.zip/content
https://ndownloader.figshare.com/files/58496362
https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/227/S-BSST2227/Files/DSV70_archaea_gtdbtk.unrooted.tree
https://ont-open-data.s3.amazonaws.com/compost_mgx_2026.04/analysis/checkm2.quality_report.tsv
https://legacy.tritrypdb.org/common/downloads/Current_Release/AdeanaiCavalhoATCCPRA-265/fasta/data/TriTrypDB-68_AdeanaiCavalhoATCCPRA-265_Genome.fasta
https://legacy.tritrypdb.org/common/downloads/Current_Release/AdeanaiCavalhoATCCPRA-265/gff/data/TriTrypDB-68_AdeanaiCavalhoATCCPRA-265.gff
https://legacy.microsporidiadb.org/common/downloads/Current_Release/AalgeraePRA109/fasta/data/MicrosporidiaDB-68_AalgeraePRA109_Genome.fasta
https://legacy.microsporidiadb.org/common/downloads/Current_Release/AalgeraePRA109/gff/data/MicrosporidiaDB-68_AalgeraePRA109.gff
https://legacy.giardiadb.org/common/downloads/Current_Release/GintestinalisAssemblageAAS175/fasta/data/GiardiaDB-68_GintestinalisAssemblageAAS175_Genome.fasta
https://legacy.giardiadb.org/common/downloads/Current_Release/GintestinalisAssemblageADH/gff/data/GiardiaDB-68_GintestinalisAssemblageADH.gff
https://zenodo.org/api/records/15788452/files/02.BGC_final.tar.xz/content
https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs/OMDBv2.0_data.tsv.gz
https://ndownloader.figshare.com/files/44053844
https://ndownloader.figshare.com/files/71316448
https://ndownloader.figshare.com/files/47726602
https://ndownloader.figshare.com/files/56190881
https://zenodo.org/api/records/22834836/files/GKS_genomes.tar.gz/content
https://zenodo.org/api/records/17735769/files/README.md/content
https://zenodo.org/api/records/16414468/files/Supplementary_Dataset.xlsx/content
https://zenodo.org/api/records/22159188/files/high_quality_vOTUs_final.zip/content
https://zenodo.org/api/records/21450953/files/GHV.fasta/content
https://zenodo.org/api/records/14684892/files/CGVD.gz/content
https://ndownloader.figshare.com/files/43845801
https://zenodo.org/api/records/23266178/files/Supplementary_Data_7.xlsx/content
https://zenodo.org/api/records/17143670/files/OPD-host.xlsx/content
https://data.mendeley.com/public-files/datasets/wdyjgcdgbm/files/70495243-7b59-45f1-846c-5641b75512d6/file_downloaded
https://zenodo.org/api/records/20399072/files/symgs_catalog_release1_fna.tar.gz/content
https://zenodo.org/api/records/16757580/files/readme/content
https://zenodo.org/api/records/16784684/files/earthhologenome/EHI_data_release_1-1.0.3.zip/content
https://zenodo.org/api/records/17468704/files/metadata_MAGs_v2_EMERGE.tsv/content
https://ndownloader.figshare.com/files/69639846
https://zenodo.org/api/records/12617226/files/geNomad-output.zip/content
https://zenodo.org/api/records/15743210/files/2060vMAGs.zip/content
https://zenodo.org/api/records/7313591/files/MAGs_genome_info.tsv/content
https://ndownloader.figshare.com/files/57783499
https://zenodo.org/api/records/21341769/files/merged_MAGs.zip/content
https://zenodo.org/api/records/13643296/files/MAGs_catalogue_v1.0.0.xlsx/content
https://ndownloader.figshare.com/files/67605795
https://zenodo.org/api/records/13266492/files/paobrien/SpongeSymbionts-SpongeSymbionts.zip/content
https://zenodo.org/api/records/21135542/files/Phyllospongia%20foliascens_metagenome_bins.csv/content
https://zenodo.org/api/records/22053915/files/WT29A_MQ_MAGs.zip/content
https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/marine/v3.0/genomes-all_metadata.tsv
""".split()


def test_every_source_of_the_2026_11_25_job_lists_is_covered_by_the_trees_sources_yaml(world):
    """One real URL per source and per URL form of `_batch1_jobs.tsv` .. `_batch5.tsv` on tower (read 2026-11-25)."""
    w = world
    assert len(REAL_URLS) == len(set(REAL_URLS)) == 92
    rows = [("S1", f"f{i}.bin", u, "-", "x") for i, u in enumerate(REAL_URLS)]
    r = subprocess.run(["bash", str(SCRIPT), "--check", "--sources", str(real_sources()), str(jobs(w, rows))],
                       env=w["env"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stdout + r.stderr
    assert f"check: {len(REAL_URLS)} URLs, all covered" in r.stdout
    assert calls(w) == []


def test_a_url_nobody_listed_is_refused_with_the_trees_sources_yaml(world):
    w = world
    rows = [("S1", "f.bin", "https://example.org/some/file.fa", "-", "x")]
    r = subprocess.run(["bash", str(SCRIPT), "--check", "--sources", str(real_sources()), str(jobs(w, rows))],
                       env=w["env"], capture_output=True, text=True, timeout=120)
    assert r.returncode == 2 and "not covered" in r.stderr


def test_the_tree_sources_yaml_parses_and_every_tower_entry_has_its_record():
    yaml = pytest.importorskip("yaml")
    doc = yaml.safe_load(real_sources().read_text())
    assert doc["version"] == 1
    names = [e["name"] for e in doc["sources"]]
    assert len(names) == len(set(names))
    tower = [e for e in doc["sources"] if e["name"].startswith("tower_")]
    assert len(tower) >= 70
    for e in tower:
        assert e["url"].startswith("https://") and "license" in e, e["name"]
        assert str(e["approved"]).startswith("by the PI, 2026-11-25 "), e["name"]
        assert "fetched" in e and "stored" in e, e["name"]
        for f in e.get("files", []):
            assert re.fullmatch(r"[0-9a-f]{64}", f["sha256"]) and isinstance(f["bytes"], int) and f["bytes"] > 0, (e["name"], f)
            assert "md5" not in f or re.fullmatch(r"[0-9a-f]{32}", f["md5"]), (e["name"], f)
    gsv = next(e for e in doc["sources"] if e["name"] == "gsv_atlas")
    assert "8ff20350e9986f2a99efccf4ae7fb5ccd7f18d8899635a4e816fb6b2ad987df0" in gsv["fetched_tower"]
    gsv_tower = next(e for e in tower if e["name"] == "tower_gsv1")
    assert gsv_tower["files"][0]["sha256"] == "8ff20350e9986f2a99efccf4ae7fb5ccd7f18d8899635a4e816fb6b2ad987df0"
    assert next(e for e in tower if e["name"] == "tower_gbr_mgd")["license"] == "{}"
