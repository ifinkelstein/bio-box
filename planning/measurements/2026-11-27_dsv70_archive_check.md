# Archive diagnosis for `DSV70_mag.tar.gz` (`tower_dsv70`, card R-24)

Measured `2026-11-27` (`19:22Z`–`19:38Z`) on `tower` (`root@tower`) by read-only inspection of
`/mnt/user/research-bulk/bio-box/sources_incoming/DSV70/`, `/mnt/user/research-bulk/bio-box/fetch_logs/`, and HTTP `HEAD`
plus byte-range (`Range: bytes=...`) requests against the exact URL already listed under `tower_dsv70` in
[`sources.yaml`](../../sources.yaml):
`https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/227/S-BSST2227/Files/DSV70_mag.tar.gz`.

All probe slices were written only under `/tmp/r24_dsv70/` on `tower` (`-A curl --proto =https --proto-redir =https`, no
identities or parameters sent, `AGENTS.md` §8) and removed (`rm -r /tmp/r24_dsv70`, verified `test ! -e /tmp/r24_dsv70`
exit `0`) before finishing. Nothing under `/mnt/user/research-bulk/bio-box/` was written, modified, or deleted, and no
script or job list was changed.

## 1. Summary diagnosis

- **Verdict (card R-24 case `(b)`/`(c)` — server-side mid-stream corruption inside a genuine `.tar.gz`, not `(a)` a
  mislabelled format and not tail truncation).** EBI BioStudies FIRE (`S-BSST2227`) serves a `4,435,296,005`-byte
  (`4.4353 GB`, *derived*) HTTP/1.1 object (`Content-Type: application/x-gzip`, no `Content-Encoding`, `Accept-Ranges:
  bytes`, `ETag: "1085d3f05-65071af6dd040"`, `Last-Modified: Mon, 27 Apr 2026 14:11:05 GMT`; *measured*, `tower`,
  `2026-11-27`, `curl -sSIL`).
- **Head (`bytes 0-65535`) is valid `gzip` + POSIX `ustar` (`(a)` ruled out).** Bytes `0..9` are `1f 8b 08 00 00 00 00 00
  00 03` (RFC 1952 `gzip`, deflate, Unix), identical to the two sibling `.tar.gz` archives of `S-BSST2227` that already
  landed on `tower` (*measured*, `od -tx1`). Decompressing the first `65,536` bytes (`tar -ztvf`, `tar -zxf`) yields
  `DSV70_mag/` followed by `DSV70_mag/444-23-R2-1-1360_bin.152.fa` (`1,997,170` bytes in the tar header; `216,064` bytes
  of clean uppercase `ACGT` FASTA extracted from the `65,536`-byte slice, starting `>lcl_k127_1399159`; *measured*).
- **Tail (`bytes 4435291909-4435296004`) ends with a complete `gzip` trailer (tail truncation ruled out).** The final
  `64` bytes hold deflate zero-run blocks (`af df bf 7e ff fa fd eb f7` repeated, compressing tar's trailing NUL blocks),
  deflate end-of-block (`ff 03`), and an `8`-byte RFC 1952 trailer `bf dc 40 74 00 48 41 69` with `0` trailing garbage
  bytes (*measured*, `od -tx1 -j 4032 -N 64`). Its little-endian `ISIZE = 0x69414800` (`1,765,885,952` bytes mod `2^32`,
  *derived*) corresponds at `3 × 2^32 + ISIZE` to `14,650,787,840` uncompressed bytes (`1,430,741 × 10,240`-byte GNU
  `tar` records, `mod 10,240 = 0`; *derived*).
- **Why `chain_r22.sh` failed `FAIL_GZIP`.** In `2026-11-27_chain_r22.log` (`12:31:39Z` to `14:02:45Z`, `5,466 s` =
  `91 min 06 s`, *derived*; `811.4 KB/s` average over `4,435,296,005` bytes, *derived*), `curl 8.16.0` (`-fsSL --retry 3`)
  wrote `0` `curl:` error lines to stderr and exited `0`, whereas every failed or retried `curl` attempt in the same log
  (`lines 1354..3625`) and in `_batch4_detached.log` (`line 23`) prints `curl: (<errno>) ...` to stderr under `-S`
  (*measured*, `grep -n -C 2 "curl:"`). Thus `curl` received all `4,435,296,005` TLS-authenticated HTTP/1.1 bytes in one
  clean pass, and `gzip 1.14 -t` rejected the complete server object (`FAIL_GZIP`) because of internal deflate/CRC32
  corruption mid-stream on EBI's server.
- **Recommendation (do not retry `18_later_rescue.tsv` as-is).** Because `18_later_rescue.tsv` has `md5 = -`,
  `fetch_to_tower.sh` requires `gzip -t` to pass and unlinks `.DSV70_mag.tar.gz.partial` on `FAIL_GZIP`; re-running it
  will spend ~1.5 h (*derived* from `5,466 s`) and delete the file again. Instead, a future card should fetch the `7,422`
  MAGs (*estimate*, job-list/study count) by accession from NCBI BioProject `PRJNA1244896` (individual `.fna.gz` files
  with NCBI `md5checksums.txt`), using the two GTDB-Tk trees (`DSV70_archaea_gtdbtk.unrooted.tree`,
  `DSV70_bacteria_gtdbtk.unrooted.tree`) and `DSV70_assemblies.tar.gz` (`8,564,570,028` bytes, `gzip -t` clean) already on
  `tower`.

---

## 2. Current state of `sources_incoming/DSV70/` and fetch logs on `tower`

Read-only check of `/mnt/user/research-bulk/bio-box/sources_incoming/DSV70/` (`ls -la`, `cat SOURCE_URLS.tsv`, `cat
SHA256SUMS`; *measured*, `tower`, `2026-11-27`):

| File | Bytes (*measured*) | SHA-256 (*measured*, `SHA256SUMS`) | Landed UTC (*measured*, `SOURCE_URLS.tsv`) | Contents / notes (*measured*, `tower`, `2026-11-27`) |
|---|---:|---|---|---|
| `DSV70_archaea_gtdbtk.unrooted.tree` | `67,354` | `944bd7810764110414ec3ea577e4cb30fda6bfa383b9f351d39a8f947c3ec3b1` | `2026-11-26T08:13:55Z` | GTDB-Tk archaeal tree |
| `DSV70_bacteria_gtdbtk.unrooted.tree` | `300,969` | `50623bb43b756bbf489d7a16d86d8f7ebf298d3721faecd4a3f8186a2fd772f5` | `2026-11-26T08:14:12Z` | GTDB-Tk bacterial tree |
| `DSV70_ntus.tar.gz` | `1,707,597` | `03d8c015049d68a4054bdba6674ab5169704a0672375b62a8495a09b0dbce9c6` | `2026-11-26T08:15:21Z` | `1` dir (`DSV70_ntus/`) + `70` `*.phyloFlash.NTUfull_abundance.csv` files (`tar -ztvf`) |
| `DSV70_assemblies.tar.gz` | `8,564,570,028` | `ef78630187e8fb50501ea52688d2119bd25b9e0f4835ce2e5dd8ccaa1555b80d` | `2026-11-26T11:37:37Z` | `gzip` (`1f 8b 08 00 00 00 00 00 00 03`), starts `DSV70_assemblies/363-7-R1-3-628_assembly.fasta` (`159,607,862` B) |

- `4` of `5` files are present; `0` `.partial` files remain; directory mtime is `Nov 27 08:02` CST (`14:02 UTC`,
  *measured*, `ls -la`), when `fetch_to_tower.sh` removed `.DSV70_mag.tar.gz.partial`. Free disk on `/mnt/user` is `16 TiB`
  of `32 TiB` (`52%` used; *measured*, `tower`, `2026-11-27`, `df -h /mnt/user/research-bulk`).

All three recorded download attempts for `DSV70_mag.tar.gz` (*measured*, `tower`, `2026-11-27`, log inspection):

1. **Batch 4 (`_batch4_detached.log` lines 24–28, `2026-11-26`):** between `08:15:21Z` (`DSV70_ntus.tar.gz`) and
   `11:37:37Z` (`DSV70_assemblies.tar.gz`), `curl` logged `curl: (28) Connection timed out after 30002 milliseconds` and
   `curl: (56) OpenSSL SSL_read: OpenSSL/3.5.4: error:0A000126:SSL routines::unexpected eof while reading, errno 0`, then
   `FAIL_CURL`. Note that on line 23 (`DSV70_ntus.tar.gz`), a transient `curl: (28)` retry that succeeded still left
   `curl: (28) Connection timed out after 30002 milliseconds` in stderr before `OK`.
2. **Lists 01–20 (`2026-11-26_lists_01_to_20.log` lines 11437–11440, `2026-11-26T23:53:02Z` to `2026-11-27T01:04:15Z`):**
   `curl: (56) OpenSSL SSL_read: OpenSSL/3.5.4: error:0A000126:SSL routines::unexpected eof while reading, errno 0`, then
   `FAIL_CURL`. (In the same log at line 12156, a premature clean HTTP close on `SMAGs_contigs_individual.fna.tar.gz`
   logged `curl: (18) end of response with 6213389723 bytes missing` and `FAIL_CURL`.)
3. **`chain_r22.sh` (`2026-11-27_chain_r22.log` lines 1–5, `2026-11-27`):**
   ```text
   2026-11-27T12:31:39Z start 18_later_rescue.tsv (jobs 1, pause 1)
   FAIL_GZIP	DSV70	DSV70_mag.tar.gz
   SKIP	EUK_Y1000	y1000p_genome_files.zip
   DONE	ok=0	skip=1	fail=1
   2026-11-27T14:02:45Z end 18_later_rescue.tsv exit 0
   ```
   Elapsed wall time is `5,466 s` (`1 h 31 m 06 s`, *derived* from `14:02:45Z - 12:31:39Z`), giving an effective rate of
   `811.4 KB/s` (`4,435,296,005 / 5,466`, *derived*, including `gzip -t` scan time), matching EBI FIRE's single-stream
   rate on `DSV70_assemblies.tar.gz` (`>= 705.7 KB/s`, *derived* from `8,564,570,028 B / 12,136 s`). No `curl:` message
   appears anywhere before line 1354 (`14:10:36Z`, inside `22_later_batch7_EEMC1.tsv`, where `curl 8.16.0` under `-fsSL
   --retry 3` logs `curl: (28)` and `curl: (19)` on stderr whenever a retry fires; *measured*, `grep -n -C 2 "curl:"`).

---

## 3. HTTP response headers

Command on `tower` (`2026-11-27T19:24:13Z`):

```sh
curl -sSIL -A curl --proto =https --proto-redir =https \
  -D /tmp/r24_dsv70/headers.txt -o /dev/null \
  https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/227/S-BSST2227/Files/DSV70_mag.tar.gz
```

Measured response (*measured*, `tower`, `2026-11-27`):

```http
HTTP/1.1 200 OK
Date: Fri, 27 Nov 2026 19:24:13 GMT
Server: Apache
Last-Modified: Mon, 27 Apr 2026 14:11:05 GMT
ETag: "1085d3f05-65071af6dd040"
Accept-Ranges: bytes
Content-Length: 4435296005
Access-Control-Allow-Origin: *
Access-Control-Allow-Methods: GET,OPTIONS
Access-Control-Allow-Headers: Authorization,Origin,User-Agent,X-Requested-With,If-Modified-Since,Cache-Control,Content-Type,Range,Accept
Content-Type: application/x-gzip
```

- `Content-Length`: `4,435,296,005` bytes (*measured*) = `4.435296005 GB` decimal / `4.130691 GiB` binary (*derived*),
  within `-0.11%` (*derived*) of the `4.44 GB` *estimate* in `18_later_rescue.tsv`.
- `ETag`: `"1085d3f05-65071af6dd040"`, where `0x1085d3f05 = 4,435,296,005` and `0x65071af6dd040 = 1,777,299,065,000,000`
  µs since Unix epoch (`2026-04-27T14:11:05Z`, *derived*), matching `Last-Modified`.
- `Content-Encoding`: absent (*measured*); Apache serves the raw stored object without on-the-fly `Content-Encoding: gzip`.
- `Accept-Ranges: bytes`: honored with `HTTP/1.1 206 Partial Content` on both `Content-Range: bytes 0-4095/4435296005` and
  `Content-Range: bytes 4435291909-4435296004/4435296005` (*measured*, `tower`, `2026-11-27`).

---

## 4. Head and tail byte-range probes

| Slice | Byte range (*measured*) | Length (B, *measured*) | SHA-256 (*measured*, `sha256sum`) |
|---|---|---:|---|
| `head.bin` | `0–4,095` | `4,096` | `762b09aa9ca78c6a067d033acc46bd91006a1d1a4098890587d08a6e394532c1` |
| `head64k.bin` | `0–65,535` | `65,536` | `ed5e6c9d56f1d8a3a03b66af0f848841179c8f411ed728f709d742597e4f303f` |
| `tail.bin` | `4,435,291,909–4,435,296,004` | `4,096` | `c49efa216639f851a663fa6c1ee77c0bbbf59cb1ec2797c2ceed9cad806d520e` |

### 4.1 Head (`bytes 0-4095` and `0-65535`)

First `64` bytes (`od -tx1 -N 64 /tmp/r24_dsv70/head.bin`, *measured*, `tower`, `2026-11-27`):

```text
0000000 1f 8b 08 00 00 00 00 00 00 03 e4 bd cf ae 26 c9
0000020 8e e4 d7 6b 7f 8a 7e 81 1e 65 66 e5 df 8d 00 82
0000040 03 70 2f 11 da 36 ee 08 90 20 a8 35 8b 19 e9 fd
0000060 95 41 fb 19 3d 4e dd ea 9e 01 a4 59 08 aa be 5d
```

- Bytes `0..9` (`1f 8b 08 00 00 00 00 00 00 03`): standard RFC 1952 `gzip` header (`ID1=0x1f`, `ID2=0x8b`, `CM=8` deflate,
  `FLG=0x00`, `MTIME=0`, `XFL=0x00`, `OS=0x03` Unix), matching `DSV70_ntus.tar.gz` and `DSV70_assemblies.tar.gz` byte-for-byte
  across the `10`-byte header (*measured*, `od -tx1 -N 16`). Byte `10` (`0xe4`) begins a non-final dynamic-Huffman deflate
  block (`BFINAL=0`, `BTYPE=2`, *derived*).
- Bytes `256..287` of the outer file (`od -tx1 -tc -j 256 -N 32`, *measured*) are compressed deflate data (`cf ff cf d7 ...`),
  not `ustar` at offset `257`. Thus `DSV70_mag.tar.gz` is **not** an uncompressed `.tar` (unlike `NX2_CAVE`), nor `bzip2`,
  `xz`, `zstd`, or an HTML error document.
- `tar -ztvf /tmp/r24_dsv70/head64k.bin` and `tar -zxf /tmp/r24_dsv70/head64k.bin` (*measured*, `tower`, `2026-11-27`)
  decompress the `65,536`-byte prefix cleanly into:
  - `drwxrwxrwx alrlab/alrlab 0 2025-09-15 18:23 DSV70_mag/`
  - `-rwxrwxrwx alrlab/alrlab 1997170 2025-07-03 17:59 DSV70_mag/444-23-R2-1-1360_bin.152.fa` (`216,064` bytes extracted
    before the `65,536`-byte slice ends, starting `>lcl_k127_1399159` and `70`-bp uppercase `ACGT` lines; *measured*,
    `head -n 6` and `wc -c`).

### 4.2 Tail (`bytes 4435291909-4435296004`)

Last `64` bytes of `tail.bin` (`od -tx1 -j 4032 -N 64 /tmp/r24_dsv70/tail.bin`, *measured*, `tower`, `2026-11-27`):

```text
0007700 af df bf 7e ff fa fd eb f7 af df bf 7e ff fa fd
0007720 eb f7 af df bf 7e ff fa fd eb f7 af df bf 7e ff
0007740 fa fd eb f7 af df bf 7e ff fa fd eb f7 af df bf
0007760 7e ff fa ff f7 af ff 03 bf dc 40 74 00 48 41 69
```

- Preceding the final `8` bytes is a repeating `9`-byte deflate zero-run (`af df bf 7e ff fa fd eb f7`, encoding the
  trailing NUL padding of the POSIX `.tar`) closed by deflate end-of-block (`ff 03`), followed at bytes
  `4,435,295,997..4,435,296,004` by the `8`-byte RFC 1952 trailer `bf dc 40 74 00 48 41 69` with `0` trailing bytes after
  it (*measured*).
- Little-endian trailer values (*derived*):
  - `CRC32 = 0x7440dcbf` (`1,950,407,871`)
  - `ISIZE = 0x69414800` (`1,765,885,952` bytes modulo `2^32`)
- Uncompressed archive size check against GNU `tar`'s default `10,240`-byte (`20 × 512` B) record size (*derived* from
  `k × 2^32 + 1,765,885,952`):
  - `k = 0`: `1,765,885,952` B (`mod 10,240 = 8,192`)
  - `k = 1`: `6,060,853,248` B (`mod 10,240 = 2,048`)
  - `k = 2`: `10,355,820,544` B (`mod 10,240 = 6,144`)
  - `k = 3`: `14,650,787,840` B = `1,430,741 × 10,240` B (`mod 10,240 = 0`, `14.6508 GB`; at `7,422` MAGs *estimate*,
    `1.974 MB` per MAG *derived*, consistent with the first member's `1,997,170` B; sibling `DSV70_assemblies.tar.gz` has
    trailer `b5 02 30 24 00 d0 8b 8e`, `ISIZE = 0x8e8bd000 = 2,391,527,424`, which at `k = 6` is `28,161,331,200` B =
    `2,750,130 × 10,240` B, *measured* and *derived*).

Because both the `gzip` header at byte `0` and the `gzip` trailer at byte `4,435,296,004` are intact on the server while a
complete, zero-error `curl 8.16.0` download (`gzip 1.14 -t` on `tower`) fails `FAIL_GZIP`, the corruption resides inside
the compressed stream stored on EBI BioStudies FIRE (a corrupt deflate block or internal upload-chunk corruption causing a
deflate error or CRC32/ISIZE mismatch).

---

## 5. Recommended next step

1. **Do not retry `DSV70_mag.tar.gz` through `fetch_to_tower.sh` (`18_later_rescue.tsv`).** With `md5 = -`, `fetch_to_tower.sh`
   enforces `gzip -t` and unlinks `.DSV70_mag.tar.gz.partial`, wasting ~`1.5 h` (`5,466 s`, *derived*) per attempt.
2. **Preferred route for a follow-up staging/fetch card:** pull the `7,422` MAGs (*estimate*, study/job-list count) by
   INSDC assembly accession under NCBI BioProject `PRJNA1244896` (adding the route to `sources.yaml` first per `AGENTS.md`
   §8), where each genome is an independent `.fna.gz` with a provider MD5 in `md5checksums.txt`. All four companion files
   of `S-BSST2227` — the two GTDB-Tk trees (`67,354` B and `300,969` B), `DSV70_ntus.tar.gz` (`1,707,597` B), and the full
   metagenome assemblies `DSV70_assemblies.tar.gz` (`8,564,570,028` B) — are already checksummed and intact in
   `tower:/mnt/user/research-bulk/bio-box/sources_incoming/DSV70/` (*measured*).
3. **Optional salvage check if any `DSV70` bin is absent from `PRJNA1244896`:** download `DSV70_mag.tar.gz` once to a
   scratch path outside `fetch_to_tower.sh` so the file is retained, then run `gzip -tv` and `tar -ztvf` to locate the
   exact corrupt byte/member offset and count recoverable `DSV70_mag/*.fa` entries.
