# A-28: External Bare-Metal & Metal-Cloud Compute for `bio-box` (2025–2026)

*Prepared: 2026-11-25 | Card: `A-28` (**H** — PI decision + operational bootstrap/sync checklist) | Scope: Read-only market and engineering evaluation of commercial bare-metal compute and UT Austin TACC allocations to parallelize `bio-box` stage builds (Profile A) and big-release clustering (Profile B) alongside `wilkcomp02` (112 threads / 1.0 TiB RAM / 1.1 TiB `/tmp` NVMe; shared 40-core day / 80-core night cap, `DEC-73`, `DEC-74`).*

*Number labeling convention (`AGENTS.md` §6): every price is labeled **published** (provider 2025–2026 list pricing) or **derived** (formula and inputs stated).*

---

## Executive Summary & Top Recommendations

Because every `bio-box` step runs inside the pinned x86_64 Singularity/Apptainer `.sif` image (`bio-box_faec1f4.sif`, verified by `sha256: d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`) and produces **byte-identical Parquet outputs on any x86_64 Linux host** (`AGENTS.md` §4.3), external bare metal is a drop-in horizontal multiplier requiring **zero core-code changes** — only a dedicated `profiles/<host>.yaml` and a read-only sync boundary.

**Three structural 2025–2026 market facts shape the decision:**
1. **2025–2026 DDR5 DRAM & enterprise NVMe price surge:** Monthly bare-metal providers (Hetzner, OVHcloud) significantly raised prices on 512 GB–1.5 TB DDR5 SKUs and charge steep month-to-month setup fees (Hetzner `AX162-3` 512 GB is `€1,597/mo + €799 setup` *published*; OVHcloud charges `1×` monthly rent setup on month-to-month).
2. **Equinix Metal is permanently defunct** (EOL **June 30, 2026**, *published*).
3. **Burst vs. Sustained split:**
   - For **bursty big-release clustering** (`24–48 h` per release, **Profile B**), **Oracle Cloud (OCI) Bare Metal** (`BM.Standard.E6.256` at **`$13.824/hr`** *published* for **512 Turin threads + 3 TB RAM**, or `BM.DenseIO.E5.128` at **`$15.60/hr`** *published* for **256 Genoa threads + 1.5 TB RAM + 81.6 TB local NVMe**) with **10 TB/month FREE egress** beats every monthly bare-metal provider on total cash outlay (`~$346–$750` per full release *derived*, `$0` idle burn, `$0` setup).
   - For **sustained multi-week stage building** (**Profile A**), **Cherry Servers** (US Chicago, `$0` setup, hourly or monthly, `100 TB` free egress) and **Hetzner Server Auction** (EU, `€0` setup, unlimited egress) crush hyperscalers on `$/CPU-hour`.

### Top 3 Winners — Profile A (Stage-Build Workhorse: 64–192t, 256–512 GB RAM, 4–8 TB NVMe)
1. **Cherry Servers — AMD EPYC 9554P (US Chicago / EU):** `64c/128t`, `384 GB` DDR5, NVMe, **`100 TB` free egress**, **`$0` setup**, **`$1.34/hr`** or **`$630–$783/mo`** (*published*; **`$8.03–$10.47` / 1,000 thread-hr**, *derived*). Best US balance of price, zero setup fee, and hourly-to-monthly flexibility.
2. **Hetzner Server Auction / AX162-R (DE/FI — EU):** `48–64c` (`96–128t`), `256–512 GB` RAM, `4+ TB` NVMe, **unlimited egress (`$0/TB`)**, **`€0` setup on Auction** (`~$165–$495/mo` *published*, **`$4.71–$7.06` / 1,000 thread-hr**, *derived*). Unbeatable raw `$/CPU-hour` when EU latency (`~120 ms` to Austin, immaterial for batch `rsync` of `20–50 GB` stage tars) is acceptable. *(Hetzner US Ashburn/Hillsboro is Cloud VM only — no bare metal.)*
3. **Oracle Cloud (OCI) — `VM.Standard.E5.Flex` (US Dallas / Ashburn / Phoenix):** `64 OCPU / 128 vCPU` (Genoa), `512 GB` RAM + `4 TB` Block Volume, **`$0` setup**, **`10 TB/mo` free egress**, **`~$3.18/hr`** (*derived* from `$2.944/hr` compute + `$0.24/hr` block *published*; **`~$76` per 24-hour stage build**, *derived*). Best for pure on-demand bursts with `$0` idle cost.

### Top 3 Winners — Profile B (Big-Release Clustering Monster: 128–256+t, 768 GB–3 TB RAM, 4–80+ TB NVMe)
1. **Oracle Cloud (OCI) — `BM.Standard.E6.256` (+ 10 TB Block Vol) or `BM.DenseIO.E5.128` (US):**
   - `BM.Standard.E6.256`: **256c/512t AMD EPYC Turin, 3,072 GiB (3 TB) RAM**, **`$13.824/hr`** compute + `~$0.60/hr` 10 TB block = **`$14.42/hr`** (*published/derived*; **`$27.00` / 1,000 thread-hr** compute, *derived*), **`10 TB/mo` free egress**, **`$0` setup**. A 24-hour full-release clustering run costs **`~$346` total** (*derived*).
   - `BM.DenseIO.E5.128`: **128c/256t AMD EPYC Genoa, 1,536 GiB (1.5 TB) RAM, 81.6 TB local NVMe (`12 × 6.8 TB`)**, **`$15.60/hr`** (*published*), **`10 TB/mo` free egress**, **`$0` setup**. Zero scratch-disk bottleneck; 24-hour run = **`$374` total** (*derived*).
2. **Cherry Servers — AMD EPYC 9654 (US Chicago / EU):** **96c/192t, 768 GB DDR5**, `4+ TB` NVMe, **`100 TB` free egress**, **`$0` setup**, **`~$1.80–$2.00/hr`** or **`$1,000–$1,100+/mo`** (*published*; **`$7.85` / 1,000 thread-hr** monthly, *derived*). Best sustained US box if running weekly releases and 768 GB RAM suffices (custom Turin 9575F scales up to 3 TB RAM at `~$1,300+/mo`).
3. **OVHcloud US — `SCALE-a7` (Vint Hill, VA / Hillsboro, OR):** Dual EPYC 9654 (**192c/384t**), configurable to **1–1.5 TB DDR5 + 8–15 TB NVMe**, **unmetered 25 Gbps egress (`$0/TB`)**, base `~$1,402/mo` (`~$2,200–$3,200/mo` configured, *published*). Sensible *only* for `>= 3`-month commitments or zero-setup promos (month-to-month normally charges `1×` monthly rent setup).

> **Free UT Austin TACC Parallel Track:** Because the PI is UT Austin faculty, submit a free **TxRAS Lonestar6** (`128c / 256 GB` EPYC nodes + `5 TB` Corral) and **NSF ACCESS Stampede3** (`nvdimm` queue: **4 TB RAM / 80c** single-node jobs) Startup Allocation in parallel — native `module load tacc-apptainer` runs our pinned x86_64 `.sif` unchanged at `$0` compute and `$0` UT campus egress (§5).

---

## 1. Profile A Comparison — Stage-Build Workhorse
*Target: 64–192 threads (AMD EPYC Genoa/Turin preferred), 256–512 GB RAM, 4–8 TB local NVMe, cheap/unmetered egress (`~20–50 GB/stage`).*

| Provider & Location | Exact SKU / CPU / RAM / NVMe | Hourly $ (*published*) | Monthly $ (*published*) | Setup Fee (*published*) | Egress $/TB (*published*) | Root / Apptainer | Effective $/1,000 Thread-Hr (*derived*) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Hetzner Server Auction** (DE/FI — EU) | Refurb EPYC 7xxx/9xxx, **64–128t**, 256–512 GB RAM, 2–4× NVMe | Monthly only | **~$165–$495** (€150–€450) | **$0** (€0) | **$0/TB** (unmetered 1 Gbps) | Yes (full root) | **$4.71–$7.06** (`$330 / (730×96)×1000`) |
| **Hetzner Dedicated** (DE/FI — EU) | `AX162-R` / `AX162-1` (EPYC 9454P, **48c/96t**, 256 GB DDR5, 2× 1.92 TB + 1× 3.84 TB NVMe add-on) | Monthly only | **~$354** (legacy) to **~$785** (2026 `AX162-1`+NVMe) | **~$87–$97** (€79–€88) | **$0/TB** (unmetered 1 Gbps) | Yes (full root) | **$5.05–$11.20** (excl. setup) |
| **Cherry Servers** (US Chicago / EU) | **EPYC 9554P** (**64c/128t**, 384 GB DDR5, 2× 1.92 TB NVMe; upgradeable to 4+ TB) | **~$1.34/hr** | **$630–$783** | **$0** | **$0/TB** (100 TB/mo free) | Yes (full root) | **$8.03** monthly (`$750/(730×128)×1000`) / **$10.47** hourly |
| **OVHcloud US** (Vint Hill, VA / Hillsboro, OR) | `SCALE-a5` (EPYC 9554/9555, **64c/128t**, base 128 GB → 256–512 GB, 4 TB NVMe) | Monthly only | **~$840–$1,043** base (~$1,150–$1,400 cfg'd) | **1× mo rent** (m-t-m; $0 on 12-mo/promo) | **$0/TB** (unmetered 25 Gbps) | Yes (full root) | **$12.31** excl. setup (`$1,150/(730×128)×1000`); ~2× for 1-mo rental |
| **OVHcloud US** (Vint Hill / Hillsboro) | `SCALE-a6` (EPYC 9654/9655, **96c/192t**, base 128 GB → 384–512 GB, 4+ TB NVMe) | Monthly only | **~$1,033–$1,249** base (~$1,500 cfg'd) | **1× mo rent** (m-t-m; $0 on 12-mo/promo) | **$0/TB** (unmetered 25 Gbps) | Yes (full root) | **$10.70** excl. setup (`$1,500/(730×192)×1000`) |
| **Latitude.sh** (US Dallas / Ashburn / Chicago) | `m4.metal.large` (EPYC 9254, **24c/48t**, 384 GB DDR5, 2× 480 GB + 2× 3.8 TB NVMe = **7.6 TB**) | **~$0.95–$1.20/hr** | **$700–$850** | **$0** | **$0/TB** (20 TB/mo free) | Yes (full root) | **$22.50** (`$1.08/48×1000`) |
| **Oracle Cloud (OCI)** (US Dallas / Phoenix / Ashburn) | `VM.Standard.E5.Flex` (**64 OCPU / 128 vCPU** Genoa, 512 GB RAM + 4 TB Balanced Block Vol) | **$2.944/hr** compute (+ $0.24/hr block = **$3.18/hr**) | ~$2,323 (if 24/7) | **$0** | **$0/TB** (10 TB/mo free; $8.70/TB after) | Yes (root VM) | **$23.00** compute / **$24.84** w/ 4 TB block |
| **Vultr Bare Metal** (US **Dallas** — <200 mi from Austin) | `vbm-32c-755gb-amd` (EPYC 9354P, **32c/64t**, 768 GB RAM, 2× 480 GB NVMe — *needs Block Storage*) | **~$1.99/hr** | **$1,450** | **$0** | **$0/TB** (10 TB/mo free; $10/TB after) | Yes (full root) | **$31.09** (`$1.99/64×1000`) |
| **phoenixNAP BMC** (US Phoenix / Ashburn) | Dual-Xeon BMC (`d`-series; EPYC is CTO/dedicated only), 64–128t, 256–512 GB | **$0.80–$2.50/hr** | ~$584–$1,825 | **$0** | **$0/TB** (15 TB/mo free) | Yes (full root) | ~$12–$25 (Intel Xeon; EPYC custom quote) |
| **Equinix Metal** | *All SKUs* | — | — | — | — | — | **DISCONTINUED / EOL June 30, 2026** |

*Derived 100 Gbp Stage-Build Cost Estimate (assuming `~4,000–6,000` thread-hours for a full `annotated` stage incl. Pfam 38.2 + `50 GB` egress):*
- **Hetzner Auction (EU):** `~$20–$42` compute-share + `$0` egress (rented monthly at `$165–$495` flat).
- **Cherry Servers 9554P (US, hourly):** `~31–47` wall-clock hr × `$1.34/hr` = **`$42–$63`** + **`$0` egress**.
- **OCI `E5.Flex` 128t (US, hourly):** `~31–47` wall-clock hr × `$3.18/hr` = **`$99–$149`** + **`$0` egress**.

---

## 2. Profile B Comparison — Big-Release Clustering Monster
*Target: 128–256+ threads, 768 GB–3 TB RAM, 4–15+ TB fast local NVMe, cheap/unmetered egress (`~100–300 GB` release tar back to UT Austin).*

| Provider & Location | Exact SKU / CPU / RAM / NVMe | Hourly $ (*published*) | Monthly $ (*published*) | Setup Fee (*published*) | Egress $/TB (500 GB Cost) (*published*) | Root / Apptainer | Effective $/1,000 Thread-Hr (*derived*) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Oracle Cloud (OCI)** (US Ashburn / Phoenix / Dallas) | **`BM.Standard.E6.256`** (Dual EPYC 9J45 **Turin**, **256c/512t**, **3,072 GiB [3 TB] RAM**, 1.92 TB local + 10 TB Block Vol ~$0.60/hr) | **$13.824/hr** (+ $0.60/hr block = **$14.42/hr**) | ~$10,092 (if 24/7) | **$0** | **$0/TB** up to 10 TB/mo; $8.70/TB after (**$0** for 500 GB) | Yes (full root) | **$27.00** compute (`$13.824/512×1000`); **~$346 / 24-hr burst** |
| **Oracle Cloud (OCI)** (US Ashburn / Phoenix / Dallas) | **`BM.Standard.E5.192`** (Dual EPYC 9J14 Genoa, **192c/384t**, **2,304 GiB [2.3 TB] RAM**, + 10 TB Block Vol ~$0.60/hr) | **$10.368/hr** (+ $0.60/hr block = **$10.97/hr**) | ~$7,569 (if 24/7) | **$0** | **$0/TB** up to 10 TB/mo (**$0** for 500 GB) | Yes (full root) | **$27.00** compute (`$10.368/384×1000`); **~$263 / 24-hr burst** |
| **Oracle Cloud (OCI)** (US Ashburn / Phoenix / Dallas) | **`BM.DenseIO.E5.128`** (Dual EPYC 9J14 Genoa, **128c/256t**, **1,536 GiB [1.5 TB] RAM**, **81.6 TB local NVMe [12× 6.8 TB]**) | **$15.60/hr** | ~$11,388 (if 24/7) | **$0** | **$0/TB** up to 10 TB/mo (**$0** for 500 GB) | Yes (full root) | **$60.94** (`$15.60/256×1000`); **~$374 / 24-hr burst** |
| **Cherry Servers** (US Chicago / EU) | **EPYC 9654** (**96c/192t**, **768 GB DDR5**, 4+ TB NVMe; Turin 9575F up to 3 TB RAM custom) | **~$1.80–$2.00+/hr** | **$1,000–$1,100+** (Turin high-RAM: $1,300+) | **$0** | **$0/TB** (100 TB/mo free; **$0** for 500 GB) | Yes (full root) | **$7.85** monthly (`$1,100/(730×192)×1000`) / **$9.38** hourly |
| **OVHcloud US** (Vint Hill, VA / Hillsboro, OR) | **`SCALE-a7`** (Dual EPYC 9654, **192c/384t**, cfg'd to **1–1.5 TB DDR5 + 8–15 TB NVMe**) | Monthly only | **~$1,402** base (~$2,200–$3,200 cfg'd) | **1× mo rent** (m-t-m; $0 on 12-mo/promo) | **$0/TB** (unmetered 25 Gbps; **$0** for 500 GB) | Yes (full root) | **$9.28** excl. setup (`$2,600/(730×384)×1000`) |
| **Hetzner Dedicated** (DE/FI — EU) | **`AX162-3`** (EPYC 9454P, **48c/96t**, 512 GB DDR5 [max 1,152 GB], **15.36 TB NVMe [4× 3.84 TB]**) | Monthly only | **~$1,757** (€1,597.30; >€2,500 at 1,152 GB) | **~$878** (€799) | **$0/TB** (unmetered 1 Gbps; **$0** for 500 GB) | Yes (full root) | **$25.07** excl. setup (`$1,757/(730×96)×1000`) — *under-cored after 2026 DDR5 hike* |
| **Vultr Bare Metal** (US **Dallas**) | `vbm-64c-1536gb-amd` (**64c/128t**, **1.5 TB RAM**, 480 GB NVMe) / `vbm-128c-2048gb-amd` (**128c/256t** Milan, **2 TB RAM**, 960 GB NVMe) | **$4.01/hr** (64c) / **$7.53/hr** (128c) | **$2,925** (64c) / **$5,500** (128c) | **$0** | **$0/TB** up to 10 TB/mo; $10/TB after (**$0** for 500 GB) | Yes (full root) | **$31.30** (64c) / **$29.43** (128c) — *small local NVMe needs Block Storage* |
| **Latitude.sh** (US Dallas / Ashburn) | High-RAM EPYC 9554P (**64c/128t**, **1.5 TB DDR5**, local NVMe) | **~$3.79/hr** (*derived*) | **$2,770+** | **$0** | **$0/TB** (20 TB/mo free; **$0** for 500 GB) | Yes (full root) | **$29.65** (`$2,770/(730×128)×1000`) |
| **AWS EC2** (US East `us-east-1`) | **`r7a.metal-48xl`** (**192 vCPU** Genoa, **1.5 TB RAM**, EBS-only) | **$14.61/hr** OD / **$3.16–$4.50/hr** Spot | $10,662 OD | **$0** | **$90/TB** ($0.09/GB; **$45** for 500 GB) | Yes (full root) | **$76.07** OD / **$19.95** Spot — *Spot reclaim risk + egress trap* |
| **AWS EC2** (US East `us-east-1`) | **`i4i.metal`** (**128 vCPU** Xeon, **1 TB RAM**, **30 TB local NVMe**) | **$10.98/hr** OD / **$1.15–$3.87/hr** Spot | $8,016 OD | **$0** | **$90/TB** (**$45** for 500 GB) | Yes (full root) | **$85.80** OD / **$19.53** Spot — *Spot reclaim risk + egress trap* |
| **AWS EC2** (US East `us-east-1`) | **`i7ie.metal-48xl`** (**192 vCPU** Xeon, **1.5 TB RAM**, up to **120 TB Gen5 NVMe**) | **$24.95/hr** OD / **$5.67–$7.33/hr** Spot | $18,214 OD | **$0** | **$90/TB** (**$45** for 500 GB) | Yes (full root) | **$129.95** OD / **$33.85** Spot |

---

## 3. Why AWS Egress & Spot Are Traps vs. OCI / Cherry / OVH / Hetzner

1. **AWS Egress Trap:** AWS charges **`$0.09/GB` (`~$90/TB`)** outbound (`10 GB/mo` free). Syncing a `300 GB` release tar back to UT Austin costs `$27/release`; `500 GB` costs `$45`; a `2 TB` mirror costs `$180`. By contrast, **OCI** gives **`10 TB/mo` free** (`$8.70/TB` after), **Cherry Servers** gives **`100 TB/mo` free**, **Latitude.sh** gives **`20 TB/mo` free**, **Vultr** gives **`10 TB/mo` free**, and **OVHcloud / Hetzner** are unmetered (`$0/TB`).
2. **AWS `.metal` Spot Interruption Trap:** `.metal` Spot capacity is scarce and subject to 2-minute reclamation notices. A monolithic MMseqs2 `hybrid-b1` clustering step on `250M–500M` proteins cannot checkpoint mid-prefilter; losing a Spot node at hour 14 wastes the run. **OCI On-Demand Bare Metal** carries zero preemption risk at **`$10.37–$15.60/hr`** with 2–4× more threads and 1.5–3× more RAM than `i4i.metal`.

---

## 4. Niche / GPU-Neocloud Check

Evaluated **Rescale, Cirrascale, Massed Compute, Vast.ai, RunPod, Prime Intellect, and TensorDock**:
- **RunPod & Prime Intellect:** GPU-centric pods (`H100/H200/B200`); granular CPU-only sandboxes run `~$0.02/vCPU-hr + $0.0125/GiB-hr` (`~$15.36/hr` for `128 vCPU + 1 TB RAM` *derived*, pricier than OCI bare metal).
- **Vast.ai / TensorDock:** Decentralized consumer/prosumer marketplaces lacking reliable `1 TB+` RAM + multi-TB NVMe bare-metal SLAs.
- **Massed Compute / Cirrascale / Rescale:** Enterprise custom quotes or managed HPC SaaS markup over AWS/Azure/OCI. None beat OCI or Cherry Servers for self-managed `bio-box` runs.

---

## 5. Free Academic Parallel Track: UT Austin TACC

Because the PI is UT Austin faculty (`wilkcomp02.ccbb.utexas.edu`), submit a Startup Allocation in parallel with any commercial rental:
- **Lonestar6 (UT System exclusive, via TxRAS):** Dual-socket AMD EPYC (`128` cores / `256 GB` RAM per node) — ideal for parallel **Profile A** stage builds (`--depth called` / `annotated` across multiple nodes), plus up to `5 TB` free Corral storage.
- **Stampede3 (NSF ACCESS) & Frontera:** Stampede3 `nvdimm` queue provides **`4 TB` RAM / `80`-core Ice Lake nodes** (single-node jobs) and Frontera `nvdimm` provides **`2.1 TB` RAM / `112`-core nodes** — directly covering **Profile B** big-release clustering. *(Skip **Vista**: ARM Grace/Grace-Hopper nodes cannot run our x86_64 `.sif` without an ARM container rebuild.)*
- **Zero cost, zero egress, native Apptainer:** `module load tacc-apptainer` runs `bio-box_faec1f4.sif` rootless and unchanged, over UT campus / Internet2 peering (`$0` egress).
- **Trade-offs vs. commercial metal:** Slurm queue wait times, `24–48 h` job wall-clock caps, and shared Lustre `$SCRATCH` rather than local 80 TB NVMe RAID.

---

## 6. Zero-Corruption Architecture: Splitting `bio-box` Work Across an External Host

Because `bio-box` enforces **stage-partitioned immutability** (`AGENTS.md` §3.10 — a stage never rewrites an earlier stage's file), **manifest-driven reads** (§3.9), and **cross-host byte-identity** (§4.3), external parallelization requires **zero core-code changes** and **zero risk to `/stor` or `wilkcomp02`**.

### One-Time Host Bootstrap (~5 min on Ubuntu/Rocky bare metal)
1. Install `apptainer` (or `singularity-ce`), `rsync`, `zstd`, `pigz`.
2. Mount local NVMe drives as a single RAID-0 scratch volume at `/nvme` (`mdadm --create` + `mkfs.xfs`).
3. Sync `bio-box_faec1f4.sif` and verify `sha256sum == d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`; clone the pinned `bio-box` commit; copy `Pfam-A.hmm` (38.2).
4. Write a dedicated host profile `profiles/<host>.yaml` (mirroring `profiles/wilkcomp02.yaml`, with full single-tenant resources):
   ```yaml
   max_threads_day: 128          # or 256 / 512 for Profile B
   max_threads_night: 128
   mem_scope: 350G               # leave ~10% headroom for OS/page cache
   tmp: /nvme/bio-box
   durable: /nvme/durable
   bulk: /nvme/bulk
   nice_monolithic: 0
   nice_fanout: 0
   day_start_hour: 8
   day_end_hour: 20
   cpu_min_per_process: 999999   # single-tenant host: no shared-host watchdog
   cluster_cpu_min_max: 999999
   small_file_floor_bytes: 67108864
   require_mem_cap: false        # single-tenant dedicated box
   ```
5. **First-use byte-identity gate (`AGENTS.md` §4.3):** Build one small fixture stage (`micro-S0` or `S0`) on both `wilkcomp02` and the external host using the identical `.sif` and commit, and verify `sha256sum` matches across every Parquet output file before building production stages.

### Workflow 1 — Parallel Stage Builds (Profile A: Cherry Servers / Hetzner Auction / TACC Lonestar6)
1. **Lightweight read-only `DEC-26` root mirror (optional for overlap-free stages, required if contigs may overlap earlier stages):** Sync only the read-only `l0_seq/contig_index`, `l0_seq/seq_locator`, `manifests/`, `runs/` (and for Pfam domain copying `C-38h`, `l1_domains/` + `l1_proteins/`) from the canonical root — **never copying raw FASTAs or `seq_shards`** — and `chmod -R a-w` those earlier stage partitions so the new stage's `DEC-26` ownership matches the canonical root identically. (For stages with verified zero overlap against earlier stages, an isolated empty root produces identical stage partitions with `0` held contigs; verify zero `contig_uid` collision before landing into canonical `$ROOT`!)
2. **Stage inputs directly onto `/nvme/inputs`:** Inbound ingress is `$0` everywhere and saves UT campus bandwidth.
3. **Build stage-only (omit `--release`)** inside `bio-box_faec1f4.sif` using `profiles/<host>.yaml`.
4. **Package & checksum:** Tar only `stage=<STAGE>/` partitions + `runs/stage=<STAGE>/` + ledger into `added_<STAGE>.tar` + `added_<STAGE>.tar.sha256` (`~20–50 GB`).
5. **Sync & atomic land:** `rsync --partial --inplace` to `wilkcomp02:/tmp/if2453/incoming/` (plus durable mirror to `tower` and `/stor/work`). Verify `sha256sum -c`, verify zero path collisions and zero `contig_uid` ownership collisions against canonical `$ROOT`, then move the new `stage=<STAGE>/` directories into `$ROOT` atomically.

### Workflow 2 — Burst Big-Release Clustering (Profile B: OCI `BM.Standard.E6.256` / `BM.DenseIO.E5.128` / TACC Stampede3 `nvdimm`)
1. **Spin up OCI bare metal** (`BM.Standard.E6.256` + 10 TB Block Vol, or `BM.DenseIO.E5.128` with 81.6 TB local NVMe).
2. **Sync read-only stage inputs:** `rsync` only the stage partitions needed for release clustering (`l1_proteins`, `l1_tapes`, `l1_domains`, `l0_seq/{contig_index,seq_locator,genome_input}`, `contig_label`, `runs/`, `manifests/` — **no `seq_shards` needed unless running sequence-oracle exports on-box**; inbound ingress is `$0`). Lock stage partitions `chmod -R a-w`.
3. **Run release build + search pack** (`CLUSTER_METHOD=hybrid-b1`, `DEC-77`) at full `256–512` threads and `1.5–3 TB` RAM.
4. **Verify worked examples on-box** (when `seq_shards` are present, or run `release_examples.sh` immediately after landing on `wilkcomp02`).
5. **Tar & sync back release-only outputs:** Tar `release=<REL>/` + `manifests/release=<REL>.json` + `runs/stage=<REL>/` + search pack (`~100–300 GB` — **`$0` egress** under OCI's `10 TB/mo` free tier) + `.sha256`, verify `sha256sum -c`, land atomically on `wilkcomp02` + `tower` + `/stor`.
6. **Terminate OCI instance immediately:** Total burst cost `~$346–$560` (*derived*); ongoing idle burn `$0`.

---

## 7. PI Decisions Needed (`A-28`, Tier H)

1. **Q1 (Profile B Burst Clustering):** Approve opening a pay-as-you-go **OCI** account for `24–48 h` big-release bursts (`BM.Standard.E6.256` @ `$14.42/hr` or `BM.DenseIO.E5.128` @ `$15.60/hr`, `$0` setup, `10 TB/mo` free egress, `~$350–$560/release`) once releases outgrow `wilkcomp02`'s `1.1 TiB` `/tmp` (expected around `r26.11.5`–`r26.11.6`, `> 350–500 Gbp`)?
2. **Q2 (Profile A Parallel Stage Builds):** If parallel stage builds beyond `wilkcomp02` are wanted, choose **(a) Cherry Servers EPYC 9554P** (US Chicago, `$1.34/hr` or `$630–$783/mo`, `$0` setup, `100 TB` free egress), **(b) Hetzner Server Auction** (EU, `~$165–$495/mo`, `€0` setup, unmetered egress), **(c) OCI `VM.Standard.E5.Flex`** on-demand (`$3.18/hr`, `$0` idle), or **(d) defer Profile A until TACC Lonestar6 responds**?
3. **Q3 (Free UT Austin TACC Track):** Submit a free **TxRAS Lonestar6** + **NSF ACCESS Stampede3 (`nvdimm` 4 TB)** Startup Allocation now?
4. **Q4 (Repo Profile & Sync Helper):** Once Q1–Q3 are decided, authorize a follow-on **L** card (`A-28a`) to commit `profiles/<host>.yaml`, a `scripts/ops/stage_sync_export.sh` tar+checksum+collision-check helper, and the `S0` cross-host byte-identity gate.
