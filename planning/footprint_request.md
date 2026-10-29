# Footprint request to the pod owner (draft for the PI to send)

Status: draft, approved in principle (DEC-10, 2026-10-29). Not sent. File the reply here when it arrives.

**To:** Claus Wilke
**Cc:** rctf-support (only for the `/tmp` question at the end)
**Subject:** wilkcomp02 footprint request: bio-box genome corpus (staged, through ~mid-2027)

Hi Claus,

I'd like to run a staged build of a genome corpus (bio-box) on wilkcomp02 and want to agree the footprint with you before anything large starts. It grows in stages, each with a pause point, so the load is predictable and can stop any time.

**What it is.** Phage, plasmid and prokaryotic genomes (GTDB R226 and public collections), annotated per genome and indexed for gene-neighborhood and fusion searches. Early stages are small; the largest stage that runs fully here is all 143,614 GTDB species representatives. Anything bigger moves to a rented server.

**What I'm asking for** (estimates; each stage measures its real use and I'll report it):

| | Request |
|---|---|
| Threads | ≤40 physical-core threads by day, ≤56 at night and weekends, never more |
| Memory | Every job in its own `systemd-run` scope with `MemoryMax=200G`; one large job at a time by day. One exception to discuss later (not now): a clustering step at the largest stage may need ~160 GB for several hours |
| Priority | Long single-process tools start at `nice 19`; fan-out workers at `nice 15` and short enough (< 300 CPU-min each) not to trip the watchdog; all I/O at `ionice -c 2 -n 7` |
| `/stor/scratch/Wilke/if2453/all-the-data` | grows to ~0.1 TB over the first months, ~0.3 TB at the complete-genomes stage, ~1.5–2 TB at the species-representatives stage. Only large files (≥64 MB); no per-genome small files |
| `/stor/work/Wilke/if2453/all-the-data` | small and durable: manifests, checksums, release tables, < 50 GB |
| Local `/tmp` | working space per job; up to ~500–900 GB at the largest stage only, cleaned after each pass |
| Run lengths | days for early stages; the species-representatives stage is CPU-bound over several weeks, in ~100 resumable shards |
| Pausing | a STOP file pauses everything at the next shard boundary within minutes; you or anyone in the lab can create it. Each directory has a README saying what it is, who owns it and when it can be deleted |

**Two further asks:**
1. Occasional GPU time on wilkcomp03 later in the project (a protein-structure step, roughly 0.7–1.2k GPU-hours in total, *estimate*), on whatever terms suit you.
2. A second account for a postdoc, so the project doesn't depend on one login.

**For rctf-support:** is `/tmp` on wilkcomp02 cleaned on reboot or by a timer? Jobs checkpoint to NFS either way; I'd like to know the policy.

Happy to adjust any of this. Thanks,
Ilya
