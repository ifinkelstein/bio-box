# Storage: a decision set for the PI (draft)

Drafted 2026-11-20 by the orchestrator at the PI's word: "draft a decision set for me to re storage. ideally, tower
holds a copy and removes pressure from pod. also can backup on tacc/corral". **Nothing here is decided.** Each item
has a recommendation; the PI answers by number, and each answer becomes a row of `planning/DECISIONS.md`. Until then
STATUS open item 19 stands. Numbers carry the labels of AGENTS §6.

## What we have

| Place | Size, free | What it is | What bio-box holds there |
|---|---|---|---|
| Pod `/tmp` | 1.1 TB, 556 GB free (*measured*, `df`, wilkcomp02, 2026-11-20 15:17) | local disk, shared, cleared at a reboot (2026-11-17) | the live data root (about 110 GB), search packs (59 GB), staged stage inputs (about 133 GB: `SPX1` 71, `PL2` 33, `V2` 9.2, `PL1` 7.0, `V1` 5.0, `GV1` 4.4), other inputs and raw downloads (`s2/phase2` 16, `s2/dl` 7.0, `s3` 57) (*measured*, `du`, 2026-11-20) |
| `/stor/work` | 68 TB, 19 TB free | the group's, backed up | 1.2 TB: the JGI archives (1.14 TB) and the root's seven tars (107 GB). The footprint draft says "< 50 GB" |
| `/stor/scratch` | 26 TB, 19 TB free | the group's, not backed up | 164 GB: earlier stage inputs, reference libraries. The footprint draft says 1.5 to 2 TB |
| `tower` | 32 TB array, 18 TB free (*measured*, 2026-11-20) | the PI's NAS at home; 8 cores, 31 GB; the pod cannot reach it, the laptop relays at 9.70 MB/s (*measured*, one stream) | the second copy of the root's tars (being mirrored today) |
| TACC Corral | first 5 TB free per PI on one project for UT System researchers; then $60 per TB per year (*published*, TACC's Corral guide, read 2026-11-20) | TACC's project storage: not purged, a tape replica in a second datacenter, yearly renewal, 1 million files per 5 TB, reached by `scp` to `data.tacc.utexas.edu`, **mounted on Lonestar6 login and compute nodes** | nothing yet; no allocation requested |

Two facts fix what is possible (card C-49, `planning/design/C-49_tiers.md`): (1) **anything a release uses must be on
the machine that computes**: the builders find earlier stages by what is on the pod's disk, and a release must verify
every file it lists; so `tower` and Corral can hold copies, inputs and archives, never the only copy of a built
stage; (2) the relay to `tower` needs the laptop awake at home, and a file interrupted mid-stream starts again (the
largest JGI archive is 525 GB, about 15 hours).

## The decisions

**S-1. What each place is for.** *Recommended*: pod `/tmp` = the live root, the packs of the newest two releases, and
the inputs of stages not yet built; `/stor/work` = one durable copy of the root's tars and nothing else;
`/stor/scratch` = reference libraries and re-fetchable inputs; `tower` = a second copy of the root's tars and a copy
of everything that otherwise exists once (source archives, staged inputs, the deferred-data lists); Corral = a third
copy off campus hardware, and later the staging area for TACC compute.

**S-2. `tower` gets a copy of everything that exists once.** In this order: (a) the root's tars (running, 3 of 7
verified at 15:15); (b) the deferred-data lists and skipped-record tables that are on `/tmp` only (`V2`'s two lists,
165 MB; `PL2`'s skipped records); (c) each stage's staged inputs before it is built; (d) the source archives.
Every file is hashed on arrival against its recorded sha256 and listed in a `SHA256SUMS` there; a re-read of
everything on `tower` twice a year. *Recommended*: yes to (a), (b), (c) now; (d) under S-7.

**S-3. Use the `tower` copy to take pressure off the pod's `/tmp`.** Once a stage is built, in a verified release,
and its inputs are verified on `tower`, its staged inputs and raw downloads are removed from `/tmp`. Over the stages
now queued that frees about 133 GB of inputs, the 51 GB of `PL2`'s raw downloads under `s3`, and `s2/phase2`'s 16 GB
after `C1` to `C6` (*derived*: about 200 GB in all). And: keep search packs for the newest two releases only (a pack
is made again from its release in hours); that frees 33.5 GB today (*derived*, card C-49) and one pack per release
after. *Recommended*: yes to both. This is the most relief the pod can get without code changes.

**S-4. Old releases' large files** (DEC-52 left this open). Keep every release's manifest and small records always;
keep the large `release=` files (clusters, tokens, postings) on the pod for the newest two releases and any signed
`rS` release; an older release's large files go to `tower` and are then removed from the pod, and the release is
marked retired (still listed, restorable, not searchable). This needs code (nine functions named in the C-49 note)
and a change to the rule that the root is append-only. Frees 61 to 76 GB now (*derived*), about 25 to 60 GB per
release after. *Recommended*: decide the number (two) now, authorise the card, and remove nothing until it is merged.

**S-5. The 1.14 TB of JGI archives on `/stor/work`.** They put us 24 times over the draft footprint there. Options:
(a) leave them and tell the pod owner; (b) the PI moves them to `/stor/scratch` (inside the 1.5 to 2 TB asked for);
(c) once two other verified copies exist (S-6, S-7), the PI removes them from `/stor`. *Recommended*: (b) now, since
it needs no transfer and they are re-fetchable with the PI's JGI login; (c) later if `/stor` space is wanted.
Agents do not move or delete under `/stor`; this is the PI's hand.

**S-6. TACC Corral.** Request the free 5 TB allocation now, in the TACC User Portal, alongside the two Lonestar6
requests already drafted (card C-37a). Put there: the root's tars (107 GB now, 25 to 60 GB more per release) and the
source archives (1.14 TB): about 1.3 TB, well inside 5 TB and its file limit (the tars and archives are tens of
files). Because Corral is mounted on Lonestar6's compute nodes, the same copy is what a TACC build would start from.
Beyond 5 TB it is $60 per TB per year: holding the built SPIRE corpus too (about 3 TB at depth `called`, *derived*)
would be about 10 TB, $300 a year. *Recommended*: yes, 5 TB now; decide on more after the TACC trial (card C-37).

**S-7. How the big files travel, and with what credentials.** For `tower`: the laptop relay, no key on `tower`
(the PI's choice of 2026-11-20), for everything up to about 100 GB per sitting. For the pod to Corral: `scp` from the
pod to `data.tacc.utexas.edu` inside a session the PI opens with the TACC token (as card C-37 plans for logins); no
TACC credential is stored, and nothing is installed on the pod. The rate is not measured; both ends are on the
university's network. For the two JGI archives too large for one laptop sitting (409 GB and 525 GB): either (i) a
restricted read-only pull key on `tower` for this one job (the exact lines are in the C-49 note), or (ii) they do not
go to `tower`, and Corral plus `/stor` are their two copies. *Recommended*: (ii).

**S-8. What this does not solve.** The live root still has to fit where the compute is. On the pod that is about
1 TB shared; the stages queued today fit with S-3 and S-4; SPIRE beyond two or three tranches does not. The answer
for the body of SPIRE is TACC (build on Lonestar6 with its scratch, keep on Corral), and it waits on the two
allocation requests and the trial of card C-37. *Recommended*: no decision now; submit the requests.

## If the PI says yes to all as recommended

1. Today: the mirror finishes; the lists of S-2 (b) go to `tower` with the new script.
2. This week: inputs of `V2`, `GV1`, `PL2`, `SPX1` to `tower` before their builds (about 118 GB, about 3.4 hours of
   relay, *derived* at 9.70 MB/s); old packs dropped (33.5 GB); the PI moves the JGI archives to `/stor/scratch` and
   sends the pod owner a revised footprint (`/stor/work` under 0.5 TB, `/stor/scratch` under 2 TB, `/tmp` up to
   600 GB with 200 GB always left free); the PI requests Corral with the Lonestar6 requests.
3. When Corral is granted: tars and archives there, in a session the PI opens.
4. A card for S-4's retired releases; card C-37's trial; then the SPIRE decision of S-8.
