# Cluster id churn between releases when each release clusters from scratch (wilkcomp02, 2026-11-10)

DEC-24 decided that each release clusters all its proteins from scratch and left open "how many
cluster ids change from one release to the next". This note measures that on the three provisional
releases of the import ladder (`rL.1` phages; `rL.2` plus 250 prokaryotes; `rL.3` plus 2,247 more;
see `2026-11-10_import_ladder/README.md`). It changes no code and recommends no design.

Every number below is *measured* on wilkcomp02 on 2026-11-10 by the script at the end, unless it is
marked *derived* (formula given). Fractions in the tables are *derived* (count divided by the
denominator named in the table title). Host: `wilkcomp02.ccbb.utexas.edu`. Data (read only):
`/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`, built by the pre-B-17 code at `f9af5b6`
inside `bio-box-test_74b3e05.sif`.

## What determines `cluster_uid` (quoted from `src/biobox/l2/cluster.py`)

```python
def cluster_uid(member_uids: list[bytes], phash: str) -> bytes:
    return blake3(b"".join(sorted(member_uids)) + phash.encode()).digest(length=16)
```

and in `_cluster_members` (members arrive sorted by (rep, member); "each cluster's sorted members are
hashed when its last one has passed"). The old code that built `rL.*` has the identical function
(`~/bio-box-next/src/biobox/l2/cluster.py` line 50). `phash` is `params_hash(level)`: a hash of the
level, identity and coverage thresholds and the MMseqs2 flags, so it is constant across these releases.

So the id is a hash of the cluster's complete member set, not of its representative. Consequences for
reading the numbers:

- A cluster keeps its id if and only if its member set is exactly unchanged. Gaining one new member
  (the normal effect of adding genomes) changes the id, even when nothing else moved. The first rows
  of every table ("cluster_uid still present") therefore are "cluster did not gain, lose or swap any
  member".
- The representative is not an input of the id. A representative can change with the id unchanged
  (the "uid same, rep changed" rows) and the id can change with the representative unchanged
  (the "uid changed, rep same" rows).
- `rep_uid` is whatever MMseqs2 chose (`--cluster-mode 2`, greedy by length, `--cov-mode 1`).

## Method

Inputs are the `cluster_member` files each manifest lists for its own release (partition
`release == R`; one file per release). Note that `manifests/release=rL.3.json` also lists a file
`release=rL.2p` (the contaminated build of the ladder README, finding 1); it is not used. Schema of
`cluster_member`: `member_uid BLOB, level VARCHAR ('nr90'|'nr50'), cluster_uid BLOB, rep_uid BLOB,
cluster_params_hash VARCHAR, release VARCHAR`. Row counts per level equal protein counts (no
duplicate `(level, member_uid)` rows in any file, checked below).

For an earlier release A and a later release B, at one level: `j` = A joined to B on `member_uid`
(the carried proteins; in all four pairs it is every protein of A). `p(X, Y, k)` = number of carried
proteins in earlier cluster X and later cluster Y. For earlier cluster X, `nc(X)` = carried members
(= its size, since every protein is carried). For later cluster Y, `mb(Y)` = carried members.

Classes of an earlier cluster, exclusive and exhaustive, by where its carried members went:

1. **same members**: one later cluster Y with `mb(Y) = nc(X)` and `size_B(Y) = nc(X)` (id kept).
2. **same carried, gained**: one Y with `mb(Y) = nc(X)` and `size_B(Y) > nc(X)` (only gained members).
3. **split**: carried members in two or more later clusters.
4. **merged**: one Y with `mb(Y) > nc(X)`, i.e. Y also holds carried members of another earlier
   cluster (not split).

The memo rows count clusters that touch members of another earlier cluster in any destination
(including after a split). Jaccard of an (X, Y) pair is `J = k / (nc(X) + mb(Y) - k)`, over carried
proteins only: new proteins in Y are ignored, so a cluster that only gained members has J = 1.0.
"Best" is the maximum J over the later clusters that hold any of X's members.

Size filters (`>= 2`, `>= 10`) use the earlier-release cluster size, for clusters and for the
proteins in them.

## Headline (derived from the tables below: one minus "cluster_uid still present", and the protein rows)

| Pair | Level | Earlier clusters losing their id | Carried proteins whose cluster_uid changed | Carried proteins whose rep changed | Proteins in their old cluster's largest fragment |
|---|---|---|---|---|---|
| rL.1 to rL.2 | NR90 | 0.21% (1,138 of 534,101) | 0.45% | 0.12% | 99.97% |
| rL.1 to rL.2 | NR50 | 2.39% (7,093 of 297,082) | 5.09% | 1.14% | 99.73% |
| rL.2 to rL.3 | NR90 | 12.44% (173,683 of 1,395,712) | 11.40% | 7.06% | 99.63% |
| rL.2 to rL.3 | NR50 | 39.94% (304,276 of 761,833) | 43.98% | 30.89% | 94.39% |

Of the earlier clusters that lost their id, the share that only gained members (class 2) is, for
rL.2 to rL.3: NR90 92.5% (160,588 of 173,683), NR50 74.8% (227,581 of 304,276); for rL.1 to rL.2:
NR90 74.5% (848 of 1,138), NR50 73.7% (5,225 of 7,093).

## Results

Everything from here to the script is the unedited output of the script (`result.md`), which
prints the same tables for every pair and level, for all clusters and for earlier size >= 2 and
>= 10. The step times printed (9 to 30 s per pair and level) are *measured* wall times.

## Schema of cluster_member (rL.1 file)

```
member_uid BLOB
level VARCHAR
cluster_uid BLOB
rep_uid BLOB
cluster_params_hash VARCHAR
release VARCHAR
```

Files per release (from the manifests, partition release == R):

- rL.1: l2/cluster_member/release=rL.1/part-0.parquet
- rL.2: l2/cluster_member/release=rL.2/part-0.parquet
- rL.3: l2/cluster_member/release=rL.3/part-0.parquet

# rL.1 -> rL.2, level nr90

Sizes (measured)

| quantity | value |
|---|---|
| proteins in rL.1 | 1,134,550 |
| proteins in rL.2 | 2,053,075 |
| carried (in both) | 1,134,550 |
| carried / earlier | 1.0000 |
| clusters in rL.1 | 534,101 |
| clusters in rL.2 | 1,395,712 |
| duplicate member rows (a, b) | 0, 0 |

## rL.1 -> rL.2 nr90: all clusters


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 534,101 | 1.0000 |
| cluster_uid still present in later release | 532,963 | 0.9979 |
| same carried members, none gained (class 1) | 532,963 | 0.9979 |
| same carried members, gained new members (class 2) | 848 | 0.0016 |
| split: carried members in >= 2 later clusters (class 3) | 124 | 0.0002 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 166 | 0.0003 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 248 | 0.0005 |
| memo: split and also merged | 82 | 0.0002 |
| memo: earlier-release proteins in these clusters | 1,134,550 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 1,134,230 | 0.9997 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 1,134,550 | 1.0000 |
| cluster_uid changed | 5,104 | 0.0045 |
| representative changed | 1,364 | 0.0012 |
| uid same, rep changed | 11 | 0.0000 |
| uid changed, rep same | 3,751 | 0.0033 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 1,138 | 1.0000 |
| best J = 1.0 | 848 | 0.7452 |
| best J >= 0.9 | 874 | 0.7680 |
| best J >= 0.5 | 1,000 | 0.8787 |
| best J < 0.5 | 138 | 0.1213 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.3333 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 1,000 | 0.8787 |

## rL.1 -> rL.2 nr90: earlier size >= 2


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 153,266 | 1.0000 |
| cluster_uid still present in later release | 152,780 | 0.9968 |
| same carried members, none gained (class 1) | 152,780 | 0.9968 |
| same carried members, gained new members (class 2) | 316 | 0.0021 |
| split: carried members in >= 2 later clusters (class 3) | 124 | 0.0008 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 46 | 0.0003 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 128 | 0.0008 |
| memo: split and also merged | 82 | 0.0005 |
| memo: earlier-release proteins in these clusters | 753,715 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 753,395 | 0.9996 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 753,715 | 1.0000 |
| cluster_uid changed | 4,452 | 0.0059 |
| representative changed | 1,000 | 0.0013 |
| uid same, rep changed | 11 | 0.0000 |
| uid changed, rep same | 3,463 | 0.0046 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 486 | 1.0000 |
| best J = 1.0 | 316 | 0.6502 |
| best J >= 0.9 | 342 | 0.7037 |
| best J >= 0.5 | 448 | 0.9218 |
| best J < 0.5 | 38 | 0.0782 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.5000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 448 | 0.9218 |

## rL.1 -> rL.2 nr90: earlier size >= 10


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 13,552 | 1.0000 |
| cluster_uid still present in later release | 13,463 | 0.9934 |
| same carried members, none gained (class 1) | 13,463 | 0.9934 |
| same carried members, gained new members (class 2) | 31 | 0.0023 |
| split: carried members in >= 2 later clusters (class 3) | 43 | 0.0032 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 15 | 0.0011 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 53 | 0.0039 |
| memo: split and also merged | 38 | 0.0028 |
| memo: earlier-release proteins in these clusters | 307,047 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 306,842 | 0.9993 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 307,047 | 1.0000 |
| cluster_uid changed | 3,061 | 0.0100 |
| representative changed | 561 | 0.0018 |
| uid same, rep changed | 0 | 0.0000 |
| uid changed, rep same | 2,500 | 0.0081 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 89 | 1.0000 |
| best J = 1.0 | 31 | 0.3483 |
| best J >= 0.9 | 57 | 0.6404 |
| best J >= 0.5 | 80 | 0.8989 |
| best J < 0.5 | 9 | 0.1011 |
| median best J | 0.9444 |  |
| 10th percentile best J | 0.4860 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 80 | 0.8989 |

Later-release clusters (rL.2) with an id not present in rL.1: best-overlap earlier cluster (Jaccard over carried proteins)

| quantity | count | fraction |
|---|---|---|
| clusters in rL.2 | 1,395,712 | 1.0000 |
| id not present in earlier release | 862,749 | 0.6181 |
| of them with no carried member (wholly new proteins) | 861,597 | 0.6173 |
| of them with carried members | 1,152 | 0.0008 |
|   best J = 1.0 | 848 | 0.7361 |
|   best J >= 0.9 | 874 | 0.7587 |
|   best J >= 0.5 | 1,013 | 0.8793 |
|   best J < 0.5 | 139 | 0.1207 |

(step time 12 s)

# rL.1 -> rL.2, level nr50

Sizes (measured)

| quantity | value |
|---|---|
| proteins in rL.1 | 1,134,550 |
| proteins in rL.2 | 2,053,075 |
| carried (in both) | 1,134,550 |
| carried / earlier | 1.0000 |
| clusters in rL.1 | 297,082 |
| clusters in rL.2 | 761,833 |
| duplicate member rows (a, b) | 0, 0 |

## rL.1 -> rL.2 nr50: all clusters


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 297,082 | 1.0000 |
| cluster_uid still present in later release | 289,989 | 0.9761 |
| same carried members, none gained (class 1) | 289,989 | 0.9761 |
| same carried members, gained new members (class 2) | 5,225 | 0.0176 |
| split: carried members in >= 2 later clusters (class 3) | 1,080 | 0.0036 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 788 | 0.0027 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 1,410 | 0.0047 |
| memo: split and also merged | 622 | 0.0021 |
| memo: earlier-release proteins in these clusters | 1,134,550 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 1,131,539 | 0.9973 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 1,134,550 | 1.0000 |
| cluster_uid changed | 57,710 | 0.0509 |
| representative changed | 12,889 | 0.0114 |
| uid same, rep changed | 209 | 0.0002 |
| uid changed, rep same | 45,030 | 0.0397 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 7,093 | 1.0000 |
| best J = 1.0 | 5,225 | 0.7366 |
| best J >= 0.9 | 5,558 | 0.7836 |
| best J >= 0.5 | 6,742 | 0.9505 |
| best J < 0.5 | 351 | 0.0495 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.5000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 6,742 | 0.9505 |

## rL.1 -> rL.2 nr50: earlier size >= 2


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 122,319 | 1.0000 |
| cluster_uid still present in later release | 118,395 | 0.9679 |
| same carried members, none gained (class 1) | 118,395 | 0.9679 |
| same carried members, gained new members (class 2) | 2,320 | 0.0190 |
| split: carried members in >= 2 later clusters (class 3) | 1,080 | 0.0088 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 524 | 0.0043 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 1,146 | 0.0094 |
| memo: split and also merged | 622 | 0.0051 |
| memo: earlier-release proteins in these clusters | 959,787 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 956,776 | 0.9969 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 959,787 | 1.0000 |
| cluster_uid changed | 54,541 | 0.0568 |
| representative changed | 10,793 | 0.0112 |
| uid same, rep changed | 209 | 0.0002 |
| uid changed, rep same | 43,957 | 0.0458 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 3,924 | 1.0000 |
| best J = 1.0 | 2,320 | 0.5912 |
| best J >= 0.9 | 2,653 | 0.6761 |
| best J >= 0.5 | 3,714 | 0.9465 |
| best J < 0.5 | 210 | 0.0535 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.5000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 3,714 | 0.9465 |

## rL.1 -> rL.2 nr50: earlier size >= 10


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 20,570 | 1.0000 |
| cluster_uid still present in later release | 19,473 | 0.9467 |
| same carried members, none gained (class 1) | 19,473 | 0.9467 |
| same carried members, gained new members (class 2) | 453 | 0.0220 |
| split: carried members in >= 2 later clusters (class 3) | 437 | 0.0212 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 207 | 0.0101 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 523 | 0.0254 |
| memo: split and also merged | 316 | 0.0154 |
| memo: earlier-release proteins in these clusters | 600,061 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 597,921 | 0.9964 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 600,061 | 1.0000 |
| cluster_uid changed | 43,916 | 0.0732 |
| representative changed | 6,652 | 0.0111 |
| uid same, rep changed | 182 | 0.0003 |
| uid changed, rep same | 37,446 | 0.0624 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 1,097 | 1.0000 |
| best J = 1.0 | 453 | 0.4129 |
| best J >= 0.9 | 782 | 0.7129 |
| best J >= 0.5 | 1,038 | 0.9462 |
| best J < 0.5 | 59 | 0.0538 |
| median best J | 0.9762 |  |
| 10th percentile best J | 0.6366 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 1,038 | 0.9462 |

Later-release clusters (rL.2) with an id not present in rL.1: best-overlap earlier cluster (Jaccard over carried proteins)

| quantity | count | fraction |
|---|---|---|
| clusters in rL.2 | 761,833 | 1.0000 |
| id not present in earlier release | 471,844 | 0.6194 |
| of them with no carried member (wholly new proteins) | 464,356 | 0.6095 |
| of them with carried members | 7,488 | 0.0098 |
|   best J = 1.0 | 5,225 | 0.6978 |
|   best J >= 0.9 | 5,558 | 0.7423 |
|   best J >= 0.5 | 6,839 | 0.9133 |
|   best J < 0.5 | 649 | 0.0867 |

(step time 9 s)

# rL.2 -> rL.3, level nr90

Sizes (measured)

| quantity | value |
|---|---|
| proteins in rL.2 | 2,053,075 |
| proteins in rL.3 | 10,620,684 |
| carried (in both) | 2,053,075 |
| carried / earlier | 1.0000 |
| clusters in rL.2 | 1,395,712 |
| clusters in rL.3 | 8,278,537 |
| duplicate member rows (a, b) | 0, 0 |

## rL.2 -> rL.3 nr90: all clusters


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 1,395,712 | 1.0000 |
| cluster_uid still present in later release | 1,222,029 | 0.8756 |
| same carried members, none gained (class 1) | 1,222,029 | 0.8756 |
| same carried members, gained new members (class 2) | 160,588 | 0.1151 |
| split: carried members in >= 2 later clusters (class 3) | 6,040 | 0.0043 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 7,055 | 0.0051 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 8,721 | 0.0062 |
| memo: split and also merged | 1,666 | 0.0012 |
| memo: earlier-release proteins in these clusters | 2,053,075 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 2,045,520 | 0.9963 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 2,053,075 | 1.0000 |
| cluster_uid changed | 233,994 | 0.1140 |
| representative changed | 144,926 | 0.0706 |
| uid same, rep changed | 61 | 0.0000 |
| uid changed, rep same | 89,129 | 0.0434 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 173,683 | 1.0000 |
| best J = 1.0 | 160,588 | 0.9246 |
| best J >= 0.9 | 160,672 | 0.9251 |
| best J >= 0.5 | 171,368 | 0.9867 |
| best J < 0.5 | 2,315 | 0.0133 |
| median best J | 1.0000 |  |
| 10th percentile best J | 1.0000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 171,368 | 0.9867 |

## rL.2 -> rL.3 nr90: earlier size >= 2


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 191,405 | 1.0000 |
| cluster_uid still present in later release | 158,237 | 0.8267 |
| same carried members, none gained (class 1) | 158,237 | 0.8267 |
| same carried members, gained new members (class 2) | 25,640 | 0.1340 |
| split: carried members in >= 2 later clusters (class 3) | 6,040 | 0.0316 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 1,488 | 0.0078 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 3,154 | 0.0165 |
| memo: split and also merged | 1,666 | 0.0087 |
| memo: earlier-release proteins in these clusters | 848,768 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 841,213 | 0.9911 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 848,768 | 1.0000 |
| cluster_uid changed | 93,479 | 0.1101 |
| representative changed | 57,795 | 0.0681 |
| uid same, rep changed | 61 | 0.0001 |
| uid changed, rep same | 35,745 | 0.0421 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 33,168 | 1.0000 |
| best J = 1.0 | 25,640 | 0.7730 |
| best J >= 0.9 | 25,724 | 0.7756 |
| best J >= 0.5 | 32,486 | 0.9794 |
| best J < 0.5 | 682 | 0.0206 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.5000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 32,486 | 0.9794 |

## rL.2 -> rL.3 nr90: earlier size >= 10


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 13,785 | 1.0000 |
| cluster_uid still present in later release | 13,260 | 0.9619 |
| same carried members, none gained (class 1) | 13,260 | 0.9619 |
| same carried members, gained new members (class 2) | 336 | 0.0244 |
| split: carried members in >= 2 later clusters (class 3) | 164 | 0.0119 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 25 | 0.0018 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 119 | 0.0086 |
| memo: split and also merged | 94 | 0.0068 |
| memo: earlier-release proteins in these clusters | 310,254 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 309,646 | 0.9980 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 310,254 | 1.0000 |
| cluster_uid changed | 11,582 | 0.0373 |
| representative changed | 2,851 | 0.0092 |
| uid same, rep changed | 41 | 0.0001 |
| uid changed, rep same | 8,772 | 0.0283 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 525 | 1.0000 |
| best J = 1.0 | 336 | 0.6400 |
| best J >= 0.9 | 408 | 0.7771 |
| best J >= 0.5 | 504 | 0.9600 |
| best J < 0.5 | 21 | 0.0400 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.6918 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 504 | 0.9600 |

Later-release clusters (rL.3) with an id not present in rL.2: best-overlap earlier cluster (Jaccard over carried proteins)

| quantity | count | fraction |
|---|---|---|
| clusters in rL.3 | 8,278,537 | 1.0000 |
| id not present in earlier release | 7,056,508 | 0.8524 |
| of them with no carried member (wholly new proteins) | 6,880,803 | 0.8312 |
| of them with carried members | 175,705 | 0.0212 |
|   best J = 1.0 | 160,588 | 0.9140 |
|   best J >= 0.9 | 160,672 | 0.9144 |
|   best J >= 0.5 | 172,248 | 0.9803 |
|   best J < 0.5 | 3,457 | 0.0197 |

(step time 30 s)

# rL.2 -> rL.3, level nr50

Sizes (measured)

| quantity | value |
|---|---|
| proteins in rL.2 | 2,053,075 |
| proteins in rL.3 | 10,620,684 |
| carried (in both) | 2,053,075 |
| carried / earlier | 1.0000 |
| clusters in rL.2 | 761,833 |
| clusters in rL.3 | 2,966,200 |
| duplicate member rows (a, b) | 0, 0 |

## rL.2 -> rL.3 nr50: all clusters


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 761,833 | 1.0000 |
| cluster_uid still present in later release | 457,557 | 0.6006 |
| same carried members, none gained (class 1) | 457,557 | 0.6006 |
| same carried members, gained new members (class 2) | 227,581 | 0.2987 |
| split: carried members in >= 2 later clusters (class 3) | 41,283 | 0.0542 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 35,412 | 0.0465 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 59,773 | 0.0785 |
| memo: split and also merged | 24,361 | 0.0320 |
| memo: earlier-release proteins in these clusters | 2,053,075 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 1,937,887 | 0.9439 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 2,053,075 | 1.0000 |
| cluster_uid changed | 902,959 | 0.4398 |
| representative changed | 634,222 | 0.3089 |
| uid same, rep changed | 817 | 0.0004 |
| uid changed, rep same | 269,554 | 0.1313 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 304,276 | 1.0000 |
| best J = 1.0 | 227,581 | 0.7479 |
| best J >= 0.9 | 229,974 | 0.7558 |
| best J >= 0.5 | 279,664 | 0.9191 |
| best J < 0.5 | 24,612 | 0.0809 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.5000 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 279,664 | 0.9191 |

## rL.2 -> rL.3 nr50: earlier size >= 2


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 238,501 | 1.0000 |
| cluster_uid still present in later release | 118,128 | 0.4953 |
| same carried members, none gained (class 1) | 118,128 | 0.4953 |
| same carried members, gained new members (class 2) | 65,841 | 0.2761 |
| split: carried members in >= 2 later clusters (class 3) | 41,283 | 0.1731 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 13,249 | 0.0556 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 37,610 | 0.1577 |
| memo: split and also merged | 24,361 | 0.1021 |
| memo: earlier-release proteins in these clusters | 1,529,743 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 1,414,555 | 0.9247 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 1,529,743 | 1.0000 |
| cluster_uid changed | 719,056 | 0.4701 |
| representative changed | 504,770 | 0.3300 |
| uid same, rep changed | 817 | 0.0005 |
| uid changed, rep same | 215,103 | 0.1406 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 120,373 | 1.0000 |
| best J = 1.0 | 65,841 | 0.5470 |
| best J >= 0.9 | 68,234 | 0.5669 |
| best J >= 0.5 | 107,131 | 0.8900 |
| best J < 0.5 | 13,242 | 0.1100 |
| median best J | 1.0000 |  |
| 10th percentile best J | 0.4333 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 107,131 | 0.8900 |

## rL.2 -> rL.3 nr50: earlier size >= 10


Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)

| quantity | count | fraction |
|---|---|---|
| earlier clusters | 31,552 | 1.0000 |
| cluster_uid still present in later release | 17,182 | 0.5446 |
| same carried members, none gained (class 1) | 17,182 | 0.5446 |
| same carried members, gained new members (class 2) | 3,469 | 0.1099 |
| split: carried members in >= 2 later clusters (class 3) | 9,269 | 0.2938 |
| merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4) | 1,632 | 0.0517 |
| memo: any destination holds members of another earlier cluster (merge, split or not) | 8,991 | 0.2850 |
| memo: split and also merged | 7,359 | 0.2332 |
| memo: earlier-release proteins in these clusters | 836,053 | 1.0000 |
| memo: of them, in their cluster's largest fragment (sum of max k) | 767,458 | 0.9180 |

Carried proteins (measured counts; fractions derived = count / carried proteins in scope)

| quantity | count | fraction |
|---|---|---|
| carried proteins | 836,053 | 1.0000 |
| cluster_uid changed | 370,405 | 0.4430 |
| representative changed | 235,687 | 0.2819 |
| uid same, rep changed | 428 | 0.0005 |
| uid changed, rep same | 135,146 | 0.1616 |

Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))

| quantity | count | fraction of id-changed |
|---|---|---|
| clusters with changed id | 14,370 | 1.0000 |
| best J = 1.0 | 3,469 | 0.2414 |
| best J >= 0.9 | 5,765 | 0.4012 |
| best J >= 0.5 | 11,312 | 0.7872 |
| best J < 0.5 | 3,058 | 0.2128 |
| median best J | 0.8000 |  |
| 10th percentile best J | 0.3636 |  |
| best match is mutual (also that later cluster's best) and J >= 0.5 | 11,312 | 0.7872 |

Later-release clusters (rL.3) with an id not present in rL.2: best-overlap earlier cluster (Jaccard over carried proteins)

| quantity | count | fraction |
|---|---|---|
| clusters in rL.3 | 2,966,200 | 1.0000 |
| id not present in earlier release | 2,508,643 | 0.8457 |
| of them with no carried member (wholly new proteins) | 2,185,833 | 0.7369 |
| of them with carried members | 322,810 | 0.1088 |
|   best J = 1.0 | 227,581 | 0.7050 |
|   best J >= 0.9 | 229,974 | 0.7124 |
|   best J >= 0.5 | 282,109 | 0.8739 |
|   best J < 0.5 | 40,701 | 0.1261 |

(step time 21 s)

## Reading notes

- "Best J >= 0.5" equals "best match is mutual and J >= 0.5" in every table. That is expected: two
  different clusters cannot both have J above 0.5 with the same partner, so a best match at or above
  0.5 is unique on both sides (an exact tie at 0.5 would be the only exception; none showed up in
  the counts).
- In the J tables the "best J = 1.0" row is the class-2 count: those clusters are exactly the ones
  that only gained members, which J over carried proteins cannot see as a change.
- For the later-release direction, 61.7% (NR90) and 61.0% (NR50) of rL.2's clusters and 83.1% / 73.7%
  of rL.3's clusters are wholly new proteins and have no earlier counterpart by construction.
- The three releases are nested (every protein of A is in B), so "carried / earlier" is 1.0000 in
  all four pairs; the restriction to proteins present in both changed nothing here.

## What this does not measure

- Nothing separates real clustering change from MMseqs2 run-to-run variation. No release was
  clustered twice, so the share of "uid same, rep changed" and of splits that is nondeterminism
  (thread count, input order) rather than new data is unknown. `rL.1` was clustered at 16 threads
  and reproduced the pilot exactly once (ladder README, finding 6); that is one agreement.
- Three releases of one nested growth path (phages, then 250, then 2,247 prokaryotes). The later
  step adds 9.3 times more proteins than the earlier one (8,567,609 against 918,525, *derived* as
  differences of the protein counts above), and the earlier release is 1.8 then 5.2 times smaller than
  the later (*derived*), so churn per added protein is not separable from churn per added fraction.
  Whether churn at S2/S6 scale follows either is not measured.
- Only identity of member sets. Nothing about whether clusters' biological content changed.
- No memory or CPU of the measurement itself was recorded; it ran under a 24 GB scope with a 16 GB
  DuckDB limit, 4 threads, and the scope was not exceeded (the unit finished normally).

## Command

Script copied to `/tmp/if2453/churn/churn.py` with `scp`; run detached on the pod:

```
ssh -o BatchMode=yes if2453@wilkcomp02.ccbb.utexas.edu "mkdir -p /tmp/if2453/churn/spill; systemd-run --user --unit=churn-1 --collect -p MemoryMax=24G -p MemorySwapMax=0 nice -n 19 singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/pilot_s1/pydeps --bind /tmp/if2453 --bind /stor/scratch/Wilke/if2453/all-the-data/provisional/rL:/stor/scratch/Wilke/if2453/all-the-data/provisional/rL:ro /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif sh -c 'python /tmp/if2453/churn/churn.py /stor/scratch/Wilke/if2453/all-the-data/provisional/rL /tmp/if2453/churn > /tmp/if2453/churn/log 2>&1'"
```

DuckDB 2.0.0 (from the image's `PYTHONPATH` dependency directory), `memory_limit` 16GB, `threads` 4,
spill under `/tmp/if2453/churn/spill`. Nothing was written under `/stor`. `/tmp/if2453/churn` is
deleted after the run.

## Script (`churn.py`)

```python
"""Cluster-id churn between consecutive releases (read only). Run inside the pod image.

usage: python churn.py <root> <outdir>
File lists come from manifests/release=<R>.json (table cluster_member, partition release == R).
"""
import json
import sys
import time
from pathlib import Path

import duckdb

ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])
PAIRS = [("rL.1", "rL.2"), ("rL.2", "rL.3")]
LEVELS = ["nr90", "nr50"]
THRESH = [1, 2, 10]  # earlier-release cluster size >= t


def files_of(rel):
    m = json.load(open(ROOT / "manifests" / f"release={rel}.json"))
    fs = [f["path"] for f in m["files"] if f["table"] == "cluster_member" and f["partitions"].get("release") == rel]
    assert len(fs) >= 1, rel
    return [str(ROOT / f) for f in fs]


con = duckdb.connect(str(OUT / "churn.duckdb"))
con.execute("SET memory_limit='16GB'")
con.execute("SET threads=4")
con.execute(f"SET temp_directory='{OUT}/spill'")
con.execute("SET preserve_insertion_order=false")

lines = []


def emit(s=""):
    print(s, flush=True)
    lines.append(s)


def table(title, cols, rows):
    emit(f"\n{title}\n")
    emit("| " + " | ".join(cols) + " |")
    emit("|" + "---|" * len(cols))
    for r in rows:
        emit("| " + " | ".join(f"{x:,}" if isinstance(x, int) else (f"{x:.4f}" if isinstance(x, float) else str(x)) for x in r) + " |")


emit("## Schema of cluster_member (rL.1 file)\n")
f1 = files_of("rL.1")
emit("```")
for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet({f1!r})").fetchall():
    emit(f"{r[0]} {r[1]}")
emit("```")
emit("\nFiles per release (from the manifests, partition release == R):\n")
for r in ("rL.1", "rL.2", "rL.3"):
    emit(f"- {r}: " + ", ".join(Path(x).relative_to(ROOT).as_posix() for x in files_of(r)))

for ra, rb in PAIRS:
    for lv in LEVELS:
        t0 = time.time()
        emit(f"\n# {ra} -> {rb}, level {lv}")
        fa, fb = files_of(ra), files_of(rb)
        con.execute("DROP TABLE IF EXISTS a")
        con.execute("DROP TABLE IF EXISTS b")
        con.execute(f"CREATE TABLE a AS SELECT member_uid m, cluster_uid cu, rep_uid rep FROM read_parquet({fa!r}) WHERE level='{lv}'")
        con.execute(f"CREATE TABLE b AS SELECT member_uid m, cluster_uid cu, rep_uid rep FROM read_parquet({fb!r}) WHERE level='{lv}'")
        con.execute("DROP TABLE IF EXISTS csa")
        con.execute("DROP TABLE IF EXISTS csb")
        con.execute("CREATE TABLE csa AS SELECT cu, count(*)::BIGINT s FROM a GROUP BY cu")
        con.execute("CREATE TABLE csb AS SELECT cu, count(*)::BIGINT s FROM b GROUP BY cu")
        con.execute("DROP TABLE IF EXISTS j")
        con.execute("CREATE TABLE j AS SELECT a.m, a.cu cua, a.rep repa, b.cu cub, b.rep repb FROM a JOIN b USING (m)")
        na, nb, nj = (con.execute(f"SELECT count(*) FROM {t}").fetchone()[0] for t in ("a", "b", "j"))
        ca, cb = (con.execute(f"SELECT count(*) FROM {t}").fetchone()[0] for t in ("csa", "csb"))
        dup = con.execute("SELECT (SELECT count(*) FROM a) - (SELECT count(DISTINCT m) FROM a), (SELECT count(*) FROM b) - (SELECT count(DISTINCT m) FROM b)").fetchone()
        table("Sizes (measured)", ["quantity", "value"], [
            [f"proteins in {ra}", na], [f"proteins in {rb}", nb], ["carried (in both)", nj],
            ["carried / earlier", nj / na], [f"clusters in {ra}", ca], [f"clusters in {rb}", cb],
            ["duplicate member rows (a, b)", f"{dup[0]}, {dup[1]}"]])

        # pair table: old cluster x new cluster, k carried members in common
        con.execute("DROP TABLE IF EXISTS p")
        con.execute("CREATE TABLE p AS SELECT cua, cub, count(*)::BIGINT k FROM j GROUP BY cua, cub")
        # carried members per old cluster (nc) and per new cluster (mb)
        con.execute("DROP TABLE IF EXISTS pa")
        con.execute("CREATE TABLE pa AS SELECT cua, sum(k)::BIGINT nc, count(*)::BIGINT ndest, max(k) maxk FROM p GROUP BY cua")
        con.execute("DROP TABLE IF EXISTS pb")
        con.execute("CREATE TABLE pb AS SELECT cub, sum(k)::BIGINT mb FROM p GROUP BY cub")
        # per pair: Jaccard over carried proteins
        con.execute("DROP TABLE IF EXISTS pj")
        con.execute("""CREATE TABLE pj AS SELECT p.cua, p.cub, p.k, pa.nc, pb.mb, p.k::DOUBLE/(pa.nc+pb.mb-p.k) jac
                       FROM p JOIN pa USING (cua) JOIN pb USING (cub)""")
        con.execute("DROP TABLE IF EXISTS single")
        con.execute("""CREATE TABLE single AS
          SELECT p.cua, pb.mb mb1, csb.s sb1 FROM p JOIN pa USING (cua) JOIN pb USING (cub) JOIN csb ON csb.cu=p.cub
          WHERE pa.ndest=1""")
        con.execute("DROP TABLE IF EXISTS oc")
        con.execute("""CREATE TABLE oc AS
          SELECT pa.cua, csa.s sa, pa.nc, pa.ndest, pa.maxk, (k.cu IS NOT NULL) id_kept, single.mb1, single.sb1,
                 -- any destination holding carried members from another old cluster
                 coalesce(f.anyforeign, false) anyforeign
          FROM pa JOIN csa ON csa.cu=pa.cua
          LEFT JOIN csb k ON k.cu=pa.cua
          LEFT JOIN single ON single.cua=pa.cua
          LEFT JOIN (SELECT cua, bool_or(mb>k) anyforeign FROM pj GROUP BY cua) f ON f.cua=pa.cua""")
        con.execute("DROP TABLE IF EXISTS ocl")
        con.execute("""CREATE TABLE ocl AS SELECT *,
          CASE WHEN ndest>=2 THEN '3 split'
               WHEN mb1>nc THEN '4 merged'
               WHEN sb1=nc THEN '1 same members'
               ELSE '2 same carried, gained new' END cls FROM oc""")

        for t in THRESH:
            lab = {1: "all clusters", 2: "earlier size >= 2", 10: "earlier size >= 10"}[t]
            r = con.execute(f"""SELECT count(*), sum(id_kept::INT), sum(sa), sum(maxk), sum(nc),
                   sum((cls='1 same members')::INT), sum((cls='2 same carried, gained new')::INT),
                   sum((cls='3 split')::INT), sum((cls='4 merged')::INT), sum(anyforeign::INT),
                   sum((cls='3 split' AND anyforeign)::INT)
                   FROM ocl WHERE sa>={t}""").fetchone()
            n, kept, sa_sum, maxk_sum, nc_sum, same, gained, split, merged, anyf, splitf = r
            emit(f"\n## {ra} -> {rb} {lv}: {lab}\n")
            table("Clusters of the earlier release (measured counts; fractions derived = count / earlier clusters)",
                  ["quantity", "count", "fraction"], [
                      ["earlier clusters", n, 1.0],
                      ["cluster_uid still present in later release", kept, kept / n],
                      ["same carried members, none gained (class 1)", same, same / n],
                      ["same carried members, gained new members (class 2)", gained, gained / n],
                      ["split: carried members in >= 2 later clusters (class 3)", split, split / n],
                      ["merged: one later cluster that also holds carried members of another earlier cluster, not split (class 4)", merged, merged / n],
                      ["memo: any destination holds members of another earlier cluster (merge, split or not)", anyf, anyf / n],
                      ["memo: split and also merged", splitf, splitf / n],
                      ["memo: earlier-release proteins in these clusters", sa_sum, 1.0],
                      ["memo: of them, in their cluster's largest fragment (sum of max k)", maxk_sum, maxk_sum / sa_sum]])
            pr = con.execute(f"""SELECT count(*), sum((j.cua<>j.cub)::INT), sum((j.repa<>j.repb)::INT),
                   sum((j.cua=j.cub AND j.repa<>j.repb)::INT), sum((j.cua<>j.cub AND j.repa=j.repb)::INT)
                   FROM j JOIN csa ON csa.cu=j.cua WHERE csa.s>={t}""").fetchone()
            n, ch, rch, both_a, both_b = pr
            table("Carried proteins (measured counts; fractions derived = count / carried proteins in scope)",
                  ["quantity", "count", "fraction"], [
                      ["carried proteins", n, 1.0],
                      ["cluster_uid changed", ch, ch / n],
                      ["representative changed", rch, rch / n],
                      ["uid same, rep changed", both_a, both_a / n],
                      ["uid changed, rep same", both_b, both_b / n]])
            # Jaccard of best match, for clusters whose id changed
            jr = con.execute(f"""WITH bm AS (
                   SELECT oc.cua, max(pj.jac) best FROM oc JOIN pj USING (cua)
                   WHERE NOT oc.id_kept AND oc.sa>={t} GROUP BY oc.cua)
                SELECT count(*), sum((best>=1.0)::INT), sum((best>=0.9)::INT), sum((best>=0.5)::INT), sum((best<0.5)::INT),
                       median(best), quantile_cont(best,0.1) FROM bm""").fetchone()
            n2, j1, j9, j5, jl, med, q10 = [x if x is not None else 0 for x in jr]
            # mutual best (one-to-one) pairs with J >= 0.5
            mr = con.execute(f"""WITH mx_a AS (SELECT cua, max(jac) m FROM pj GROUP BY cua),
                 mx_b AS (SELECT cub, max(jac) m FROM pj GROUP BY cub)
                SELECT count(DISTINCT pj.cua) FROM pj JOIN mx_a USING (cua) JOIN mx_b USING (cub) JOIN oc USING (cua)
                WHERE NOT oc.id_kept AND oc.sa>={t} AND pj.jac=mx_a.m AND pj.jac=mx_b.m AND pj.jac>=0.5""").fetchone()[0]
            table("Earlier clusters whose id changed: Jaccard of the best-overlap later cluster over carried proteins (J = k / (nc + mb - k))",
                  ["quantity", "count", "fraction of id-changed"], [
                      ["clusters with changed id", n2, 1.0],
                      ["best J = 1.0", j1, j1 / n2 if n2 else 0.0],
                      ["best J >= 0.9", j9, j9 / n2 if n2 else 0.0],
                      ["best J >= 0.5", j5, j5 / n2 if n2 else 0.0],
                      ["best J < 0.5", jl, jl / n2 if n2 else 0.0],
                      ["median best J", float(med), ""],
                      ["10th percentile best J", float(q10), ""],
                      ["best match is mutual (also that later cluster's best) and J >= 0.5", mr, mr / n2 if n2 else 0.0]])

        # direction new -> old (all sizes): later clusters whose id is new
        nr = con.execute("""WITH bm AS (SELECT pb.cub, max(pj.jac) best FROM pb JOIN pj USING (cub)
                                         LEFT JOIN csa ON csa.cu=pb.cub WHERE csa.cu IS NULL GROUP BY pb.cub)
                            SELECT count(*), sum((best>=1.0)::INT), sum((best>=0.9)::INT), sum((best>=0.5)::INT), sum((best<0.5)::INT)
                            FROM bm""").fetchone()
        nnew = con.execute("SELECT count(*) FROM csb LEFT JOIN csa USING (cu) WHERE csa.cu IS NULL").fetchone()[0]
        nwhole = con.execute("SELECT count(*) FROM csb LEFT JOIN pb ON pb.cub=csb.cu WHERE pb.cub IS NULL").fetchone()[0]
        n2, j1, j9, j5, jl = [x if x is not None else 0 for x in nr]
        table(f"Later-release clusters ({rb}) with an id not present in {ra}: best-overlap earlier cluster (Jaccard over carried proteins)",
              ["quantity", "count", "fraction"], [
                  [f"clusters in {rb}", cb, 1.0],
                  ["id not present in earlier release", nnew, nnew / cb],
                  ["of them with no carried member (wholly new proteins)", nwhole, nwhole / cb],
                  ["of them with carried members", n2, n2 / cb],
                  ["  best J = 1.0", j1, j1 / n2 if n2 else 0.0],
                  ["  best J >= 0.9", j9, j9 / n2 if n2 else 0.0],
                  ["  best J >= 0.5", j5, j5 / n2 if n2 else 0.0],
                  ["  best J < 0.5", jl, jl / n2 if n2 else 0.0]])
        emit(f"\n(step time {time.time() - t0:.0f} s)")
        for t in ("a", "b", "j", "p", "pa", "pb", "pj", "oc", "ocl", "single", "csa", "csb"):
            con.execute(f"DROP TABLE IF EXISTS {t}")

(OUT / "result.md").write_text("\n".join(lines) + "\n")
(OUT / "done").write_text("ok\n")
```

## What the data shows

Ids churn little while a release grows by a factor of about 1.8 and a lot when it grows by a factor of about 5: between rL.2 and rL.3, 12.4% of NR90 clusters and 39.9% of NR50 clusters lose their id (11.4% and 44.0% of carried proteins change cluster_uid), against 0.2% and 2.4% of clusters between rL.1 and rL.2. Members mostly stay together even when ids change: 99.6% (NR90) and 94.4% (NR50) of earlier proteins remain in the largest fragment of their old cluster in rL.2 to rL.3, and most id changes (92.5% NR90, 74.8% NR50) are clusters that only gained members, while splits and merges together affect 0.9% (NR90) and 10.1% (NR50) of earlier clusters (*derived*: (class 3 + class 4) / earlier clusters; 29% of NR50 clusters of at least 10 members are split). An overlap-based mapping of new ids onto old ones looks workable by this measure: for 98.7% (NR90) and 91.9% (NR50) of id-changed earlier clusters the best-overlap later cluster has J of at least 0.5 and is mutually best, but the J here ignores new members, covers only two steps, and the data does not show what should happen to the second fragment of a split or to the two sides of a merge (about 8% of NR50 id-changed clusters, and 21% of those with at least 10 members, have best J below 0.5).
