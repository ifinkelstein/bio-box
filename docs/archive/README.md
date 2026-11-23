# Archived Initial Specifications and Pre-Build Audits

These documents are the **October–November 2026 pre-build specifications, stage plans, prior-art survey, and external audit rounds** that launched `bio-box`.

> **Authority rule ([`../../AGENTS.md`](../../AGENTS.md) §1.2):** Several sections of these initial plans (including `mmseqs clusterupdate`, materialized `slot_pairs`, ordinal keys, and full-Rfam whole-genome scans) were superseded after empirical stress testing. **[`../../planning/DECISIONS.md`](../../planning/DECISIONS.md)** and **[`../../planning/CATALOG_GOAL.md`](../../planning/CATALOG_GOAL.md)** always outrank these files where they differ. For a full map of live documentation, start at **[`../../DOCS.md`](../../DOCS.md)**.

## Files

| Path | Role |
|---|---|
| [`00-spec.md`](00-spec.md) | Original system specification; §3 invariants remain authoritative unless amended in [`../../planning/DECISIONS.md`](../../planning/DECISIONS.md) |
| [`PLAN.md`](PLAN.md) | Original pipeline and L0–L5 schema design (token bit layout §4.1 and tape geometry §4.2 remain referenced by R0/R1/R2 tests) |
| [`STAGES.md`](STAGES.md) | Original stage ladder `S0`..`S7+` and the 2026-11-11 status mapping |
| [`DEPLOYMENT-AND-MALLEABILITY.md`](DEPLOYMENT-AND-MALLEABILITY.md) | Original container, `wilkcomp02`, TACC, and bare-metal portability plan |
| [`PRIOR-ART.md`](PRIOR-ART.md) | Survey of prior microbial genomics, synteny, and homology platforms |
| [`audits/SYNTHESIS.md`](audits/SYNTHESIS.md), [`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md) | Synthesis of the two pre-build adversarial review rounds (`D-01`..`D-28`) |
| [`audits/critique_1_biology.md`](audits/critique_1_biology.md) .. [`audits/critique_8_data_quality.md`](audits/critique_8_data_quality.md) | Eight domain-specific pre-build critiques (biology, systems, ML/stats, bioinformatics, data science, architecture, database, data quality) |
