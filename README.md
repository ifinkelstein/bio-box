# bio-box (all-the-data)

Staged, content-addressed corpus of prokaryotic and viral genomes for neighborhood, fusion and search workflows. Planning docs are authoritative; code is starting at S0.

| Path | What |
|---|---|
| `00-spec.md`, `PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md` | specification, plan, stage ladder, deployment |
| `planning/` | data layer and tests, sizing, preflight gaps, search plan, test architecture |
| `src/biobox/registry/` | schema registry (Arrow schemas, keys, enums, file locations) |
| `src/biobox/l0seq/`, `l1/`, `l2/` | sequence store (bgzip shards, DEC-17), gene calling, domain scans, clustering (S0) |
| `src/biobox/catalog.py` | release manifest, generated read-only DuckDB catalog, sequence macros |
| `tests/registry/` | operation registry, metamorphic relations, failure modes, pending tables |
| `tests/r0`, `tests/r1`, `tests/r2`, `tests/r3` | rings R0 (registry lint), R1 (unit), R2 (micro-S0 fixture), R3 (golden set; needs `BIOBOX_GOLDEN_DIR`) |
| `audits/`, `presentations/` | audit record and the group-meeting deck |

Run the tests (laptop or CI host; rings R0–R2):

```
uv sync
uv run pytest
```

On wilkcomp02 the checkout lives at `~/bio-box`; data lives under `/stor/scratch/Wilke/if2453/all-the-data` and `/stor/work/Wilke/if2453/all-the-data` (`DEPLOYMENT-AND-MALLEABILITY.md` §1).
