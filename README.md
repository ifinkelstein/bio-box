# bio-box (all-the-data)

Staged, content-addressed corpus of prokaryotic and viral genomes for neighborhood, fusion and search workflows. Planning docs are authoritative; code is starting at S0.

| Path | What |
|---|---|
| `00-spec.md`, `PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md` | specification, plan, stage ladder, deployment |
| `planning/` | data layer and tests, sizing, preflight gaps, search plan, test architecture |
| `src/biobox/registry/` | schema registry (Arrow schemas, keys, enums) |
| `tests/registry/` | operation registry, metamorphic relations, failure modes, pending tables |
| `tests/r0`, `tests/r1` | ring R0 (registry lint) and R1 (unit) tests |
| `audits/`, `presentations/` | audit record and the group-meeting deck |

Run the tests (laptop or CI host; ring R0–R1):

```
uv sync
uv run pytest
```

On wilkcomp02 the checkout lives at `~/bio-box`; data lives under `/stor/scratch/Wilke/if2453/all-the-data` and `/stor/work/Wilke/if2453/all-the-data` (`DEPLOYMENT-AND-MALLEABILITY.md` §1).
