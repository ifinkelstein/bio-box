# presentations/

`all-the-data-deck.html` is a 46-slide group-meeting / pitch deck for the staged bio-box plan, built with reveal.js 5.1.0 (loaded from cdnjs) and pre-rendered mermaid diagrams.

## Files

| Path | What |
|---|---|
| `all-the-data-deck.html` | The deck. Open in a browser. Arrow keys move between slides; `s` opens speaker notes; `?print-pdf` in the URL gives the print layout (then Chrome → Print → Save as PDF, landscape, no margins). |
| `diagrams/d*.mmd` | Mermaid sources: layers, lineage, per-genome pipeline, barriers, containers, run-step, catalog, ladder, statistics, accretion. |
| `diagrams/d*.svg` | Rendered diagrams referenced by the deck. |
| `diagrams/mermaid-config.json`, `puppeteer-config.json` | Theme, plus the path to the system Chrome that mermaid-cli uses. |
| `render.sh` | Re-renders every `.mmd` to `.svg` with mermaid-cli 11.12; exits non-zero if any diagram fails to parse. |
| `all-the-data-deck.pptx`, `diagrams/png/d*.png` | PowerPoint version (16:9, 51 slides): the 46 HTML slides in order with titles, bullets, tables, source lines and speaker notes, plus a 5-slide section "Update, 12 November 2026" before the appendix, built only from `planning/tasks/STATUS.md` and `planning/DECISIONS.md` DEC-24 to DEC-29. The PNGs are the diagrams rendered at 3× from the SVGs. Built by `build_pptx.py`. |
| `build_pptx.py` | Builds the `.pptx` from the HTML (the update slides are in the script as `UPDATE_HTML`) and renders the diagram PNGs. |

## Rebuilding the PowerPoint deck

From the repository root (Python packages come through `uv run --with`; nothing is added to `pyproject.toml` or `uv.lock`):

```sh
# Diagrams: SVG -> PNG with the installed Google Chrome (the SVGs have HTML labels, which
# rsvg-convert and Inkscape do not draw). Only needed when a diagram changes.
uv run --no-project --with playwright==1.48.0 python presentations/build_pptx.py --diagrams
# Deck
uv run --no-project --with python-pptx==1.0.2 --with beautifulsoup4==4.12.3 python presentations/build_pptx.py
```

Text heights are estimated from the Arial metrics on the host; each slide takes the largest sizes that fit (body 18 to 16 pt, tables 18 to 12 pt, notes and source lines 14 and 12 pt). Check a rebuild by rendering it (`soffice --headless --convert-to pdf`, then `pdftoppm -png`) and looking at every slide.

## How this brief was rewritten

**What a top professional would do that was not asked:**
- Build the deck around the decision the meeting has to make (approve preflight and S0, the pod request, owners), with the ask on slide 2.
- Give every slide a one-sentence claim as its title.
- Put detail in an appendix.
- Mark every number as measured, derived or estimate.
- Draw diagrams that show the actual mechanism, not generic boxes.
- Include the pre-mortem and kill criteria, not just the happy path.

**Rawest material:** the plan files themselves, not summaries of them:
- `STAGES.md`, `planning/*.md`, `PLAN.md` v3, `DEPLOYMENT-AND-MALLEABILITY.md` v3, `audits/SYNTHESIS_2.md`
- the wilkcomp02 envelope notes
- the pre-mortem in `planning/preflight_gaps.md`

**Methods and tools:**
- Assertion-evidence slides (pyramid principle).
- Mermaid flowcharts and sequence diagrams, rendered offline with mermaid-cli so they work without network and can be parse-checked.
- reveal.js with speaker notes.
- Headless-Chrome screenshots of every slide to catch overflow.

**What shows it is right, other than us:**
- `render.sh` parses and renders every diagram without error.
- Every slide was screenshotted and checked for overflow (done 2026-10-28).
- Every number on a slide traces to a named section of a plan file; each slide carries a source line.
- Still to do: a dry run with a lab member who has not read the plan, who should be able to restate the ask and the stage ladder afterwards.

**Rewritten brief:**
> Before executing anything, list remaining gaps (pre-mortem). Then build a group-meeting deck in `presentations/`. It leads with the decision needed and uses claim-style titles and terse bullets. Mermaid diagrams show the data layers, the per-genome pipeline and global barriers, container composition and the tools in it, the wilkcomp02 run-step wrapper, the DuckDB catalog, and accretion equivalence. One card per stage gives what it adds, the capability, the science check, the gate and estimates. It also covers portability, risks, kill criteria and the first-10 questions, with an appendix. Every number traces to a plan file, and every diagram renders from source.

## Known limits

- reveal.js loads from cdnjs, so the deck needs network on first open. The diagrams are local SVGs.
- Headless `--print-to-pdf` produced an empty file; export to PDF manually as described above.
- Numbers are estimates until S0 measures real rates; the deck says so on the slides that use them.
- In the `.pptx`, the run-step sequence diagram (slide 17) shares its slide with eight boxes and is about 2.2 in tall, so its labels are small; the full-size render is `diagrams/png/d06_runstep.png`. Only slides 1, 2 and 42 of the HTML deck have speaker notes, so only those (and the update slides) have notes in the `.pptx`.
- The deck predates the 2026-10-29 search round: it has no stage S2C (complete prokaryotic genomes), no search tiers ST0–ST5 and no search pack. See `planning/search/SEARCH_PLAN.md`.
