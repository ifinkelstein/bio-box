"""Card C-48b: JGI GOLD + IMG/M environmental metadata crosswalk and DuckDB sidecar tables.

Constructed synthetic fixtures (no external network, no top-level repo files outside the image):
- Streaming .xlsx reader (inline strings & sharedStrings) with PII/contact column guard
- GOLD Gs/Gb/Gp/Ga <-> IMG taxon OID <-> SRA run / BioSample linking & SPIRE/MGnify/INSDC crosswalk
- GOLD 5-level ecosystem harmonization to canonical top-level habitats
- JGI Data Utilization Status (Unrestricted vs Restricted) & redistributable flag (DEC-51)
- Deterministic byte-identical Parquet output across shuffled input orders
- Bounded RSS on 100x synthetic inputs in an isolated subprocess
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

import duckdb
import pyarrow.parquet as pq
import pytest

pytestmark = pytest.mark.covers("C-48b")

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "stress" / "c48b_env_crosswalk.py"


def _mod():
    spec = importlib.util.spec_from_file_location("c48b_env_crosswalk", SCRIPT)
    assert spec is not None and spec.loader is not None
    m = importlib.util.module_from_spec(spec)
    sys.modules["c48b_env_crosswalk"] = m
    spec.loader.exec_module(m)
    return m


def _col_letter(idx: int) -> str:
    out = []
    n = idx + 1
    while n > 0:
        n, rem = divmod(n - 1, 26)
        out.append(chr(65 + rem))
    return "".join(reversed(out))


def _write_xlsx(path: Path, sheets: list[tuple[str, list[list[str]]]], *, use_shared_strings: bool = False) -> None:
    """Write a minimal valid OOXML .xlsx file for testing (inlineStr or sharedStrings)."""
    shared: list[str] = []
    shared_idx: dict[str, int] = {}

    def _sid(val: str) -> int:
        if val not in shared_idx:
            shared_idx[val] = len(shared)
            shared.append(val)
        return shared_idx[val]

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        wb_sheets = []
        rels = []
        for s_i, (name, rows) in enumerate(sheets, start=1):
            wb_sheets.append(
                f'<sheet name="{escape(name)}" sheetId="{s_i}" r:id="rId{s_i}"/>'
            )
            rels.append(
                f'<Relationship Id="rId{s_i}" '
                f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                f'Target="worksheets/sheet{s_i}.xml"/>'
            )
            xml_rows = []
            for r_i, row in enumerate(rows, start=1):
                cells = []
                for c_i, val in enumerate(row):
                    ref = f"{_col_letter(c_i)}{r_i}"
                    if val == "":
                        continue
                    if use_shared_strings:
                        cells.append(f'<c r="{ref}" t="s"><v>{_sid(val)}</v></c>')
                    else:
                        cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(val)}</t></is></c>')
                xml_rows.append(f'<row r="{r_i}">{"".join(cells)}</row>')
            sheet_xml = (
                '<?xml version="1.0" encoding="UTF-8"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>'
            )
            zf.writestr(f"xl/worksheets/sheet{s_i}.xml", sheet_xml)

        wb_xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f'<sheets>{"".join(wb_sheets)}</sheets></workbook>'
        )
        rels_xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{"".join(rels)}</Relationships>'
        )
        zf.writestr("xl/workbook.xml", wb_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", rels_xml)
        if use_shared_strings:
            si_xml = "".join(f"<si><t>{escape(s)}</t></si>" for s in shared)
            zf.writestr(
                "xl/sharedStrings.xml",
                '<?xml version="1.0" encoding="UTF-8"?>'
                f'<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">{si_xml}</sst>',
            )


def _write_tsv_gz(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def _write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


IMG_M_HEADER = [
    "NCBI Domain",
    "IMG Genome ID",
    "GOLD Analysis Project Type",
    "GOLD Analysis Project ID",
    "GOLD Sequencing Project ID",
    "GOLD Study ID",
    "Sequencing Center",
    "NCBI Assembly Accession",
    "NCBI Bioproject Accession",
    "NCBI Biosample Accession",
    "NCBI GenBank ID",
    "Ecosystem",
    "Ecosystem Category",
    "Ecosystem Type",
    "Ecosystem Subtype",
    "Specific Ecosystem",
    "Habitat",
    "pH",
    "Sample Collection Temperature",
    "Temperature Range",
    "Salinity",
    "Salinity Concentration",
    "Oxygen Requirement",
    "Oxygen Concentration",
    "Depth In Meters",
    "Elevation In Meters",
    "Subsurface In Meters",
    "Latitude",
    "Longitude",
    "Longhurst Code",
    "Host Name",
    "Sample Body Site",
    "Sample Body Subsite",
    "Geographic Location",
    "Genome Size * assembled",
    "SRA ID",
    "SRA Run",
    "JGI Data Utilization Status",
]


def _make_synthetic_bundle(base: Path, *, reverse_order: bool = False) -> dict[str, Path]:
    """Create a complete synthetic GOLD + IMG/M + SPIRE + MGnify fixture with known truth."""
    gold_dir = base / "gold"
    img_dir = base / "img_m"
    spire_dir = base / "spire_study_tsv"
    gold_dir.mkdir(parents=True, exist_ok=True)
    img_dir.mkdir(parents=True, exist_ok=True)
    spire_dir.mkdir(parents=True, exist_ok=True)

    # GOLD Biosample
    bs_rows = [
        ["Gb0001", "Marine sample 1", "408172", "marine metagenome", "Pacific", "2024-01-01",
         "North Pacific", "35.5", "-140.2", "101", "Environmental", "Aquatic", "Marine", "Pelagic", "Photic zone"],
        ["Gb0002", "Soil sample 2", "410658", "soil metagenome", "Forest", "2024-02-01",
         "Oregon", "44.2", "-122.1", "201", "Environmental", "Terrestrial", "Soil", "Forest soil", "Rhizosphere"],
        ["Gb0003", "Gut sample 3", "408170", "human gut metagenome", "Gut", "2024-03-01",
         "USA", "30.3", "-97.7", "301", "Host-associated", "Human", "Digestive system", "Large intestine", "Fecal"],
        ["Gb0004", "WWTP sample 4", "527639", "wastewater metagenome", "WWTP", "2024-04-01",
         "Texas", "30.2", "-97.8", "401", "Engineered", "Wastewater", "Activated sludge", "Municipal", "Unclassified"],
    ]
    # GOLD Sequencing Project (includes PII column PROJECT CONTACTS that must be stripped)
    sp_rows = [
        ["Gp0001", "Marine SP", "WGS", "PRJNA1", "SAMN00000001", "SRX1001", "JGI", "Gs0001", "", "Gb0001",
         "secret@example.org", "Unrestricted"],
        ["Gp0002", "Soil SP", "WGS", "PRJNA2", "SAMN00000002", "SRX1002", "JGI", "Gs0002", "", "Gb0002",
         "secret@example.org", "Unrestricted"],
        ["Gp0003", "Gut SP", "WGS", "PRJNA3", "SAMN00000003", "SRX1003", "JGI", "Gs0003", "", "Gb0003",
         "secret@example.org", "Unrestricted"],
        ["Gp0004", "WWTP SP", "WGS", "PRJNA4", "SAMN00000004", "SRX1004", "JGI", "Gs0004", "", "Gb0004",
         "secret@example.org", "Restricted"],
    ]
    # GOLD Analysis Project (includes PII column AP CONTACTS that must be stripped)
    ap_rows = [
        ["Ga0001", "Marine AP", "Metagenome Analysis", "3300000001", "", "SRR100001", "Gs0001", "", "Gp0001",
         "pi@example.org", "Unrestricted", "50000", "25000000"],
        ["Ga0002", "Soil AP", "Metagenome Analysis", "3300000002", "", "SRR100002", "Gs0002", "", "Gp0002",
         "pi@example.org", "Unrestricted", "80000", "40000000"],
        ["Ga0003", "Gut AP", "Metagenome Analysis", "3300000003", "", "SRR100003", "Gs0003", "", "Gp0003",
         "pi@example.org", "Unrestricted", "60000", "30000000"],
        ["Ga0004", "WWTP AP", "Metagenome Analysis", "3300000004", "", "SRR100004", "Gs0004", "", "Gp0004",
         "pi@example.org", "Restricted", "90000", "50000000"],
        ["Ga0005", "Isolate AP", "Genome Analysis", "2500000005", "GCA_000123456.1", "", "Gs0005", "Go0005", "Gp0005",
         "pi@example.org", "Unrestricted", "4000", "4200000"],
    ]
    org_rows = [
        ["Go0005", "Escherichia coli K-12", "511145", "Human gut", "USA", "30.3", "-97.7", "301",
         "Host-associated", "Human", "Digestive system", "Large intestine", "Fecal",
         "Homo sapiens", "Gut", "Fecal", "Non-saline", "", "Facultative", "7.0", "Mesophile", "37"],
    ]
    sra_biome_rows = [
        ["Ga0001", "Marine AP", "Metagenome Analysis", "3300000001", "SRR100001", "Gs0001", "Gp0001",
         "101", "Environmental", "Aquatic", "Marine", "Pelagic", "Photic zone"],
    ]
    eco_paths_rows = [
        ["101", "Environmental", "Aquatic", "Marine", "Pelagic", "Photic zone"],
        ["201", "Environmental", "Terrestrial", "Soil", "Forest soil", "Rhizosphere"],
        ["301", "Host-associated", "Human", "Digestive system", "Large intestine", "Fecal"],
        ["401", "Engineered", "Wastewater", "Activated sludge", "Municipal", "Unclassified"],
    ]

    # IMG/M microbiome rows:
    # 3300000001: Marine, Unrestricted, SRR100001 / SAMN00000001 -> COVERED by SPIRE (SAMN00000001)
    # 3300000002: Soil, Unrestricted, SRR100002 / SAMN00000002 -> NOVEL (not in SPIRE or MGnify) -> TAKE
    # 3300000003: Gut, Unrestricted, SRR100003 / SAMN00000003 -> COVERED by MGnify (SRR100003)
    # 3300000004: WWTP, Restricted, SRR100004 / SAMN00000004 -> RESTRICTED (excluded from unrestricted take-list)
    # 3300000006: Air, Unrestricted, no SRA/BioSample -> NOVEL UNLINKED -> TAKE_UNLINKED
    img_rows = [
        ["*Microbiome", "3300000001", "Metagenome Analysis", "Ga0001", "Gp0001", "Gs0001", "JGI",
         "", "PRJNA1", "SAMN00000001", "", "Environmental", "Aquatic", "Marine", "Pelagic", "Photic zone",
         "Marine pelagic", "8.1", "18.5", "Mesophile", "Marine", "35 PSU", "Aerobe", "6.2 mg/L",
         "15", "", "", "35.5", "-140.2", "NPTE", "", "", "", "North Pacific", "25000000", "SRX1001", "SRR100001",
         "Unrestricted"],
        ["*Microbiome", "3300000002", "Metagenome Analysis", "Ga0002", "Gp0002", "Gs0002", "JGI",
         "", "PRJNA2", "SAMN00000002", "", "Environmental", "Terrestrial", "Soil", "Forest soil", "Rhizosphere",
         "Forest soil", "5.8", "12.0", "Mesophile", "Non-saline", "", "Aerobe", "",
         "", "450", "0.2", "44.2", "-122.1", "", "", "", "", "Oregon", "40000000", "SRX1002", "",
         "Unrestricted"],
        ["*Microbiome", "3300000003", "Metagenome Analysis", "Ga0003", "Gp0003", "Gs0003", "JGI",
         "", "PRJNA3", "SAMN00000003", "", "Host-associated", "Human", "Digestive system", "Large intestine", "Fecal",
         "Human gut", "6.8", "37.0", "Mesophile", "", "", "Anaerobe", "",
         "", "", "", "30.3", "-97.7", "", "Homo sapiens", "Digestive system", "Fecal", "USA", "30000000",
         "SRX1003", "SRR100003", "Unrestricted"],
        ["*Microbiome", "3300000004", "Metagenome Analysis", "Ga0004", "Gp0004", "Gs0004", "JGI",
         "", "PRJNA4", "SAMN00000004", "", "Engineered", "Wastewater", "Activated sludge", "Municipal", "Unclassified",
         "Wastewater", "7.4", "22.0", "Mesophile", "", "", "Aerobe", "",
         "", "", "", "30.2", "-97.8", "", "", "", "", "Texas", "50000000", "SRX1004", "SRR100004",
         "Restricted"],
        ["*Microbiome", "3300000006", "Metagenome Analysis", "Ga0006", "Gp0006", "Gs0006", "JGI",
         "", "", "", "", "Environmental", "Air", "Outdoor Air", "Unclassified", "Unclassified",
         "Outdoor air", "", "20.0", "Mesophile", "", "", "Aerobe", "",
         "", "100", "", "40.0", "-105.0", "", "", "", "", "Colorado", "12000000", "", "",
         "Unrestricted"],
    ]

    if reverse_order:
        bs_rows = list(reversed(bs_rows))
        sp_rows = list(reversed(sp_rows))
        ap_rows = list(reversed(ap_rows))
        img_rows = list(reversed(img_rows))

    _write_xlsx(
        gold_dir / "goldData.xlsx",
        [
            ("Biosample", [
                ["BIOSAMPLE GOLD ID", "BIOSAMPLE NAME", "BIOSAMPLE NCBI TAX ID", "BIOSAMPLE NCBI TAX NAME",
                 "BIOSAMPLE SAMPLE COLLECTION SITE", "BIOSAMPLE SAMPLE COLLECTION DATE",
                 "BIOSAMPLE GEOGRAPHIC LOCATION", "BIOSAMPLE LATITUDE", "BIOSAMPLE LONGITUDE",
                 "BIOSAMPLE ECOSYSTEM PATH ID", "BIOSAMPLE ECOSYSTEM", "BIOSAMPLE ECOSYSTEM CATEGORY",
                 "BIOSAMPLE ECOSYSTEM TYPE", "BIOSAMPLE ECOSYSTEM SUBTYPE", "BIOSAMPLE SPECIFIC ECOSYSTEM"],
                *bs_rows,
            ]),
            ("Organism", [
                ["ORGANISM GOLD ID", "ORGANISM NAME", "ORGANISM NCBI TAX ID", "ORGANISM HABITAT",
                 "ORGANISM GEOGRAPHIC LOCATION", "ORGANISM LATITUDE", "ORGANISM LONGITUDE",
                 "ORGANISM ECOSYSTEM PATH ID", "ORGANISM ECOSYSTEM", "ORGANISM ECOSYSTEM CATEGORY",
                 "ORGANISM ECOSYSTEM TYPE", "ORGANISM ECOSYSTEM SUBTYPE", "ORGANISM SPECIFIC ECOSYSTEM",
                 "ORGANISM ISOLATION HOST NAME", "ORGANISM HOST BODY SITE", "ORGANISM HOST BODY SUBSITE",
                 "ORGANISM SALINITY", "ORGANISM SALINITY CONCENTRATION", "ORGANISM OXYGEN REQUIREMENT",
                 "ORGANISM PH", "ORGANISM TEMPERATURE RANGE", "ORGANISM SAMPLE COLLECTION TEMPERATURE"],
                *org_rows,
            ]),
            ("Sequencing Project", [
                ["PROJECT GOLD ID", "PROJECT NAME", "SEQUENCING STRATEGY", "NCBI BIOPROJECT ACCESSION",
                 "NCBI BIOSAMPLE ACCESSION", "SRA EXPERIMENT IDS", "SEQUENCING CENTERS",
                 "STUDY GOLD ID", "ORGANISM GOLD ID", "BIOSAMPLE GOLD ID", "PROJECT CONTACTS",
                 "PROJECT JGI DATA UTILIZATION STATUS\n (For projects sequenced at JGI only)"],
                *sp_rows,
            ]),
            ("Analysis Project", [
                ["AP GOLD ID", "AP NAME", "AP TYPE", "AP IMG TAXON ID", "AP GENBANK", "AP NCBI SRA IDS",
                 "AP STUDY GOLD ID", "AP ORGANISM GOLD ID", "AP PROJECT GOLD IDS", "AP CONTACTS",
                 "AP JGI DATA UTILIZATION STATUS\n (For projects sequenced at JGI only)",
                 "AP GENE COUNT", "AP ESTIMATED SIZE"],
                *ap_rows,
            ]),
        ],
    )
    _write_xlsx(
        gold_dir / "goldSraBiomeImg.xlsx",
        [
            ("SRA biomes with IMG annotations", [
                ["AP GOLD ID", "AP NAME", "AP TYPE", "AP IMG TAXON ID", "AP NCBI SRA IDS",
                 "AP STUDY GOLD ID", "AP PROJECT GOLD IDS", "BIOSAMPLE ECOSYSTEM PATH ID",
                 "BIOSAMPLE ECOSYSTEM", "BIOSAMPLE ECOSYSTEM CATEGORY", "BIOSAMPLE ECOSYSTEM TYPE",
                 "BIOSAMPLE ECOSYSTEM SUBTYPE", "BIOSAMPLE SPECIFIC ECOSYSTEM"],
                *sra_biome_rows,
            ]),
        ],
    )
    _write_xlsx(
        gold_dir / "GOLDs5levelEcosystemClassificationPaths.xlsx",
        [
            ("site data", [
                ["ECOSYSTEM PATH ID", "ECOSYSTEM", "ECOSYSTEM CATEGORY", "ECOSYSTEM TYPE",
                 "ECOSYSTEM SUBTYPE", "SPECIFIC ECOSYSTEM"],
                *eco_paths_rows,
            ]),
        ],
    )

    _write_tsv_gz(img_dir / "img_m_microbiomes_72col.tsv.gz", IMG_M_HEADER, img_rows)

    # SPIRE study TSV covering SAMN00000001 (marine) and a SPIRE-only sample SAMN00000099
    spire_rows = [
        ["SAMN00000001", "5", "35.5", "-140.2", "['aquatic', 'aquatic:marine']"],
        ["SAMN00000099", "2", "40.1", "-74.0", "['host-associated', 'host-associated:animal host']"],
    ]
    if reverse_order:
        spire_rows = list(reversed(spire_rows))
    _write_tsv(spire_dir / "study_alpha.tsv", ["sample_id", "mags", "lat", "lon", "microntology"], spire_rows)

    # Optional C-43 units TSV covering MGnify SRR100003
    c43_units = base / "c43_units.tsv"
    c43_rows = [
        ["mgnify", "ERZ9001", "take", "run", "1", "SRR100003", "SAMN00000003", "", "1000000", "0",
         "host-associated", "root:Host-associated:Human:Digestive system", "", "500000"],
    ]
    _write_tsv(
        c43_units,
        ["resource", "uid", "decision", "linked_by", "n_runs", "runs", "samples", "asm_bases",
         "read_bases", "runs_without_bases", "habitat", "habitat_detail", "note", "gz_bytes"],
        c43_rows,
    )

    return {
        "gold_data": gold_dir / "goldData.xlsx",
        "gold_sra_biome": gold_dir / "goldSraBiomeImg.xlsx",
        "gold_eco_paths": gold_dir / "GOLDs5levelEcosystemClassificationPaths.xlsx",
        "img_microbiomes": img_dir / "img_m_microbiomes_72col.tsv.gz",
        "spire_dir": spire_dir,
        "c43_units": c43_units,
    }


def test_streaming_xlsx_reader_inline_and_shared_strings(tmp_path):
    m = _mod()
    for use_ss in (False, True):
        p = tmp_path / f"sample_{'ss' if use_ss else 'inline'}.xlsx"
        _write_xlsx(
            p,
            [("Analysis Project", [
                ["AP GOLD ID", "AP IMG TAXON ID", "AP CONTACTS", "AP NCBI SRA IDS"],
                ["Ga0001", "3300000001", "alice@example.org", "SRR100001"],
                ["Ga0002", "", "bob@example.org", "SRR100002"],
            ])],
            use_shared_strings=use_ss,
        )
        rows = list(m.iter_xlsx_rows(p, "Analysis Project"))
        assert len(rows) == 2
        assert rows[0]["AP GOLD ID"] == "Ga0001"
        assert rows[0]["AP IMG TAXON ID"] == "3300000001"
        assert rows[1]["AP IMG TAXON ID"] == ""
        # PII/contact column must be excluded automatically
        assert "AP CONTACTS" not in rows[0]


@pytest.mark.parametrize("eco5, expected", [
    (("Host-associated", "Human", "Digestive system", "Large intestine", "Fecal"), "host-associated"),
    (("Host-associated", "Plants", "Rhizoplane", "Root", "Unclassified"), "host-associated"),
    (("Engineered", "Wastewater", "Activated sludge", "Municipal", "Unclassified"), "anthropogenic"),
    (("Engineered", "Bioreactor", "Anaerobic", "Unclassified", "Unclassified"), "anthropogenic"),
    (("Environmental", "Terrestrial", "Soil", "Forest soil", "Rhizosphere"), "terrestrial"),
    (("Environmental", "Aquatic", "Marine", "Pelagic", "Photic zone"), "aquatic"),
    (("Environmental", "Aquatic", "Freshwater", "Lake", "Sediment"), "aquatic"),
    (("Environmental", "Air", "Outdoor Air", "Unclassified", "Unclassified"), "air"),
    (("", "", "", "", ""), "unlabelled"),
])
def test_gold_5level_ecosystem_harmonization(eco5, expected):
    m = _mod()
    assert m.harmonize_gold_ecosystem(*eco5) == expected


def test_gold_img_sra_biosample_linking_and_flags(tmp_path):
    m = _mod()
    fx = _make_synthetic_bundle(tmp_path / "in")
    out_dir = tmp_path / "out"
    summary = m.build_crosswalk(
        gold_data=fx["gold_data"],
        gold_sra_biome=fx["gold_sra_biome"],
        gold_eco_paths=fx["gold_eco_paths"],
        img_microbiomes=fx["img_microbiomes"],
        spire_dir=fx["spire_dir"],
        c43_units=fx["c43_units"],
        out_dir=out_dir,
    )

    # Check crosswalk Parquet truth
    cw = pq.read_table(out_dir / "img_m_crosswalk.parquet").to_pylist()
    by_oid = {r["img_taxon_oid"]: r for r in cw}
    assert set(by_oid) == {"3300000001", "3300000002", "3300000003", "3300000004", "3300000006"}

    # 3300000001: marine, unrestricted, covered by SPIRE via SAMN00000001
    assert by_oid["3300000001"]["in_spire"] is True
    assert by_oid["3300000001"]["decision"] == "covered"
    assert by_oid["3300000001"]["is_unrestricted"] is True
    assert by_oid["3300000001"]["redistributable"] is False
    assert by_oid["3300000001"]["harmonized_habitat"] == "aquatic"

    # 3300000002: soil, unrestricted, linked to SRR100002 via GOLD AP (IMG/M had empty SRA Run!) -> novel take
    assert by_oid["3300000002"]["sra_runs"] == "SRR100002"
    assert by_oid["3300000002"]["in_spire"] is False
    assert by_oid["3300000002"]["in_mgnify"] is False
    assert by_oid["3300000002"]["decision"] == "take"
    assert by_oid["3300000002"]["asm_bases"] == 40_000_000
    assert by_oid["3300000002"]["harmonized_habitat"] == "terrestrial"

    # 3300000003: gut, unrestricted, covered by MGnify via SRR100003
    assert by_oid["3300000003"]["in_mgnify"] is True
    assert by_oid["3300000003"]["decision"] == "covered"

    # 3300000004: WWTP, Restricted -> decision == 'restricted'
    assert by_oid["3300000004"]["is_unrestricted"] is False
    assert by_oid["3300000004"]["decision"] == "restricted"
    assert by_oid["3300000004"]["harmonized_habitat"] == "anthropogenic"

    # 3300000006: Air, unrestricted, no SRA or BioSample -> take_unlinked
    assert by_oid["3300000006"]["decision"] == "take_unlinked"
    assert by_oid["3300000006"]["harmonized_habitat"] == "air"

    # Summary checks: novel unrestricted = 3300000002 (40M) + 3300000006 (12M) = 2 datasets, 52M bp
    assert summary["img_m_novel_unrestricted_units"] == 2
    assert summary["img_m_novel_unrestricted_bases"] == 52_000_000

    # Sidecar env_sample_dim redistributable flags: SPIRE True, IMG/M False
    con = duckdb.connect(str(out_dir / "env_sidecar.duckdb"), read_only=True)
    con.execute("SET threads=1")
    flags = dict(con.execute(
        "SELECT source, bool_and(redistributable) FROM env_sample_dim GROUP BY 1"
    ).fetchall())
    assert flags["img_m"] is False
    assert flags["spire"] is True
    # Check environmental attributes preserved
    ph, lat, lh = con.execute(
        "SELECT ph, latitude, longhurst_code FROM env_sample_dim WHERE source='img_m' AND source_id='3300000001'"
    ).fetchone()
    assert ph == pytest.approx(8.1)
    assert lat == pytest.approx(35.5)
    assert lh == "NPTE"
    con.close()


def test_deterministic_parquet_bytes_across_input_orders(tmp_path):
    m = _mod()
    fx_a = _make_synthetic_bundle(tmp_path / "in_a", reverse_order=False)
    fx_b = _make_synthetic_bundle(tmp_path / "in_b", reverse_order=True)
    out_a = tmp_path / "out_a"
    out_b = tmp_path / "out_b"

    for fx, out in ((fx_a, out_a), (fx_b, out_b)):
        m.build_crosswalk(
            gold_data=fx["gold_data"],
            gold_sra_biome=fx["gold_sra_biome"],
            gold_eco_paths=fx["gold_eco_paths"],
            img_microbiomes=fx["img_microbiomes"],
            spire_dir=fx["spire_dir"],
            c43_units=fx["c43_units"],
            out_dir=out,
        )

    for fname in ("env_sample_dim.parquet", "img_m_crosswalk.parquet"):
        ha = hashlib.sha256((out_a / fname).read_bytes()).hexdigest()
        hb = hashlib.sha256((out_b / fname).read_bytes()).hexdigest()
        assert ha == hb, f"{fname} bytes differed across input orders"


def test_bounded_rss_on_100x_synthetic_inputs(tmp_path):
    """100x synthetic inputs (100,000 rows across IMG/M + SPIRE + GOLD) stay under 512 MiB peak RSS."""
    in_dir = tmp_path / "scale_in"
    out_dir = tmp_path / "scale_out"
    gold_dir = in_dir / "gold"
    img_dir = in_dir / "img_m"
    spire_dir = in_dir / "spire"
    gold_dir.mkdir(parents=True)
    img_dir.mkdir(parents=True)
    spire_dir.mkdir(parents=True)

    n_img = 50_000
    n_spire = 50_000

    # Minimal GOLD xlsx with 5,000 AP rows (streaming XML)
    ap_rows = [
        ["AP GOLD ID", "AP NAME", "AP TYPE", "AP IMG TAXON ID", "AP GENBANK", "AP NCBI SRA IDS",
         "AP STUDY GOLD ID", "AP ORGANISM GOLD ID", "AP PROJECT GOLD IDS",
         "AP JGI DATA UTILIZATION STATUS\n (For projects sequenced at JGI only)",
         "AP GENE COUNT", "AP ESTIMATED SIZE"]
    ]
    for i in range(5_000):
        ap_rows.append([
            f"Ga{i:07d}", f"AP {i}", "Metagenome Analysis", str(3300000000 + i), "", f"SRR{i:07d}",
            f"Gs{i % 100:05d}", "", f"Gp{i:07d}", "Unrestricted", "1000", "1000000",
        ])
    _write_xlsx(gold_dir / "goldData.xlsx", [("Analysis Project", ap_rows)])

    # 50,000 IMG/M rows
    with gzip.open(img_dir / "img_m_microbiomes_72col.tsv.gz", "wt", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(IMG_M_HEADER) + "\n")
        for i in range(n_img):
            row = [""] * len(IMG_M_HEADER)
            row[0] = "*Microbiome"
            row[1] = str(3300000000 + i)
            row[2] = "Metagenome Analysis"
            row[3] = f"Ga{i:07d}"
            row[4] = f"Gp{i:07d}"
            row[5] = f"Gs{i % 500:05d}"
            row[9] = f"SAMN{i:08d}"
            row[11] = "Environmental"
            row[12] = "Aquatic" if i % 2 == 0 else "Terrestrial"
            row[13] = "Marine" if i % 2 == 0 else "Soil"
            row[14] = "Pelagic" if i % 2 == 0 else "Forest soil"
            row[15] = "Unclassified"
            row[34] = "1000000"
            row[36] = f"SRR{i:07d}"
            row[37] = "Unrestricted"
            fh.write("\t".join(row) + "\n")

    # 50,000 SPIRE rows (first 10,000 overlap IMG/M SAMN00000000..SAMN00009999)
    with open(spire_dir / "study_big.tsv", "w", encoding="utf-8", newline="") as fh:
        fh.write("sample_id\tmags\tlat\tlon\tmicrontology\n")
        for i in range(n_spire):
            sid = f"SAMN{i:08d}" if i < 10_000 else f"SAMEA{i:08d}"
            fh.write(f"{sid}\t1\t0.0\t0.0\t['aquatic', 'aquatic:marine']\n")

    # Run in isolated subprocess so peak RSS measures only the builder. The child's own peak is VmHWM of
    # /proc/self/status (the high-water mark since its exec), as tests/r2/test_w39_panel_fixture.py reads it since card
    # R-01: on Linux a child's `ru_maxrss` carries the peak of the pytest process it was started from, so it failed this
    # bound on GitHub Actions once an earlier test had raised the suite's peak (575.82 MiB, run 45754246039; 2,893.64
    # MiB on wilkcomp02 under a parent that had peaked at 3 GB, while the child's VmHWM was 440.45 MiB; measured, card
    # R-10, 2026-11-25). `ru_maxrss` is used only where there is no /proc (macOS), and is printed beside VmHWM.
    probe = (
        "import json, os, resource, sys, importlib.util; "
        f"spec = importlib.util.spec_from_file_location('c48b', {str(SCRIPT)!r}); "
        "m = importlib.util.module_from_spec(spec); sys.modules['c48b'] = m; spec.loader.exec_module(m); "
        f"res = m.build_crosswalk(gold_data={str(gold_dir / 'goldData.xlsx')!r}, "
        f"img_microbiomes={str(img_dir / 'img_m_microbiomes_72col.tsv.gz')!r}, "
        f"spire_dir={str(spire_dir)!r}, out_dir={str(out_dir)!r}); "
        "rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss; "
        "maxrss_mb = (rss / (1024 * 1024)) if sys.platform == 'darwin' else (rss / 1024.0); "
        "st = '/proc/self/status'; "
        "hwm_mb = (next(int(x.split()[1]) / 1024.0 for x in open(st) if x.startswith('VmHWM:')) "
        "if os.path.exists(st) else None); "
        "print(json.dumps({'hwm_mb': hwm_mb if hwm_mb is None else round(hwm_mb, 2), "
        "'maxrss_mb': round(maxrss_mb, 2), 'summary': res}))"
    )
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    proc = subprocess.run([sys.executable, "-I", "-c", probe], capture_output=True, text=True, env=env, check=True)
    payload = json.loads(proc.stdout.strip().splitlines()[-1])
    peak_mb = payload["hwm_mb"] if payload["hwm_mb"] is not None else payload["maxrss_mb"]
    print(f"\nc48b 100x crosswalk in a child: VmHWM {payload['hwm_mb']} MiB, "
          f"ru_maxrss {payload['maxrss_mb']} MiB (with its parent's peak)")
    assert peak_mb < 512.0, (f"Peak RSS too high: {peak_mb} MiB (VmHWM {payload['hwm_mb']}, "
                             f"ru_maxrss {payload['maxrss_mb']})")
    assert payload["summary"]["img_m_total_units"] == n_img
    assert payload["summary"]["img_m_novel_unrestricted_units"] == n_img - 10_000
