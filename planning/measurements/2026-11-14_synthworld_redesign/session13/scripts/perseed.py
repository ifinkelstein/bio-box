"""perseed.py SCAN_TSV: the per-seed Markdown table of one scan (card F-02 session 13, scratch; session 12's §15.4
columns plus the decoy columns of DEC-48 (1)). Formatting only: every number is read from the scan's TSV."""
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
cols = ["seed", "crown age (MRCA) / from first speciation", "genes per genome mean (min-max)",
        "chromosome mean (min-max)",
        "planted / profile-only / decoy / vehicle / hub / plasmid / background share %",
        "plasmid genomes pls0, pls1 (0/1/2 per genome)", "modules all members in >= 9 genomes",
        "copies per present member", "realized transfers per member, level medians", "N0 background size ratio, r",
        "N0 root KS D full, strict (bound)", "module transfers into a holder",
        "decoy transfers into a partner carrier, per decoy", "partner share", "decoy genes per genome",
        "decoy copies per present genome (max)"]
print("| " + " | ".join(cols) + " |")
print("|" + "---|" * len(cols))
for r in rows:
    sh = " / ".join(f"{100 * float(r['share_' + k]):.1f}" for k in
                    ("planted", "profile_only", "decoy", "vehicle", "hub", "plasmid", "background"))
    print("| " + " | ".join([
        r["seed"], f"{r['crown_age_mrca']} / {r['crown_age']}", f"{r['genes_mean']} ({r['genes_min']}-{r['genes_max']})",
        f"{r['chrom_mean']} ({r['chrom_min']}-{r['chrom_max']})", sh,
        f"{r['pls0_genomes']}, {r['pls1_genomes']} ({r['genomes_0_1_2_plasmids']})", r["modules_all_members_ge9"],
        r["module_copies_mean"], f"{r['tpm_level0_median']}, {r['tpm_level1_median']}, {r['tpm_level2_median']}",
        f"{r['n0_bg_ratio']}, {r['n0_bg_pearson']}",
        f"{r['n0full_root_ks_d']}, {r['n0strict_root_ks_d']} ({r['n0full_root_ks_crit']})",
        f"0 of {r['module_transfer_rows']}" if r["module_transfers_into_first_holder"] == "0"
        else f"{r['module_transfers_into_first_holder']} of {r['module_transfer_rows']}",
        r["decoy_transfers_per_decoy"], r["decoy_partner_share_mean"], r["per_genome_decoy"],
        f"{r['decoy_copies_per_present_genome']} ({r['decoy_copies_max']})"]) + " |")
