#!/usr/bin/env python3
"""
merge_immune_annotations.py
===========================
Merge the three evidence lines (InterProScan, OrthoFinder, KofamScan) at the
gene level (outer join) and assign a confidence Tier.

Tier 1: all three evidence sources    Tier 2: any two    Tier 3: a single source

Inputs : interproscan_immune_results.csv, PinkPigeon_Final_Filtered_List.csv,
         final_kofam_annotated.tsv, GFF3 (for transcript->gene mapping)
Output : PinkPigeon_Immune_Gene_Master_List.csv
"""

import argparse
import csv
import os

import pandas as pd


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def parse_gff_transcript_to_gene(gff_path):
    """Tolerant GFF parser: transcript ID -> gene ID (prefixes stripped)."""
    print(f"Parsing GFF: {gff_path}")
    tx_to_gene = {}
    feature_counts = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_counts[parts[2]] = feature_counts.get(parts[2], 0) + 1
            attributes = parts[8]
            if "ID=" in attributes and "Parent=" in attributes:
                attr = {k: v for k, v in (item.split("=", 1)
                        for item in attributes.split(";") if "=" in item)}
                if "ID" in attr and "Parent" in attr:
                    tx = attr["ID"].replace("transcript:", "").replace("mRNA:", "").replace("transcript_", "")
                    gene = attr["Parent"].replace("gene:", "").replace("gene_", "")
                    tx_to_gene[tx] = gene
    print(f"    transcript->gene mappings: {len(tx_to_gene)}")
    return tx_to_gene


def clean_and_aggregate(df, group_col, agg_cols):
    """Group by gene and join unique non-empty values of each column with '; '."""
    rules = {}
    for col in agg_cols:
        if col in df.columns:
            rules[col] = lambda x: "; ".join(sorted(set(
                str(v).strip() for v in x if pd.notna(v) and str(v).strip() != "")))
    return df.groupby(group_col).agg(rules).reset_index()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--gff", help="GFF3 annotation")
    p.add_argument("--interpro", help="interproscan_immune_results.csv")
    p.add_argument("--ortho", help="PinkPigeon_Final_Filtered_List.csv")
    p.add_argument("--kofam", help="final_kofam_annotated.tsv")
    p.add_argument("--output", help="Master list CSV")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    def stage_path(section, key, default):
        node = cfg.get(section, {})
        work = node.get("work_dir", "")
        name = node.get(key, default)
        return os.path.join(work, name) if work else name

    gff_file = args.gff or cfg.get("reference", {}).get("gff3_raw", "annotation.gff3")
    interpro_file = args.interpro or stage_path("interproscan", "immune_results",
                                               "interproscan_immune_results.csv")
    ortho_file = args.ortho or stage_path("orthofinder", "filtered_list",
                                         "PinkPigeon_Final_Filtered_List.csv")
    kofam_file = args.kofam or stage_path("kofamscan", "annotated", "final_kofam_annotated.tsv")
    output_file = args.output or stage_path("tiering", "master_list",
                                           "PinkPigeon_Immune_Gene_Master_List.csv")

    tx2gene = parse_gff_transcript_to_gene(gff_file)

    # --- InterProScan ---
    print("Processing InterProScan...")
    ipr_df = pd.read_csv(interpro_file)
    ipr_df["Clean_ID"] = ipr_df["ID"].astype(str).str.replace("transcript:", "", regex=False).str.replace("transcript_", "", regex=False)
    ipr_df["GeneID"] = ipr_df["Clean_ID"].map(tx2gene)
    unmapped = ipr_df["GeneID"].isna().sum()
    if unmapped:
        print(f"    [WARN] {unmapped} InterProScan entries failed to map; dropping them.")
        ipr_df = ipr_df.dropna(subset=["GeneID"])
    ipr_gene = clean_and_aggregate(ipr_df, "GeneID",
                                   ["signature_description", "interpro_annotations_description"])
    ipr_gene["Source_Interproscan"] = True

    # --- OrthoFinder ---
    print("Processing OrthoFinder...")
    ortho_df = pd.read_csv(ortho_file)
    ortho_df["Clean_ID"] = ortho_df["PinkPigeon_ProteinID"].astype(str).str.replace("transcript:", "", regex=False).str.replace("transcript_", "", regex=False)
    ortho_df["GeneID"] = ortho_df["Clean_ID"].map(tx2gene)
    unmapped = ortho_df["GeneID"].isna().sum()
    if unmapped:
        print(f"    [WARN] {unmapped} OrthoFinder entries failed to map; dropping them.")
        ortho_df = ortho_df.dropna(subset=["GeneID"])
    ortho_cols = [
        "Total_Score", "Filter_Reason",
        "Mouse_GeneSymbols", "Mouse_Category1", "Mouse_Subcategory",
        "Mouse_Description", "Mouse_UniProt_Function", "Mouse_Immune_Source",
        "Chicken_GeneSymbols", "Chicken_Description", "Chicken_GO_Terms",
        "ZebraFinch_GeneSymbols", "ZebraFinch_Description", "ZebraFinch_GO_Terms"]
    ortho_gene = clean_and_aggregate(ortho_df, "GeneID", ortho_cols)
    ortho_gene["Source_Orthofinder"] = True

    # --- KofamScan (already gene-level) ---
    print("Processing KofamScan...")
    kofam_df = pd.read_csv(kofam_file, sep="\t")
    kofam_df.columns = [c.strip() for c in kofam_df.columns]
    kofam_cols = ["KO_Gene_Symbol", "KO_Gene_Name", "KO", "Immune_Pathway", "KO_definition"]
    kofam_gene = clean_and_aggregate(kofam_df, "GeneID", kofam_cols)
    kofam_gene["Source_Kofam"] = True

    # --- Merge ---
    print("Merging...")
    merged = pd.merge(ipr_gene, ortho_gene, on="GeneID", how="outer")
    merged = pd.merge(merged, kofam_gene, on="GeneID", how="outer")

    for col in ("Source_Interproscan", "Source_Orthofinder", "Source_Kofam"):
        merged[col] = merged[col].fillna(False)

    def sources(row):
        s = []
        if row["Source_Interproscan"]:
            s.append("Interproscan")
        if row["Source_Orthofinder"]:
            s.append("Orthofinder")
        if row["Source_Kofam"]:
            s.append("Kofam")
        return ";".join(s)

    merged["Evidence_Sources"] = merged.apply(sources, axis=1)
    merged["Evidence_Count"] = (merged["Source_Interproscan"].astype(int)
                                + merged["Source_Orthofinder"].astype(int)
                                + merged["Source_Kofam"].astype(int))

    tier_map = {3: "Tier 1", 2: "Tier 2", 1: "Tier 3"}
    merged["Tier"] = merged["Evidence_Count"].map(lambda c: tier_map.get(c, "Unknown"))

    lead = ["GeneID", "Tier", "Evidence_Sources", "Evidence_Count"]
    drop = set(lead) | {"Source_Interproscan", "Source_Orthofinder", "Source_Kofam"}
    cols = lead + [c for c in merged.columns if c not in drop]
    merged = merged[cols]

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    merged.to_csv(output_file, index=False, encoding="utf-8-sig",
                  quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved {len(merged)} genes to {output_file}")
    print("Tier distribution:")
    print(merged["Tier"].value_counts())


if __name__ == "__main__":
    main()
