#!/usr/bin/env python3
"""
merge_og_stats_to_final_gene_table.py
=====================================
Aggregate the protein/OG-level enrichment statistics to the gene level and
merge them into the symbol-prediction table. Species-agnostic (driven by
`species:` in the config).

Per gene, a representative OG (Primary_OG) is chosen by: smallest FDR, then
highest immune_enrichment, then highest ref_total_genes.

Inputs : Final_Filtered_List_with_OG_stats.csv,
         Immune_Predict_Result_Final.csv,
         GFF3 (transcript -> gene mapping)
Outputs: Immune_Predict_Result_Final_with_OG_stats.csv,
         Unmapped_OG_ProteinIDs.txt
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species


def parse_gff_transcript_to_gene(gff_path):
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


def clean_pp_id(x):
    if pd.isna(x):
        return ""
    return str(x).strip().replace("transcript:", "").replace("mRNA:", "").replace("transcript_", "")


def unique_join(values, sep="; "):
    seen, vals = set(), []
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if s == "" or s.lower() == "nan" or s in seen:
            continue
        seen.add(s)
        vals.append(s)
    return sep.join(vals)


def choose_primary_row(group):
    g = group.copy()
    for col in ("FDR", "immune_enrichment", "ref_total_genes"):
        g[col] = pd.to_numeric(g[col], errors="coerce")
    g = g.sort_values(by=["FDR", "immune_enrichment", "ref_total_genes"],
                      ascending=[True, False, False], na_position="last")
    return g.iloc[0]


def make_aggregate_gene(reference_species):
    """Build aggregate_gene(group) with per-species immune columns from config."""
    def aggregate_gene(group):
        primary = choose_primary_row(group)
        out = {
            # Newer pandas excludes the grouping column from each group, so read
            # the GeneID from the group key instead of a row of the group.
            "GeneID": group.name,
            "OG_List": unique_join(group["Orthogroup"].dropna().astype(str).unique()),
            "OG_Count": int(group["Orthogroup"].dropna().astype(str).nunique()),
            "Primary_OG": primary.get("Orthogroup", np.nan),
            "OG_Filter_Reason_List": unique_join(
                group["Filter_Reason"] if "Filter_Reason" in group.columns else []),
            "Primary_immune_enrichment": primary.get("immune_enrichment", np.nan),
            "Primary_FDR": primary.get("FDR", np.nan),
            "Primary_immune_fraction": primary.get("immune_fraction", np.nan),
            "Primary_ref_total_genes": primary.get("ref_total_genes", np.nan),
            "Primary_dup_count": primary.get("dup_count", np.nan),
            "Primary_dup_density": primary.get("dup_density", np.nan),
            "Primary_copy_variance": primary.get("copy_variance", np.nan),
            "Primary_immune_entropy": primary.get("immune_entropy", np.nan),
            "max_immune_enrichment": pd.to_numeric(group["immune_enrichment"], errors="coerce").max(),
            "min_FDR": pd.to_numeric(group["FDR"], errors="coerce").min(),
            "max_dup_density": pd.to_numeric(group["dup_density"], errors="coerce").max(),
            "max_copy_variance": pd.to_numeric(group["copy_variance"], errors="coerce").max(),
        }
        for species in reference_species:
            out[f"Primary_immune_{species}"] = primary.get(f"immune_{species}", np.nan)
        return pd.Series(out)
    return aggregate_gene


def insert_columns_after(df, after_col, new_cols):
    existing = list(df.columns)
    if after_col not in existing:
        return df
    present_new = [c for c in new_cols if c in df.columns
                   and c not in existing[:existing.index(after_col) + 1]]
    base_cols = [c for c in existing if c not in present_new]
    insert_idx = base_cols.index(after_col) + 1
    return df[base_cols[:insert_idx] + present_new + base_cols[insert_idx:]]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--og-stats", help="..._with_OG_stats.csv")
    p.add_argument("--final", help="..._Predict_Result_Final.csv")
    p.add_argument("--gff", help="GFF3 annotation")
    p.add_argument("--output", help="Output CSV with OG stats")
    p.add_argument("--unmapped", help="Unmapped protein IDs TXT")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    reference_species = sp["ref_names"]
    target_id_col = f"{sp['target']}_ProteinID"

    orth = cfg.get("orthofinder", {})
    tier = cfg.get("tiering", {})
    orth_work = orth.get("work_dir", "")
    tier_work = tier.get("work_dir", "")

    def owd(name):
        return os.path.join(orth_work, name) if orth_work else name

    def twd(name):
        return os.path.join(tier_work, name) if tier_work else name

    og_stats_file = args.og_stats or owd(orth.get("filtered_with_og_stats",
                                                 "Final_Filtered_List_with_OG_stats.csv"))
    final_file = args.final or twd(tier.get("predict_result",
                                           "Immune_Predict_Result_Final.csv"))
    gff_file = args.gff or cfg.get("reference", {}).get("gff3_raw", "annotation.gff3")
    output_file = args.output or twd(tier.get("predict_with_og_stats",
                                             "Immune_Predict_Result_Final_with_OG_stats.csv"))
    unmapped_file = args.unmapped or twd(tier.get("unmapped_og_proteins", "Unmapped_OG_ProteinIDs.txt"))

    print("Loading files...")
    df_og = pd.read_csv(og_stats_file, encoding="utf-8-sig")
    df_final = pd.read_csv(final_file, encoding="utf-8-sig")
    print(f"    OG rows: {len(df_og)}; final gene rows: {len(df_final)}")

    tx2gene = parse_gff_transcript_to_gene(gff_file)
    df_og["Clean_ID"] = df_og[target_id_col].apply(clean_pp_id)
    df_og["GeneID"] = df_og["Clean_ID"].map(tx2gene)

    unmapped = df_og[df_og["GeneID"].isna()].copy()
    mapped = df_og.dropna(subset=["GeneID"]).copy()
    print(f"    mapped OG rows: {len(mapped)} / {len(df_og)} (unmapped: {len(unmapped)})")

    if len(unmapped):
        with open(unmapped_file, "w") as fh:
            fh.write(f"{target_id_col}\tClean_ID\tOrthogroup\n")
            for _, r in unmapped[[target_id_col, "Clean_ID", "Orthogroup"]].drop_duplicates().iterrows():
                fh.write(f"{r[target_id_col]}\t{r['Clean_ID']}\t{r['Orthogroup']}\n")
        print(f"    wrote unmapped list: {unmapped_file}")

    print("Aggregating OG stats to gene level...")
    aggregate_gene = make_aggregate_gene(reference_species)
    gene_og_stats = mapped.groupby("GeneID", group_keys=False).apply(aggregate_gene).reset_index(drop=True)
    print(f"    gene-level rows: {len(gene_og_stats)}")

    merged = df_final.merge(gene_og_stats, on="GeneID", how="left")

    og_cols = [
        "OG_List", "OG_Count", "Primary_OG", "OG_Filter_Reason_List",
        "Primary_immune_enrichment", "Primary_FDR", "Primary_immune_fraction",
        "Primary_ref_total_genes", "Primary_dup_count", "Primary_dup_density",
        "Primary_copy_variance", "Primary_immune_entropy"]
    og_cols += [f"Primary_immune_{s}" for s in reference_species]
    og_cols += ["max_immune_enrichment", "min_FDR", "max_dup_density", "max_copy_variance"]

    if "Evidence_Count" in merged.columns:
        merged = insert_columns_after(merged, "Evidence_Count", og_cols)
    elif "Tier" in merged.columns:
        merged = insert_columns_after(merged, "Tier", og_cols)

    merged.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved to {output_file}")
    print(f"Genes with OG stats: {merged['Primary_OG'].notna().sum()} / {len(merged)}")


if __name__ == "__main__":
    main()
