#!/usr/bin/env python3
"""
merge_orthology_to_final_gene_table.py
======================================
Aggregate the transcript-level orthology analysis to the gene level, merge it
into the final immune-gene prediction table, and add the Orth_Sub_Tier column.

Per gene (multiple transcripts) the best value is kept per species:
  - Orth_{Species}_type: highest priority type (1:1 > 1:N > N:1 > N:M)
  - counts: maximum across transcripts
Orth_Sub_Tier = number of distinct species providing any high-weight evidence
(1:1 immune, or 1:N all-immune); union, deduplicated, range 0-3.

Inputs : PinkPigeon_Final_Filtered_List_with_OG_stats_and_orthology.csv,
         PinkPigeon_Immune_Predict_Result_Final_with_OG_stats.csv,
         GFF3 (transcript -> gene mapping)
Output : PinkPigeon_Immune_Predict_Result_Final_with_OG_and_Orthology.csv
"""

import argparse
import csv
import os

import numpy as np
import pandas as pd

SPECIES_LIST = ["Mouse", "Chicken", "ZebraFinch"]
ORTH_TYPE_PRIORITY = {"1:1": 0, "1:N": 1, "N:1": 2, "N:M": 3}


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def parse_gff_transcript_to_gene(gff_path):
    print(f"Parsing GFF: {gff_path}")
    tx_to_gene = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
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


def orth_type_rank(t):
    if pd.isna(t):
        return 99
    return ORTH_TYPE_PRIORITY.get(str(t), 99)


def aggregate_orthology_to_gene(group):
    """Collapse a gene's transcripts into one gene-level orthology summary."""
    out = {"Gene_Transcript_Count": len(group)}

    if len(group) == 1:
        row = group.iloc[0]
        for sp in SPECIES_LIST:
            prefix = f"Orth_{sp}"
            for suffix in ("type", "ref_total", "ref_immune", "all_immune", "rows"):
                out[f"{prefix}_{suffix}"] = row.get(f"{prefix}_{suffix}", np.nan)
        out["Orth_1to1_immune_count"] = row.get("Orth_1to1_immune_count", 0)
        out["Orth_1to1_immune_species"] = row.get("Orth_1to1_immune_species", np.nan)
        out["Orth_1toN_allimmune_count"] = row.get("Orth_1toN_allimmune_count", 0)
        out["Orth_1toN_allimmune_species"] = row.get("Orth_1toN_allimmune_species", np.nan)
        return pd.Series(out)

    for sp in SPECIES_LIST:
        prefix = f"Orth_{sp}"
        type_col = f"{prefix}_type"
        rows_col = f"{prefix}_rows"
        valid = group[group[type_col].notna()].copy()
        if len(valid) == 0:
            out[f"{prefix}_type"] = np.nan
            out[f"{prefix}_ref_total"] = np.nan
            out[f"{prefix}_ref_immune"] = np.nan
            out[f"{prefix}_all_immune"] = np.nan
            out[f"{prefix}_rows"] = 0
        else:
            valid["_type_rank"] = valid[type_col].apply(orth_type_rank)
            valid[f"{prefix}_ref_immune"] = pd.to_numeric(valid[f"{prefix}_ref_immune"], errors="coerce")
            valid = valid.sort_values(by=["_type_rank", f"{prefix}_ref_immune"],
                                      ascending=[True, False])
            best = valid.iloc[0]
            out[f"{prefix}_type"] = best[type_col]
            out[f"{prefix}_ref_total"] = best.get(f"{prefix}_ref_total", np.nan)
            out[f"{prefix}_ref_immune"] = best.get(f"{prefix}_ref_immune", np.nan)
            out[f"{prefix}_all_immune"] = best.get(f"{prefix}_all_immune", np.nan)
            out[f"{prefix}_rows"] = pd.to_numeric(group[rows_col], errors="coerce").max()

    count_1to1 = pd.to_numeric(group["Orth_1to1_immune_count"], errors="coerce")
    count_1toN = pd.to_numeric(group["Orth_1toN_allimmune_count"], errors="coerce")
    max_1to1_idx = count_1to1.idxmax() if count_1to1.notna().any() else None
    max_1toN_idx = count_1toN.idxmax() if count_1toN.notna().any() else None

    out["Orth_1to1_immune_count"] = int(count_1to1.max()) if count_1to1.notna().any() else 0
    out["Orth_1to1_immune_species"] = (group.loc[max_1to1_idx, "Orth_1to1_immune_species"]
                                       if max_1to1_idx is not None else np.nan)
    out["Orth_1toN_allimmune_count"] = int(count_1toN.max()) if count_1toN.notna().any() else 0
    out["Orth_1toN_allimmune_species"] = (group.loc[max_1toN_idx, "Orth_1toN_allimmune_species"]
                                         if max_1toN_idx is not None else np.nan)
    return pd.Series(out)


def compute_sub_tier(row):
    """Distinct species providing any high-weight orthology evidence (0-3)."""
    species_set = set()
    for col in ("Orth_1to1_immune_species", "Orth_1toN_allimmune_species"):
        val = row.get(col, None)
        if pd.notna(val) and str(val).strip() != "":
            for sp in str(val).split(";"):
                sp = sp.strip()
                if sp:
                    species_set.add(sp)
    return len(species_set)


def insert_columns_after(df, after_col, new_cols):
    existing = list(df.columns)
    if after_col not in existing:
        return df
    present_new = [c for c in new_cols if c in df.columns]
    base_cols = [c for c in existing if c not in present_new]
    insert_idx = base_cols.index(after_col) + 1
    return df[base_cols[:insert_idx] + present_new + base_cols[insert_idx:]]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--orthology", help="..._with_OG_stats_and_orthology.csv (transcript level)")
    p.add_argument("--final", help="..._Predict_Result_Final_with_OG_stats.csv (gene level)")
    p.add_argument("--gff", help="GFF3 annotation")
    p.add_argument("--output", help="Final table CSV")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    orth = cfg.get("orthofinder", {})
    tier = cfg.get("tiering", {})
    orth_work = orth.get("work_dir", "")
    tier_work = tier.get("work_dir", "")

    def owd(name):
        return os.path.join(orth_work, name) if orth_work else name

    def twd(name):
        return os.path.join(tier_work, name) if tier_work else name

    orthology_file = args.orthology or owd(orth.get("filtered_with_orthology",
                                                   "PinkPigeon_Final_Filtered_List_with_OG_stats_and_orthology.csv"))
    final_file = args.final or twd(tier.get("predict_with_og_stats",
                                           "PinkPigeon_Immune_Predict_Result_Final_with_OG_stats.csv"))
    gff_file = args.gff or cfg.get("reference", {}).get("gff3_raw", "annotation.gff3")
    output_file = args.output or twd(tier.get("final_table",
                                             "PinkPigeon_Immune_Predict_Result_Final_with_OG_and_Orthology.csv"))

    print("Loading data...")
    df_orth = pd.read_csv(orthology_file, encoding="utf-8-sig")
    df_final = pd.read_csv(final_file, encoding="utf-8-sig")
    print(f"    orthology transcript rows: {len(df_orth)}; final gene rows: {len(df_final)}")

    tx2gene = parse_gff_transcript_to_gene(gff_file)
    df_orth["Clean_ID"] = df_orth["PinkPigeon_ProteinID"].apply(clean_pp_id)
    df_orth["GeneID"] = df_orth["Clean_ID"].map(tx2gene)
    mapped = df_orth.dropna(subset=["GeneID"]).copy()
    print(f"    mapped: {len(mapped)} / {len(df_orth)}")

    print("Aggregating orthology to gene level...")
    gene_orth = (mapped.groupby("GeneID", group_keys=False)
                 .apply(aggregate_orthology_to_gene).reset_index())
    if "GeneID" not in gene_orth.columns and gene_orth.index.name == "GeneID":
        gene_orth = gene_orth.reset_index()

    for col in ("Orth_1to1_immune_count", "Orth_1toN_allimmune_count"):
        gene_orth[col] = pd.to_numeric(gene_orth[col], errors="coerce").fillna(0).astype(int)
    gene_orth["Orth_Sub_Tier"] = gene_orth.apply(compute_sub_tier, axis=1).astype(int)

    print("Merging into final gene table...")
    orth_merge_cols = ["GeneID", "Gene_Transcript_Count"]
    for sp in SPECIES_LIST:
        prefix = f"Orth_{sp}"
        orth_merge_cols += [f"{prefix}_type", f"{prefix}_ref_total", f"{prefix}_ref_immune",
                            f"{prefix}_all_immune", f"{prefix}_rows"]
    orth_merge_cols += ["Orth_1to1_immune_count", "Orth_1to1_immune_species",
                        "Orth_1toN_allimmune_count", "Orth_1toN_allimmune_species", "Orth_Sub_Tier"]

    merged = df_final.merge(gene_orth[orth_merge_cols], on="GeneID", how="left")

    # Genes with no OrthoFinder evidence get Orth_Sub_Tier = NaN
    if "Evidence_Sources" in merged.columns:
        no_ortho = ~merged["Evidence_Sources"].str.contains("Orthofinder", na=True)
        merged.loc[no_ortho, "Orth_Sub_Tier"] = np.nan

    insert_cols = [c for c in orth_merge_cols if c != "GeneID"]
    if "max_copy_variance" in merged.columns:
        merged = insert_columns_after(merged, "max_copy_variance", insert_cols)
    elif "Primary_OG" in merged.columns:
        merged = insert_columns_after(merged, "Primary_OG", insert_cols)

    merged.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved to {output_file}  ({len(merged)} rows)")
    if "Orth_Sub_Tier" in merged.columns:
        print("Orth_Sub_Tier distribution:")
        print(merged["Orth_Sub_Tier"].value_counts(dropna=False).sort_index().to_string())


if __name__ == "__main__":
    main()
