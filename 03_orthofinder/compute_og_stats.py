#!/usr/bin/env python3
"""
compute_og_stats.py
===================
Per-orthogroup immune-enrichment and structural statistics for the candidate
orthogroups. Species-agnostic: reference species and the target species come
from `species:` in the config.

For each candidate OG it computes, from the reference species only:
  - immune fraction and enrichment vs a genome-wide background
  - one-sided Fisher exact test (hypergeometric, no scipy) + Benjamini-Hochberg FDR
  - normalized Shannon entropy of immune counts across the reference species
  - gene-copy variance / SD and duplication density

Inputs : Orthogroups.tsv, Orthogroups.GeneCount.tsv, Duplications.tsv,
         master_lookup_table.csv, PinkPigeon_Final_Filtered_List.csv
Outputs: OG_stats_summary.tsv,
         PinkPigeon_Final_Filtered_List_with_OG_stats.csv
"""

import argparse
import csv
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species


def clean_protein_id(pid):
    if pd.isna(pid):
        return ""
    pid = str(pid).strip()
    return pid.split(".")[0] if "." in pid else pid


def split_ids(cell):
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    return [x.strip() for x in s.split(",") if x.strip()] if s else []


def safe_div(a, b):
    return a / b if b not in (0, 0.0) else np.nan


def compute_shannon_entropy(counts, n_categories):
    """Normalized Shannon entropy H/log(n_categories); NaN when total is 0."""
    total = sum(counts)
    if total == 0 or n_categories < 2:
        return np.nan
    probs = [c / total for c in counts if c > 0]
    H = -sum(p * math.log(p) for p in probs)
    return H / math.log(n_categories)


def fisher_exact_right_tail(a, b, c, d):
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    n = row1 + row2
    if b == 0 or c == 0:
        odds_ratio = np.nan if (a == 0 or d == 0) else np.inf
    else:
        odds_ratio = (a * d) / (b * c)
    x_max = min(row1, col1)

    def log_comb(n_, k_):
        if k_ < 0 or k_ > n_:
            return -np.inf
        return math.lgamma(n_ + 1) - math.lgamma(k_ + 1) - math.lgamma(n_ - k_ + 1)

    def hypergeom_p(x):
        return math.exp(log_comb(col1, x) + log_comb(col2, row1 - x) - log_comb(n, row1))

    p_value = sum(hypergeom_p(x) for x in range(a, x_max + 1))
    return min(p_value, 1.0), odds_ratio


def benjamini_hochberg(pvalues):
    pvals = np.array(pvalues, dtype=float)
    qvals = np.full(len(pvals), np.nan)
    valid = ~np.isnan(pvals)
    pv = pvals[valid]
    m = len(pv)
    if m == 0:
        return qvals
    order = np.argsort(pv)
    ranked = pv[order]
    q = np.empty(m, dtype=float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ranked[i] * m / (i + 1))
        q[i] = min(prev, 1.0)
    q_original = np.empty(m, dtype=float)
    q_original[order] = q
    qvals[valid] = q_original
    return qvals


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--orthogroups", help="Orthogroups.tsv")
    p.add_argument("--genecount", help="Orthogroups.GeneCount.tsv")
    p.add_argument("--duplications", help="Duplications.tsv")
    p.add_argument("--master", help="master_lookup_table.csv")
    p.add_argument("--filtered", help="PinkPigeon_Final_Filtered_List.csv")
    p.add_argument("--out-summary", help="OG_stats_summary.tsv")
    p.add_argument("--out-merged", help="..._with_OG_stats.csv")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    reference_species = sp["ref_names"]
    target = sp["target"]
    all_species = [target] + reference_species
    n_ref = len(reference_species)

    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    orthogroups_tsv = args.orthogroups or wd(orth.get("orthogroups_tsv", "Orthogroups.tsv"))
    genecount_tsv = args.genecount or wd(orth.get("genecount_tsv", "Orthogroups.GeneCount.tsv"))
    dup_tsv = args.duplications or wd(orth.get("duplications_tsv", "Duplications.tsv"))
    master_csv = args.master or wd(orth.get("master_lookup", "master_lookup_table.csv"))
    filtered_csv = args.filtered or wd(orth.get("filtered_list", "PinkPigeon_Final_Filtered_List.csv"))
    out_summary = args.out_summary or wd(orth.get("og_stats_summary", "OG_stats_summary.tsv"))
    out_merged = args.out_merged or wd(orth.get("filtered_with_og_stats",
                                               "PinkPigeon_Final_Filtered_List_with_OG_stats.csv"))

    print(">>> Loading files...")
    df_filtered = pd.read_csv(filtered_csv, encoding="utf-8-sig")
    df_og = pd.read_csv(orthogroups_tsv, sep="\t", dtype=str)
    df_gc = pd.read_csv(genecount_tsv, sep="\t")
    df_dup = pd.read_csv(dup_tsv, sep="\t", dtype=str)
    try:
        df_master = pd.read_csv(master_csv, encoding="utf-8-sig")
    except UnicodeDecodeError:
        df_master = pd.read_csv(master_csv, encoding="utf-8")

    candidate_ogs = set(df_filtered["Orthogroup"].dropna().astype(str).unique())
    print(f"    Candidate OGs: {len(candidate_ogs)}")

    df_master["ProteinID_Clean"] = df_master["ProteinID"].apply(clean_protein_id)
    immune_protein_set = set(df_master["ProteinID_Clean"].dropna().astype(str))
    print(f"    Immune protein IDs: {len(immune_protein_set)}")

    # --- Global background across all OGs (reference species only) ---
    print(">>> Computing global background...")
    total_ref_genes_background = 0
    total_immune_genes_background = 0
    for _, row in df_og.iterrows():
        for species in reference_species:
            ids_clean = [clean_protein_id(x) for x in split_ids(row.get(species, ""))]
            total_ref_genes_background += len(ids_clean)
            total_immune_genes_background += sum(1 for x in ids_clean if x in immune_protein_set)
    background_fraction = safe_div(total_immune_genes_background, total_ref_genes_background)
    print(f"    Background immune fraction: {background_fraction:.6f}")

    dup_count_df = df_dup.groupby("Orthogroup").size().reset_index(name="dup_count")

    # --- Gene-copy counts (candidate OGs only) ---
    needed = ["Orthogroup"] + all_species + ["Total"]
    missing = [c for c in needed if c not in df_gc.columns]
    if missing:
        raise ValueError(f"Missing columns in GeneCount.tsv: {missing} "
                         f"(expected the config species names)")
    df_gc_sub = df_gc[df_gc["Orthogroup"].astype(str).isin(candidate_ogs)].copy()
    for col in all_species + ["Total"]:
        df_gc_sub[col] = pd.to_numeric(df_gc_sub[col], errors="coerce").fillna(0).astype(int)
    df_gc_sub["total_genes_all"] = df_gc_sub[all_species].sum(axis=1)

    def copy_array(row):
        return np.array([row[s] for s in all_species], dtype=float)

    df_gc_sub["copy_variance"] = df_gc_sub.apply(lambda r: float(np.var(copy_array(r))), axis=1)
    df_gc_sub["copy_sd"] = df_gc_sub.apply(lambda r: float(np.std(copy_array(r))), axis=1)
    count_cols = {s: f"{s}_count" for s in all_species}
    df_gc_sub = df_gc_sub.rename(columns=count_cols)

    # --- Per-OG immune enrichment / entropy / Fisher test ---
    print(">>> Computing enrichment and entropy for candidate OGs...")
    rows = []
    df_og_sub = df_og[df_og["Orthogroup"].astype(str).isin(candidate_ogs)].copy()
    for _, row in df_og_sub.iterrows():
        og = str(row["Orthogroup"])
        ref_ids, immune_counts = {}, {}
        for species in reference_species:
            ids_clean = [clean_protein_id(x) for x in split_ids(row.get(species, ""))]
            ref_ids[species] = ids_clean
            immune_counts[species] = sum(1 for x in ids_clean if x in immune_protein_set)

        ref_total_genes = sum(len(ref_ids[s]) for s in reference_species)
        immune_genes = sum(immune_counts[s] for s in reference_species)
        non_immune_genes = ref_total_genes - immune_genes

        immune_fraction = safe_div(immune_genes, ref_total_genes)
        immune_enrichment = safe_div(immune_fraction, background_fraction)

        a = immune_genes
        b = non_immune_genes
        c = total_immune_genes_background - a
        d = (total_ref_genes_background - total_immune_genes_background) - b
        if min(a, b, c, d) < 0:
            p_value, odds_ratio = np.nan, np.nan
        else:
            p_value, odds_ratio = fisher_exact_right_tail(a, b, c, d)

        entropy = compute_shannon_entropy([immune_counts[s] for s in reference_species], n_ref)

        record = {
            "Orthogroup": og,
            "ref_total_genes": ref_total_genes,
            "immune_genes": immune_genes,
            "non_immune_genes": non_immune_genes,
            "immune_fraction": immune_fraction,
            "background_fraction": background_fraction,
            "immune_enrichment": immune_enrichment,
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "immune_entropy": entropy,
        }
        for species in reference_species:
            record[f"immune_{species}"] = immune_counts[species]
        rows.append(record)

    df_stats = pd.DataFrame(rows)
    df_stats["FDR"] = benjamini_hochberg(df_stats["p_value"].values)

    # --- Merge structural stats ---
    print(">>> Merging duplication and structural stats...")
    struct_cols = ["Orthogroup"] + [f"{s}_count" for s in all_species] + \
                  ["total_genes_all", "copy_variance", "copy_sd"]
    df_summary = df_stats.merge(df_gc_sub[struct_cols], on="Orthogroup", how="left")
    df_summary = df_summary.merge(dup_count_df, on="Orthogroup", how="left")
    df_summary["dup_count"] = df_summary["dup_count"].fillna(0).astype(int)
    df_summary["dup_density"] = df_summary.apply(
        lambda r: safe_div(r["dup_count"], r["total_genes_all"]), axis=1)

    immune_cols = [f"immune_{s}" for s in reference_species]
    count_out_cols = [f"{s}_count" for s in all_species]
    summary_cols = (["Orthogroup", "ref_total_genes", "immune_genes", "non_immune_genes",
                     "immune_fraction", "background_fraction", "immune_enrichment", "odds_ratio",
                     "p_value", "FDR"] + immune_cols + ["immune_entropy"] + count_out_cols +
                    ["total_genes_all", "copy_variance", "copy_sd", "dup_count", "dup_density"])
    df_summary = df_summary[summary_cols].sort_values("Orthogroup")
    df_summary.to_csv(out_summary, sep="\t", index=False)
    print(f"    Saved OG summary: {out_summary}")

    df_merged = df_filtered.merge(df_summary, on="Orthogroup", how="left")
    df_merged.to_csv(out_merged, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"    Saved merged table: {out_merged}")
    print(f"=== Done === candidate OGs: {df_summary['Orthogroup'].nunique()}, rows: {len(df_merged)}")


if __name__ == "__main__":
    main()
