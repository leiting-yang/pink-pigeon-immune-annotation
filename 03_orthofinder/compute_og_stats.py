#!/usr/bin/env python3
"""
compute_og_stats.py
===================
Per-orthogroup immune-enrichment and structural statistics for the candidate
orthogroups.

For each candidate OG it computes, from the reference species only:
  - immune fraction and enrichment vs a genome-wide background
  - one-sided Fisher exact test (hypergeometric, no scipy) + Benjamini-Hochberg FDR
  - normalized Shannon entropy of immune counts across the three species
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

import numpy as np
import pandas as pd

REFERENCE_SPECIES = ["Mouse", "Chicken", "ZebraFinch"]


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


# ------------------------- statistics utilities ------------------------------
def clean_protein_id(pid):
    """Remove version suffix, e.g. ENSMUSP00000008830.9 -> ENSMUSP00000008830."""
    if pd.isna(pid):
        return ""
    pid = str(pid).strip()
    if pid == "":
        return ""
    return pid.split(".")[0] if "." in pid else pid


def split_ids(cell):
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    if s == "":
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def safe_div(a, b):
    return a / b if b not in (0, 0.0) else np.nan


def compute_shannon_entropy(counts):
    """Normalized Shannon entropy H/log(3); NaN when total is 0."""
    total = sum(counts)
    if total == 0:
        return np.nan
    probs = [c / total for c in counts if c > 0]
    H = -sum(p * math.log(p) for p in probs)
    return H / math.log(3)


def fisher_exact_right_tail(a, b, c, d):
    """One-sided (enrichment) Fisher exact test via the hypergeometric tail."""
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
    p.add_argument("--out-merged", help="PinkPigeon_Final_Filtered_List_with_OG_stats.csv")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
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
        for sp in REFERENCE_SPECIES:
            ids_clean = [clean_protein_id(x) for x in split_ids(row.get(sp, ""))]
            total_ref_genes_background += len(ids_clean)
            total_immune_genes_background += sum(1 for x in ids_clean if x in immune_protein_set)
    background_fraction = safe_div(total_immune_genes_background, total_ref_genes_background)
    print(f"    Background immune fraction: {background_fraction:.6f}")

    # --- Duplication counts per OG ---
    dup_count_df = df_dup.groupby("Orthogroup").size().reset_index(name="dup_count")

    # --- Gene-copy counts (candidate OGs only) ---
    needed = ["Orthogroup", "Chicken", "Mouse", "PinkPigeon", "ZebraFinch", "Total"]
    missing = [c for c in needed if c not in df_gc.columns]
    if missing:
        raise ValueError(f"Missing columns in GeneCount.tsv: {missing}")

    df_gc_sub = df_gc[df_gc["Orthogroup"].astype(str).isin(candidate_ogs)].copy()
    for col in ("Chicken", "Mouse", "PinkPigeon", "ZebraFinch", "Total"):
        df_gc_sub[col] = pd.to_numeric(df_gc_sub[col], errors="coerce").fillna(0).astype(int)
    df_gc_sub["total_genes_all"] = (
        df_gc_sub["PinkPigeon"] + df_gc_sub["Mouse"] + df_gc_sub["Chicken"] + df_gc_sub["ZebraFinch"])

    def copy_array(row):
        return np.array([row["PinkPigeon"], row["Mouse"], row["Chicken"], row["ZebraFinch"]], dtype=float)

    df_gc_sub["copy_variance"] = df_gc_sub.apply(lambda r: float(np.var(copy_array(r))), axis=1)
    df_gc_sub["copy_sd"] = df_gc_sub.apply(lambda r: float(np.std(copy_array(r))), axis=1)
    df_gc_sub = df_gc_sub.rename(columns={
        "PinkPigeon": "pink_count", "Mouse": "mouse_count",
        "Chicken": "chicken_count", "ZebraFinch": "zebra_count"})

    # --- Per-OG immune enrichment / entropy / Fisher test ---
    print(">>> Computing enrichment and entropy for candidate OGs...")
    rows = []
    df_og_sub = df_og[df_og["Orthogroup"].astype(str).isin(candidate_ogs)].copy()
    for _, row in df_og_sub.iterrows():
        og = str(row["Orthogroup"])
        ref_ids, immune_counts = {}, {}
        for sp in REFERENCE_SPECIES:
            ids_clean = [clean_protein_id(x) for x in split_ids(row.get(sp, ""))]
            ref_ids[sp] = ids_clean
            immune_counts[sp] = sum(1 for x in ids_clean if x in immune_protein_set)

        ref_total_genes = sum(len(ref_ids[sp]) for sp in REFERENCE_SPECIES)
        immune_genes = sum(immune_counts[sp] for sp in REFERENCE_SPECIES)
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

        immune_entropy = compute_shannon_entropy(
            [immune_counts["Mouse"], immune_counts["Chicken"], immune_counts["ZebraFinch"]])

        rows.append({
            "Orthogroup": og,
            "ref_total_genes": ref_total_genes,
            "immune_genes": immune_genes,
            "non_immune_genes": non_immune_genes,
            "immune_fraction": immune_fraction,
            "background_fraction": background_fraction,
            "immune_enrichment": immune_enrichment,
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "immune_mouse": immune_counts["Mouse"],
            "immune_chicken": immune_counts["Chicken"],
            "immune_zebra": immune_counts["ZebraFinch"],
            "immune_entropy": immune_entropy,
        })

    df_stats = pd.DataFrame(rows)
    df_stats["FDR"] = benjamini_hochberg(df_stats["p_value"].values)

    # --- Merge structural stats ---
    print(">>> Merging duplication and structural stats...")
    df_summary = df_stats.merge(
        df_gc_sub[["Orthogroup", "pink_count", "mouse_count", "chicken_count",
                   "zebra_count", "total_genes_all", "copy_variance", "copy_sd"]],
        on="Orthogroup", how="left")
    df_summary = df_summary.merge(dup_count_df, on="Orthogroup", how="left")
    df_summary["dup_count"] = df_summary["dup_count"].fillna(0).astype(int)
    df_summary["dup_density"] = df_summary.apply(
        lambda r: safe_div(r["dup_count"], r["total_genes_all"]), axis=1)

    summary_cols = [
        "Orthogroup", "ref_total_genes", "immune_genes", "non_immune_genes",
        "immune_fraction", "background_fraction", "immune_enrichment", "odds_ratio",
        "p_value", "FDR", "immune_mouse", "immune_chicken", "immune_zebra",
        "immune_entropy", "pink_count", "mouse_count", "chicken_count", "zebra_count",
        "total_genes_all", "copy_variance", "copy_sd", "dup_count", "dup_density"]
    df_summary = df_summary[summary_cols].sort_values("Orthogroup")
    df_summary.to_csv(out_summary, sep="\t", index=False)
    print(f"    Saved OG summary: {out_summary}")

    # --- Merge back to the filtered Pink Pigeon table ---
    df_merged = df_filtered.merge(df_summary, on="Orthogroup", how="left")
    df_merged.to_csv(out_merged, index=False, encoding="utf-8-sig",
                     quoting=csv.QUOTE_NONNUMERIC)
    print(f"    Saved merged table: {out_merged}")
    print(f"=== Done === candidate OGs: {df_summary['Orthogroup'].nunique()}, "
          f"rows: {len(df_merged)}")


if __name__ == "__main__":
    main()
