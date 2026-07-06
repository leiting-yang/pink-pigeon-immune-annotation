#!/usr/bin/env python3
"""
visualization.py
================
Generate the QC / summary figures from the final annotated table.
Species-agnostic: the functional-category plot uses the curated reference
species from `species:` in the config, and prediction sources are mapped
generically.

Input : PinkPigeon_Immune_Predict_Result_Final_with_OG_and_Orthology.csv
Output: 7 PNG figures (300 dpi) in the figure directory.
"""

import argparse
import os
import sys
from collections import Counter
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species

sns.set_theme(style="whitegrid")
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

TITLE_SIZE, LABEL_SIZE, TICK_SIZE, LEGEND_SIZE = 20, 18, 16, 18
plt.rcParams["axes.titlesize"] = TITLE_SIZE
plt.rcParams["axes.labelsize"] = LABEL_SIZE
plt.rcParams["xtick.labelsize"] = TICK_SIZE
plt.rcParams["ytick.labelsize"] = TICK_SIZE
plt.rcParams["legend.fontsize"] = LEGEND_SIZE

OUTDIR = "."


def save_plot(fig, filename):
    fig.tight_layout()
    out = os.path.join(OUTDIR, filename)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_venn(df):
    ortho_idx, kofam_idx, interpro_idx = set(), set(), set()
    for idx, row in df.iterrows():
        if pd.isna(row["Evidence_Sources"]):
            continue
        for s in str(row["Evidence_Sources"]).split(";"):
            s_clean = s.strip().lower()
            if "orthofinder" in s_clean:
                ortho_idx.add(idx)
            if "kofam" in s_clean:
                kofam_idx.add(idx)
            if "interpro" in s_clean:
                interpro_idx.add(idx)
    plt.figure(figsize=(10, 10))
    out = venn3([ortho_idx, kofam_idx, interpro_idx],
                set_labels=("Orthofinder", "KofamScan", "InterProScan"))
    for text in out.set_labels:
        if text:
            text.set_fontsize(LABEL_SIZE)
    for text in out.subset_labels:
        if text:
            text.set_fontsize(LABEL_SIZE)
    plt.title(f"Evidence Sources Intersection (Total Genes: {len(ortho_idx | kofam_idx | interpro_idx)})")
    save_plot(plt.gcf(), "1_Evidence_Venn_Diagram.png")


def plot_tier_distribution(df):
    if "Tier" not in df.columns:
        print("Column 'Tier' not found.")
        return
    plt.figure(figsize=(10, 8))
    tier_counts = df["Tier"].value_counts().sort_index()
    ax = sns.barplot(x=tier_counts.index, y=tier_counts.values, palette="viridis")
    plt.title(f"Tier Distribution (Total: {len(df)})")
    plt.xlabel("Tier")
    plt.ylabel("Gene Count")
    for i, v in enumerate(tier_counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=TICK_SIZE)
    save_plot(plt.gcf(), "2_Tier_Distribution.png")


def plot_immune_pathways(df):
    pathways = []
    for entry in df["Immune_Pathway"].dropna():
        for p in str(entry).split(";"):
            p = p.strip()
            p = re.sub(r"\s*\[.*?\]\s*$", "", p)
            p = re.sub(r"^\d+\s*", "", p)
            if p:
                pathways.append(p)
    if not pathways:
        print("No immune pathways found.")
        return
    counts = Counter(pathways).most_common(10)
    y_labels = [x[0] for x in counts]
    x_values = [x[1] for x in counts]
    plt.figure(figsize=(14, 10))
    ax = sns.barplot(x=x_values, y=y_labels, palette="magma")
    plt.title("Top 10 Immune Pathways (genes may map to multiple pathways)")
    plt.xlabel("Count")
    plt.ylabel("Pathway")
    for i, v in enumerate(x_values):
        ax.text(v + 0.5, i, str(v), va="center", fontsize=TICK_SIZE)
    save_plot(plt.gcf(), "3_Top10_Immune_Pathways.png")


def plot_functional_category(df, category_col):
    """Functional category of the curated reference species (single-category genes)."""
    if not category_col or category_col not in df.columns:
        print(f"[Plot 4] no curated category column ({category_col}); skipped.")
        return
    valid = df[category_col].dropna()
    multi = valid[valid.str.contains(";")]
    single = valid[~valid.str.contains(";")]
    print(f"[Plot 4] {len(multi)} genes have multiple categories; excluded.")
    cat_counts = single.value_counts()
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(x=cat_counts.values, y=cat_counts.index, palette="coolwarm")
    plt.title(f"Functional Category (Single Category Genes: {len(single)})")
    plt.xlabel("Count")
    plt.ylabel("Category")
    for i, v in enumerate(cat_counts.values):
        ax.text(v + 0.5, i, str(v), va="center", fontsize=TICK_SIZE)
    save_plot(plt.gcf(), "4_Functional_Category_Distribution.png")


def plot_symbol_comparison(df):
    def map_symbol(s):
        if pd.isna(s) or str(s) in ("No_Info", "Not_Predictable", "nan"):
            return "No Prediction"
        s_str = str(s)
        if s_str == "Match":
            return "Match"
        if s_str == "Novel_Annotation":
            return "Novel Annotation"
        if "Mismatch" in s_str:
            return "Mismatch"
        return "Other Match"

    df["Symbol_Group"] = df["Symbol_Comparison"].apply(map_symbol)
    counts = df["Symbol_Group"].value_counts()
    colors = sns.color_palette("Set2", len(counts))
    plt.figure(figsize=(11, 11))
    plt.pie(counts, labels=counts.index,
            autopct=lambda pct: f"{int(round(pct * sum(counts) / 100))}\n({pct:.1f}%)",
            textprops=dict(color="black", fontsize=LABEL_SIZE),
            startangle=140, pctdistance=0.8, colors=colors)
    plt.title(f"Symbol Comparison Distribution (Total: {len(df)})")
    save_plot(plt.gcf(), "5_Symbol_Comparison_Pie.png")
    return df


def plot_predict_sources(df):
    def map_source(s):
        if pd.isna(s):
            return "Not predictable"
        s = str(s).replace(" (Ambiguous)", "")
        if "Consensus" in s or "Avian" in s:
            return "Consensus"
        if s == "Kofam":
            return "KEGG"
        return s  # reference species name as-is

    counts = df["Predict_Sources"].apply(map_source).value_counts()
    plt.figure(figsize=(11, 11))
    plt.pie(counts, labels=counts.index,
            autopct=lambda pct: f"{int(round(pct * sum(counts) / 100))}\n({pct:.1f}%)",
            textprops=dict(color="black", fontsize=LABEL_SIZE),
            startangle=140, pctdistance=0.8, colors=sns.color_palette("pastel"))
    plt.title(f"Prediction Source Distribution (Total: {len(df)})")
    save_plot(plt.gcf(), "6_Prediction_Source_Pie.png")


def plot_ambiguity(df):
    if "Symbol_Group" not in df.columns:
        plot_symbol_comparison(df)
    predicted = df[df["Symbol_Group"] != "No Prediction"].copy()
    predicted["Ambiguity_Str"] = predicted["Ambiguity_Flag"].astype(str)

    ambig_counts = predicted["Ambiguity_Str"].value_counts()
    plt.figure(figsize=(8, 8))
    ax = sns.barplot(x=ambig_counts.index, y=ambig_counts.values, palette="Set3")
    plt.title(f"Ambiguity in Predicted Gene Symbols (Total: {len(predicted)})")
    plt.xlabel("Ambiguous Flag")
    plt.ylabel("Count")
    for i, v in enumerate(ambig_counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=TICK_SIZE)
    save_plot(plt.gcf(), "7a_Ambiguity_Overall.png")

    target_groups = ["Match", "Novel Annotation", "Mismatch"]
    subset = predicted[predicted["Symbol_Group"].isin(target_groups)]
    if len(subset) == 0:
        print("No data for the Match/Novel/Mismatch ambiguity plot.")
        return
    plt.figure(figsize=(12, 8))
    ax = sns.countplot(data=subset, x="Symbol_Group", hue="Ambiguity_Str",
                       palette="Set3", order=target_groups)
    plt.title("Predicted Gene Symbols Ambiguity Distribution by Category")
    plt.xlabel("Category")
    plt.ylabel("Count")
    plt.legend(title="Ambiguous", title_fontsize=LEGEND_SIZE, fontsize=TICK_SIZE)
    for container in ax.containers:
        ax.bar_label(container, fontsize=TICK_SIZE)
    save_plot(plt.gcf(), "7b_Ambiguity_By_Category.png")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--input", help="Final annotated table CSV")
    p.add_argument("--figure-dir", help="Directory for the output PNGs")
    return p.parse_args()


def main():
    global OUTDIR
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    tier = cfg.get("tiering", {})
    work = tier.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    input_file = args.input or wd(tier.get("final_table",
                                          "PinkPigeon_Immune_Predict_Result_Final_with_OG_and_Orthology.csv"))
    OUTDIR = args.figure_dir or wd(tier.get("figure_dir", "figures"))
    os.makedirs(OUTDIR, exist_ok=True)

    # Curated reference species (if any) drives the functional-category plot.
    curated_species = sp["curated"][0] if sp["curated"] else None
    category_col = f"{curated_species}_Category1" if curated_species else None

    print(f"Loading {input_file} ...")
    df = pd.read_csv(input_file)

    print("Generating plots...")
    plot_venn(df)
    plot_tier_distribution(df)
    plot_immune_pathways(df)
    plot_functional_category(df, category_col)
    df = plot_symbol_comparison(df)
    plot_predict_sources(df)
    plot_ambiguity(df)
    print(f"Done. All figures in {OUTDIR}")


if __name__ == "__main__":
    main()
