#!/usr/bin/env python3
"""
data_merge_new.py
=================
Build the reference immune-protein lookup table from three species' BioMart
exports plus a curated mouse list.

Steps:
  1. clean text fields (newlines, quotes, odd punctuation) to protect the CSV
  2. left-join mouse BioMart with the curated list (functional annotation)
  3. tag each mouse entry's source (Curated_List vs BioMart_GO)
  4. concatenate all three species into one master table

Outputs:
  - reference_mouse_final.csv  (mouse BioMart + curated)
  - master_lookup_table.csv    (all species; core input to OrthoFinder filtering)
"""

import argparse
import csv
import os
import re

import pandas as pd

COMMON_COLS = [
    "ProteinID", "GeneID", "Species", "GeneSymbol", "Description", "GO_Terms",
    "Category1", "Subcategory", "UniProt_function", "Immune_Source",
]


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def clean_text_field(text):
    """Flatten a text field so it cannot break the CSV layout."""
    if pd.isna(text):
        return ""
    text = str(text)
    text = text.replace("\n", " ").replace("\r", " ")
    text = text.replace('"', "'")
    text = text.replace("’", "'").replace("‘", "'")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--chicken", help="processed chicken BioMart CSV")
    p.add_argument("--zebrafinch", help="processed zebrafinch BioMart CSV")
    p.add_argument("--mouse", help="processed mouse BioMart CSV")
    p.add_argument("--curated", help="curated mouse immune-gene CSV")
    p.add_argument("--out-mouse", help="output reference_mouse_final.csv")
    p.add_argument("--out-master", help="output master_lookup_table.csv")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    ext = cfg.get("external_data", {})
    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    chicken_f = args.chicken or ext.get("chicken_biomart", "processed_chicken_biomart.csv")
    zebra_f = args.zebrafinch or ext.get("zebrafinch_biomart", "processed_zebrafinch_biomart.csv")
    mouse_f = args.mouse or ext.get("mouse_biomart", "processed_mouse_biomart.csv")
    curated_f = args.curated or ext.get("curated_list", "ImmuneGeneFunction_20240520.csv")
    out_mouse = args.out_mouse or wd(orth.get("reference_mouse", "reference_mouse_final.csv"))
    out_master = args.out_master or wd(orth.get("master_lookup", "master_lookup_table.csv"))

    print("Loading input tables...")
    df_chicken = pd.read_csv(chicken_f)
    df_zebrafinch = pd.read_csv(zebra_f)
    df_mouse = pd.read_csv(mouse_f)
    try:
        df_curated = pd.read_csv(curated_f, encoding="utf-8")
    except UnicodeDecodeError:
        print("UTF-8 read failed for curated list; retrying with cp1252...")
        df_curated = pd.read_csv(curated_f, encoding="cp1252")

    # --- Clean long text fields ---
    print("Cleaning text fields...")
    for col in ("UniProt_function", "Description", "Comments", "GeneName"):
        if col in df_curated.columns:
            df_curated[col] = df_curated[col].apply(clean_text_field)
    for df in (df_mouse, df_chicken, df_zebrafinch):
        if "Description" in df.columns:
            df["Description"] = df["Description"].apply(clean_text_field)

    # --- Merge mouse BioMart with curated list ---
    print("Merging mouse BioMart with curated list...")
    df_mouse["GeneID"] = df_mouse["GeneID"].astype(str).str.strip()
    df_curated["GeneStableID"] = df_curated["GeneStableID"].astype(str).str.strip()

    merge_cols = [c for c in ("GeneStableID", "Category1", "Subcategory", "UniProt_function")
                  if c in df_curated.columns]
    df_mouse_final = pd.merge(
        df_mouse, df_curated[merge_cols],
        left_on="GeneID", right_on="GeneStableID", how="left",
    )
    if "GeneStableID" in df_mouse_final.columns:
        df_mouse_final.drop(columns=["GeneStableID"], inplace=True)

    df_mouse_final["Immune_Source"] = df_mouse_final["Category1"].apply(
        lambda x: "Curated_List" if pd.notnull(x) else "BioMart_GO"
    )
    df_mouse_final.to_csv(out_mouse, index=False, encoding="utf-8-sig",
                          quoting=csv.QUOTE_NONNUMERIC)
    print(f"Saved: {out_mouse}")

    # --- Build the all-species master table ---
    print("Building master lookup table...")
    for col in ("Category1", "Subcategory", "UniProt_function", "Immune_Source"):
        df_chicken[col] = None
        df_zebrafinch[col] = None

    df_chicken = df_chicken.reindex(columns=COMMON_COLS)
    df_zebrafinch = df_zebrafinch.reindex(columns=COMMON_COLS)
    df_mouse_final = df_mouse_final.reindex(columns=COMMON_COLS)

    master_df = pd.concat([df_mouse_final, df_chicken, df_zebrafinch], ignore_index=True)
    master_df.to_csv(out_master, index=False, encoding="utf-8-sig",
                     quoting=csv.QUOTE_NONNUMERIC)
    print(f"Saved: {out_master}")


if __name__ == "__main__":
    main()
