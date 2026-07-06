#!/usr/bin/env python3
"""
data_merge_new.py
=================
Build the reference immune-protein lookup table from every reference species'
BioMart export (plus an optional curated functional list per species).

Species-agnostic: the reference species, their BioMart files and any curated
lists all come from `species:` in the config.

Steps per reference species:
  1. read its BioMart CSV and clean text fields
  2. if it has a curated_list, left-join it (Category1 / Subcategory /
     UniProt_function) and tag Immune_Source (Curated_List vs BioMart_GO)
Then concatenate all species into one master table.

Outputs:
  - reference_<species>_final.csv  (per curated species)
  - master_lookup_table.csv        (all species; core input to OG filtering)
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species

COMMON_COLS = [
    "ProteinID", "GeneID", "Species", "GeneSymbol", "Description", "GO_Terms",
    "Category1", "Subcategory", "UniProt_function", "Immune_Source",
]
CURATED_COLS = ["GeneStableID", "Category1", "Subcategory", "UniProt_function"]


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
    p.add_argument("--out-master", help="output master_lookup_table.csv")
    return p.parse_args()


def read_reference_species(entry, base_dir):
    """Read one reference species' BioMart table (+ optional curated list)."""
    name = entry["name"]
    biomart_path = resolve(entry.get("biomart"), base_dir)
    print(f"  - {name}: {biomart_path}")
    df = pd.read_csv(biomart_path)
    if "Description" in df.columns:
        df["Description"] = df["Description"].apply(clean_text_field)

    curated = entry.get("curated_list")
    if curated:
        curated_path = resolve(curated, base_dir)
        print(f"      curated list: {curated_path}")
        try:
            df_curated = pd.read_csv(curated_path, encoding="utf-8")
        except UnicodeDecodeError:
            df_curated = pd.read_csv(curated_path, encoding="cp1252")
        for col in ("UniProt_function", "Description", "Comments", "GeneName"):
            if col in df_curated.columns:
                df_curated[col] = df_curated[col].apply(clean_text_field)

        df["GeneID"] = df["GeneID"].astype(str).str.strip()
        df_curated["GeneStableID"] = df_curated["GeneStableID"].astype(str).str.strip()
        merge_cols = [c for c in CURATED_COLS if c in df_curated.columns]
        df = pd.merge(df, df_curated[merge_cols], left_on="GeneID",
                      right_on="GeneStableID", how="left")
        if "GeneStableID" in df.columns:
            df.drop(columns=["GeneStableID"], inplace=True)
        df["Immune_Source"] = df["Category1"].apply(
            lambda x: "Curated_List" if pd.notnull(x) else "BioMart_GO")

    df["Species"] = name
    return df


def resolve(path, base_dir):
    """Resolve a possibly-relative data path against base_dir."""
    if not path:
        return path
    return path if os.path.isabs(path) or os.path.exists(path) else os.path.join(base_dir, path)


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    out_master = args.out_master or (os.path.join(work, orth.get("master_lookup", "master_lookup_table.csv"))
                                     if work else orth.get("master_lookup", "master_lookup_table.csv"))

    print("Loading reference species...")
    species_frames = []
    for name in sp["ref_names"]:
        entry = sp["ref_map"][name]
        df = read_reference_species(entry, work)
        # Save a per-species reference file for the curated species.
        if entry.get("curated_list"):
            ref_out = os.path.join(work, f"reference_{name}_final.csv") if work else f"reference_{name}_final.csv"
            df.reindex(columns=COMMON_COLS).to_csv(
                ref_out, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
            print(f"      saved: {ref_out}")
        species_frames.append(df.reindex(columns=COMMON_COLS))

    print("Building master lookup table...")
    master_df = pd.concat(species_frames, ignore_index=True)
    master_df.to_csv(out_master, index=False, encoding="utf-8-sig",
                     quoting=csv.QUOTE_NONNUMERIC)
    print(f"Saved: {out_master}  ({len(master_df)} rows, "
          f"{len(sp['ref_names'])} reference species)")


if __name__ == "__main__":
    main()
