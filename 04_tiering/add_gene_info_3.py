#!/usr/bin/env python3
"""
add_gene_info_3.py
==================
Add the target-species native gene symbol to the master list, then (optionally)
enrich it with functional annotation from the conserved immune-gene database
(gene_info.csv), matched via an "all-in" symbol strategy over the native symbol,
every reference species' symbols and the Kofam KO symbol.

Set `species.gene_info_enrichment: false` in the config to skip the gene_info
step (e.g. for non-avian work where that database does not apply); the native
symbol is still added so downstream symbol comparison keeps working.

Inputs : PinkPigeon_Immune_Gene_Master_List.csv, gene_info.csv (optional),
         Nesoenas_genes_by_name.tsv (target native symbols)
Outputs: PinkPigeon_Immune_Annotated_Final.csv,
         Unmapped_Gene_Info_Symbols.txt  (only when enrichment runs)
"""

import argparse
import csv
import os
import sys
import warnings

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species, symbol_columns

warnings.simplefilter(action="ignore", category=FutureWarning)

KEY_COLUMN = "entrezgene_accession"


def standardize_symbol(symbol):
    if pd.isna(symbol) or str(symbol).strip().lower() == "nan":
        return ""
    return str(symbol).strip().upper()


def load_ppg_names(filepath):
    """Native TSV: column 0 = Gene Symbol, column 1 = Gene ID -> {id: SYMBOL}."""
    print(f"Loading native symbols from {filepath}...")
    id_to_symbol = {}
    df = pd.read_csv(filepath, sep="\t", header=None, dtype=str)
    for index, row in df.iterrows():
        if len(row) < 2:
            continue
        raw_symbol, raw_id = str(row[0]).strip(), str(row[1]).strip()
        if index == 0 and ("symbol" in raw_symbol.lower() or "id" in raw_id.lower()):
            continue  # skip header
        symbol = standardize_symbol(raw_symbol)
        if raw_id and symbol:
            id_to_symbol[raw_id] = symbol
    print(f"    mapped {len(id_to_symbol)} native gene IDs to symbols")
    return id_to_symbol


def lookup_function(row, info_dict, target_columns, symbol_sources):
    """Collect candidate symbols from all sources, then aggregate every hit."""
    candidates = []
    ppg_sym = str(row.get("PPG_Native_Symbol", "")).strip()
    if ppg_sym and ppg_sym.lower() != "nan":
        candidates.append(standardize_symbol(ppg_sym))
    for source in symbol_sources:
        for s in str(row.get(source, "")).split(";"):
            clean = standardize_symbol(s)
            if clean:
                candidates.append(clean)

    result_sets = {col: set() for col in target_columns}
    match_sources = set()
    for symbol in candidates:
        data = info_dict.get(symbol)
        if not data:
            continue
        match_sources.add(symbol)
        for col in target_columns:
            val = data.get(col)
            if pd.notna(val) and str(val).strip() != "" and str(val).lower() != "nan":
                result_sets[col].add(str(val))

    result = {col: "; ".join(sorted(vals)) for col, vals in result_sets.items()}
    result["Matched_Gene_Symbol"] = "; ".join(sorted(match_sources)) if match_sources else None
    return result


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--master", help="PinkPigeon_Immune_Gene_Master_List.csv")
    p.add_argument("--gene-info", help="gene_info.csv")
    p.add_argument("--ppg-names", help="target native symbol TSV")
    p.add_argument("--output", help="Annotated output CSV")
    p.add_argument("--unmapped", help="Unmapped symbols TXT")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    # Symbol sources: every reference species' symbols + the Kofam KO symbol.
    symbol_sources = symbol_columns(sp) + ["KO_Gene_Symbol"]

    tier = cfg.get("tiering", {})
    work = tier.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    master_file = args.master or wd(tier.get("master_list", "PinkPigeon_Immune_Gene_Master_List.csv"))
    gene_info_file = args.gene_info or cfg.get("external_data", {}).get("gene_info", "gene_info.csv")
    ppg_file = args.ppg_names or cfg.get("external_data", {}).get("ppg_native_symbols",
                                                                 "Nesoenas_genes_by_name.tsv")
    output_file = args.output or wd(tier.get("annotated_final", "PinkPigeon_Immune_Annotated_Final.csv"))
    unmapped_file = args.unmapped or wd(tier.get("unmapped_symbols", "Unmapped_Gene_Info_Symbols.txt"))

    print("Loading master list...")
    df_master = pd.read_csv(master_file)

    ppg_id_map = load_ppg_names(ppg_file)
    df_master["PPG_Native_Symbol"] = df_master["GeneID"].astype(str).str.strip().map(ppg_id_map)

    if not sp["gene_info_enrichment"]:
        print("gene_info enrichment disabled (species.gene_info_enrichment: false); "
              "writing table with native symbol only.")
        df_master.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
        print(f"Done. Saved to {output_file}")
        return

    print("Loading gene_info database...")
    try:
        df_info = pd.read_csv(gene_info_file)
    except Exception:
        print("CSV read failed; trying Excel (gene_info.xlsx)...")
        df_info = pd.read_excel("gene_info.xlsx")

    info_columns = [c for c in df_info.columns if c != KEY_COLUMN]
    print(f"Columns added from gene_info: {info_columns}")

    info_dict = {}
    all_info_symbols = set()
    for _, row in df_info.iterrows():
        key = standardize_symbol(row[KEY_COLUMN])
        if not key:
            continue
        info_dict[key] = row[info_columns].to_dict()
        all_info_symbols.add(key)
    print(f"Loaded {len(info_dict)} gene_info entries.")

    native_symbols = set(df_master["PPG_Native_Symbol"].dropna().apply(standardize_symbol))
    overlap = native_symbols & all_info_symbols
    print(f"Native symbols in master: {len(native_symbols)}; overlap with gene_info: {len(overlap)}")

    print("Mapping annotations...")
    mapped = df_master.apply(
        lambda r: lookup_function(r, info_dict, info_columns, symbol_sources), axis=1).tolist()
    df_mapped = pd.DataFrame(mapped)
    matched = df_mapped["Matched_Gene_Symbol"].notna().sum()
    print(f"Mapped annotations to {matched} / {len(df_master)} genes.")

    df_final = pd.concat([df_master, df_mapped], axis=1)

    used = set()
    for item in df_final["Matched_Gene_Symbol"].dropna():
        for s in str(item).split(";"):
            if s.strip():
                used.add(s.strip())
    with open(unmapped_file, "w") as fh:
        fh.write("Gene_Symbol\n")
        for sym in sorted(all_info_symbols - used):
            fh.write(f"{sym}\n")

    base_cols = list(df_master.columns)
    new_cols = ["Matched_Gene_Symbol"] + info_columns
    base_cols = [c for c in base_cols if c not in new_cols]
    insert_idx = base_cols.index("Tier") + 1 if "Tier" in base_cols else 1
    final_cols = base_cols[:insert_idx] + new_cols + base_cols[insert_idx:]
    final_cols = [c for c in final_cols if c in df_final.columns]
    df_final = df_final[final_cols]

    df_final.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved to {output_file}")


if __name__ == "__main__":
    main()
