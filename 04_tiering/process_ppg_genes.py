#!/usr/bin/env python3
"""
process_ppg_genes.py
====================
Predict a Pink Pigeon gene symbol from the reference symbols and systematically
compare it against the native symbol produced by the genome annotation.

Prediction priority: avian consensus (Chicken n ZebraFinch) > Chicken >
ZebraFinch > Mouse > Kofam. Comparison categories: Match, Alternative_Match,
Group_Match, Broad_Match, Paralog_Likely, Mismatch, Novel_Annotation,
Not_Predictable, No_Info.

Input : PinkPigeon_Immune_Annotated_Final.csv
Output: PinkPigeon_Immune_Predict_Result_Final.csv
        (adds All_Predict_Symbols, PPG_Predict_Symbol, Predict_Sources,
         Ambiguity_Flag, Symbol_Comparison)
"""

import argparse
import csv
import difflib
import os
import re

import pandas as pd


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def expand_kegg_symbol(symbol):
    """Expand a grouped symbol like 'GRK4_5_6' into ['GRK4', 'GRK5', 'GRK6']."""
    if "_" not in symbol:
        return [symbol]
    parts = symbol.split("_")
    match = re.match(r"^([A-Za-z0-9\-]+?)(\d+)$", parts[0])
    prefix = match.group(1) if match else ""
    expanded = []
    for i, p in enumerate(parts):
        p = p.strip()
        if i == 0:
            expanded.append(p)
        elif p.isdigit() and prefix:
            expanded.append(prefix + p)
        else:
            expanded.append(p)
    return expanded


def clean_split_genes(gene_str):
    """Split a multi-symbol string into a clean list."""
    if pd.isna(gene_str) or str(gene_str).strip() == "":
        return []
    cleaned = []
    for p in re.split(r"[;,]\s*", str(gene_str)):
        p = p.strip()
        if p and p.lower() != "nan":
            cleaned.append(p)
    return cleaned


def get_all_source_symbols(row):
    """Union of all reference symbols for a gene (original case), '/'-joined."""
    all_symbols = set()
    for source in ("Chicken_GeneSymbols", "ZebraFinch_GeneSymbols",
                   "Mouse_GeneSymbols", "KO_Gene_Symbol"):
        for gene in clean_split_genes(row[source]):
            all_symbols.add(gene)
    return "/".join(sorted(all_symbols)) if all_symbols else None


def get_prediction(row):
    """Predict a symbol with the documented source priority."""
    s_ck = clean_split_genes(row["Chicken_GeneSymbols"])
    s_zf = clean_split_genes(row["ZebraFinch_GeneSymbols"])
    s_mm = clean_split_genes(row["Mouse_GeneSymbols"])
    s_ko = clean_split_genes(row["KO_Gene_Symbol"])

    set_ck = {x.upper() for x in s_ck}
    set_zf = {x.upper() for x in s_zf}
    set_mm = {x.upper() for x in s_mm}
    set_ko = {x.upper() for x in s_ko}
    validation_pool = set_mm | set_ko

    predict_list = []
    source = ""
    avian_intersect = set_ck & set_zf

    if avian_intersect:
        predict_list = [g for g in s_ck if g.upper() in avian_intersect]
        if not predict_list:
            predict_list = [g for g in s_zf if g.upper() in avian_intersect]
        source = "Avian (Consensus)"
    elif s_ck:
        if len(s_ck) == 1:
            predict_list, source = s_ck, "Chicken"
        else:
            valid = set_ck & validation_pool
            predict_list = [g for g in s_ck if g.upper() in valid] if valid else s_ck
            source = "Chicken"
    elif s_zf:
        if len(s_zf) == 1:
            predict_list, source = s_zf, "ZebraFinch"
        else:
            valid = set_zf & validation_pool
            predict_list = [g for g in s_zf if g.upper() in valid] if valid else s_zf
            source = "ZebraFinch"
    elif s_mm:
        predict_list, source = s_mm, "Mouse"
    elif s_ko:
        predict_list, source = s_ko, "Kofam"

    if not predict_list:
        return None, None, False

    seen = set()
    final_list = [x for x in predict_list if not (x.upper() in seen or seen.add(x.upper()))]
    is_ambiguous = len(final_list) > 1
    if is_ambiguous:
        source = f"{source} (Ambiguous)"
    return "/".join(final_list), source, is_ambiguous


def sophisticated_compare_final(row):
    """Compare predicted symbol against the native symbol; return a category."""
    native = str(row["PPG_Native_Symbol"]).strip()
    predict = str(row["PPG_Predict_Symbol"]).strip()
    all_potential = str(row["All_Predict_Symbols"]).strip()

    native_empty = native.lower() == "nan" or native == ""
    predict_empty = predict.lower() == "nan" or predict == "None" or predict == ""

    if native_empty and predict_empty:
        return "No_Info"
    if not native_empty and predict_empty:
        return "Not_Predictable"
    if native_empty and not predict_empty:
        return "Novel_Annotation"

    native_upper = native.upper()
    predict_candidates = [p.strip().upper() for p in predict.split("/")]
    all_candidates = ([p.strip().upper() for p in all_potential.split("/")]
                      if all_potential and all_potential.lower() != "nan" else [])

    # Priority 1: exact hit in the prediction
    if native_upper in predict_candidates:
        return "Match"
    # Priority 2: exact hit somewhere in the source pool
    if native_upper in all_candidates:
        return "Alternative_Match"

    # Priority 3: fuzzy matches
    best_fuzzy_status = "Mismatch"
    for cand in predict_candidates:
        # A. grouped family match, e.g. cand='GRK4_5_6' native='GRK5'
        if "_" in cand:
            if native_upper in [x.upper() for x in expand_kegg_symbol(cand)]:
                return "Group_Match"
        # B. broad / prefix containment
        if best_fuzzy_status in ("Mismatch", "Paralog_Likely"):
            if (cand.startswith(native_upper) or native_upper.startswith(cand)) \
                    and len(cand) > 3 and len(native_upper) > 3:
                best_fuzzy_status = "Broad_Match"
        # C. paralog likely (similarity), lowest priority
        if best_fuzzy_status == "Mismatch":
            ratio = difflib.SequenceMatcher(None, cand, native_upper).ratio()
            is_high_sim = ratio > 0.85
            is_prefix_sim = (cand[:3] == native_upper[:3] and len(cand) > 3
                             and abs(len(cand) - len(native_upper)) <= 1)
            if is_high_sim or is_prefix_sim:
                best_fuzzy_status = "Paralog_Likely"
    return best_fuzzy_status


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--input", help="PinkPigeon_Immune_Annotated_Final.csv")
    p.add_argument("--output", help="PinkPigeon_Immune_Predict_Result_Final.csv")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    tier = cfg.get("tiering", {})
    work = tier.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    input_file = args.input or wd(tier.get("annotated_final", "PinkPigeon_Immune_Annotated_Final.csv"))
    output_file = args.output or wd(tier.get("predict_result", "PinkPigeon_Immune_Predict_Result_Final.csv"))

    print(f"Reading {input_file} ...")
    df = pd.read_csv(input_file)

    print("1. building source pool...")
    df["All_Predict_Symbols"] = df.apply(get_all_source_symbols, axis=1)

    print("2. predicting symbols...")
    pred = df.apply(get_prediction, axis=1)
    df["PPG_Predict_Symbol"] = [x[0] for x in pred]
    df["Predict_Sources"] = [x[1] for x in pred]
    df["Ambiguity_Flag"] = [x[2] for x in pred]

    print("3. comparing against native symbols...")
    df["Symbol_Comparison"] = df.apply(sophisticated_compare_final, axis=1)

    df.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved to {output_file}")


if __name__ == "__main__":
    main()
