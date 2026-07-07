#!/usr/bin/env python3
"""
process_symbols.py
==================
Predict a target-species gene symbol from the reference symbols and
systematically compare it against the native symbol from the genome annotation.
Species-agnostic: the reference symbol sources, the consensus group and the
fallback priority all come from `species:` in the config.

Prediction order: consensus (intersection of `consensus_group`) > each species
in `prediction_priority` > Kofam KO. Comparison categories: Match,
Alternative_Match, Group_Match, Broad_Match, Paralog_Likely, Mismatch,
Novel_Annotation, Not_Predictable, No_Info.

Input : Immune_Annotated_Final.csv
Output: Immune_Predict_Result_Final.csv
"""

import argparse
import csv
import difflib
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species, symbol_columns


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
    if pd.isna(gene_str) or str(gene_str).strip() == "":
        return []
    cleaned = []
    for p in re.split(r"[;,]\s*", str(gene_str)):
        p = p.strip()
        if p and p.lower() != "nan":
            cleaned.append(p)
    return cleaned


def get_all_source_symbols(row, symbol_sources):
    """Union of all reference symbols for a gene (original case), '/'-joined."""
    all_symbols = set()
    for source in symbol_sources:
        for gene in clean_split_genes(row.get(source)):
            all_symbols.add(gene)
    return "/".join(sorted(all_symbols)) if all_symbols else None


def make_get_prediction(sp):
    """Build a get_prediction(row) closure from the species config."""
    ref_names = sp["ref_names"]
    consensus = [s for s in sp["consensus_group"] if s in ref_names]
    priority = [s for s in sp["prediction_priority"] if s in ref_names]
    # References that validate an ambiguous consensus-species pick: everything
    # not in the consensus group, plus the Kofam KO symbols.
    non_consensus = [s for s in ref_names if s not in consensus]

    def get_prediction(row):
        lists = {s: clean_split_genes(row.get(f"{s}_GeneSymbols")) for s in ref_names}
        sets = {s: {x.upper() for x in lists[s]} for s in ref_names}
        ko_list = clean_split_genes(row.get("KO_Gene_Symbol"))
        set_ko = {x.upper() for x in ko_list}

        validation_pool = set_ko.union(*[sets[s] for s in non_consensus]) if non_consensus else set_ko

        predict_list, source = [], ""

        # 1. Consensus: intersection of the consensus-group species
        if len(consensus) >= 2:
            inter = set.intersection(*[sets[s] for s in consensus])
            if inter:
                for cs in consensus:
                    predict_list = [g for g in lists[cs] if g.upper() in inter]
                    if predict_list:
                        break
                source = "Consensus"

        # 2. Single-species fallback in priority order
        if not predict_list:
            for spn in priority:
                s_list = lists[spn]
                if not s_list:
                    continue
                if len(s_list) == 1:
                    predict_list, source = s_list, spn
                else:
                    valid = sets[spn] & validation_pool
                    predict_list = [g for g in s_list if g.upper() in valid] if valid else s_list
                    source = spn
                break

        # 3. Kofam KO fallback
        if not predict_list and ko_list:
            predict_list, source = ko_list, "Kofam"

        if not predict_list:
            return None, None, False

        seen = set()
        final_list = [x for x in predict_list if not (x.upper() in seen or seen.add(x.upper()))]
        is_ambiguous = len(final_list) > 1
        if is_ambiguous:
            source = f"{source} (Ambiguous)"
        return "/".join(final_list), source, is_ambiguous

    return get_prediction


def sophisticated_compare_final(row):
    """Compare predicted symbol against the native symbol; return a category."""
    native = str(row["Native_Symbol"]).strip()
    predict = str(row["Predicted_Symbol"]).strip()
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

    if native_upper in predict_candidates:
        return "Match"
    if native_upper in all_candidates:
        return "Alternative_Match"

    best_fuzzy_status = "Mismatch"
    for cand in predict_candidates:
        if "_" in cand:
            if native_upper in [x.upper() for x in expand_kegg_symbol(cand)]:
                return "Group_Match"
        if best_fuzzy_status in ("Mismatch", "Paralog_Likely"):
            if (cand.startswith(native_upper) or native_upper.startswith(cand)) \
                    and len(cand) > 3 and len(native_upper) > 3:
                best_fuzzy_status = "Broad_Match"
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
    p.add_argument("--input", help="Immune_Annotated_Final.csv")
    p.add_argument("--output", help="Immune_Predict_Result_Final.csv")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    symbol_sources = symbol_columns(sp) + ["KO_Gene_Symbol"]
    get_prediction = make_get_prediction(sp)

    tier = cfg.get("tiering", {})
    work = tier.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    input_file = args.input or wd(tier.get("annotated_final", "Immune_Annotated_Final.csv"))
    output_file = args.output or wd(tier.get("predict_result", "Immune_Predict_Result_Final.csv"))

    print(f"Reading {input_file} ...")
    df = pd.read_csv(input_file)

    print("1. building source pool...")
    df["All_Predict_Symbols"] = df.apply(lambda r: get_all_source_symbols(r, symbol_sources), axis=1)

    print("2. predicting symbols...")
    pred = df.apply(get_prediction, axis=1)
    df["Predicted_Symbol"] = [x[0] for x in pred]
    df["Predict_Sources"] = [x[1] for x in pred]
    df["Ambiguity_Flag"] = [x[2] for x in pred]

    print("3. comparing against native symbols...")
    df["Symbol_Comparison"] = df.apply(sophisticated_compare_final, axis=1)

    df.to_csv(output_file, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f"Done. Saved to {output_file}")


if __name__ == "__main__":
    main()
