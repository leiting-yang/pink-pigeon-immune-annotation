#!/usr/bin/env python3
"""
compute_orthology_stats.py
==========================
From the OrthoFinder Orthologues output, classify each target-species
transcript's orthology relationship (1:1, 1:N, N:1, N:M) to each reference
species and combine it with reference immune-gene status, then merge onto the
OG-stats table. Species-agnostic (driven by `species:` in the config).

Inputs : {target}__v__{species}.tsv (Orthologues), master_lookup_table.csv,
         Final_Filtered_List_with_OG_stats.csv
Output : Final_Filtered_List_with_OG_stats_and_orthology.csv
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species

ORTH_TYPE_PRIORITY = {"1:1": 0, "1:N": 1, "N:1": 2, "N:M": 3}


def clean_protein_id(pid):
    if pd.isna(pid):
        return ""
    pid = str(pid).strip()
    return pid.split(".")[0] if "." in pid else pid


def clean_pp_id(pid):
    if pd.isna(pid):
        return ""
    return str(pid).strip().replace("transcript:", "").replace("mRNA:", "").replace("transcript_", "")


def split_ids(cell):
    if pd.isna(cell):
        return []
    s = str(cell).strip()
    return [x.strip() for x in s.split(",") if x.strip()] if s else []


def classify_orthology(pp_count, ref_count):
    if pp_count == 1 and ref_count == 1:
        return "1:1"
    if pp_count == 1 and ref_count > 1:
        return "1:N"
    if pp_count > 1 and ref_count == 1:
        return "N:1"
    return "N:M"


def pick_best_row(rows):
    def sort_key(r):
        priority = ORTH_TYPE_PRIORITY.get(r["orth_type"], 99)
        immune_frac = r["ref_immune"] / r["ref_total"] if r["ref_total"] > 0 else 0
        return (priority, -immune_frac, -r["ref_immune"])
    return sorted(rows, key=sort_key)[0]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--orthologues-dir", help="Directory holding {target}__v__*.tsv")
    p.add_argument("--master", help="master_lookup_table.csv")
    p.add_argument("--filtered", help="..._with_OG_stats.csv")
    p.add_argument("--output", help="Output CSV with orthology columns")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    reference_species = sp["ref_names"]
    target = sp["target"]
    target_id_col = f"{target}_ProteinID"

    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    orthologues_dir = args.orthologues_dir or wd(orth.get("orthologues_dir",
                                                          f"Orthologues/Orthologues_{target}"))
    master_csv = args.master or wd(orth.get("master_lookup", "master_lookup_table.csv"))
    filtered_csv = args.filtered or wd(orth.get("filtered_with_og_stats",
                                               "Final_Filtered_List_with_OG_stats.csv"))
    output_csv = args.output or wd(orth.get("filtered_with_orthology",
                                           "Final_Filtered_List_with_OG_stats_and_orthology.csv"))

    orthologues_files = {
        species: os.path.join(orthologues_dir, f"{target}__v__{species}.tsv")
        for species in reference_species
    }

    print(">>> Loading reference immune-protein list...")
    df_master = pd.read_csv(master_csv, encoding="utf-8-sig")
    df_master["ProteinID_Clean"] = df_master["ProteinID"].apply(clean_protein_id)
    immune_set = set(df_master["ProteinID_Clean"].dropna().astype(str))
    print(f"    immune proteins: {len(immune_set)}")

    print(">>> Loading candidate transcript list...")
    df_filtered = pd.read_csv(filtered_csv, encoding="utf-8-sig")
    candidate_pp_ids = set(df_filtered[target_id_col].dropna().astype(str))
    print(f"    candidate transcripts: {len(candidate_pp_ids)}")

    all_records = {}
    for species, filepath in orthologues_files.items():
        print(f">>> Parsing Orthologues: {species} ({filepath})")
        df_orth = pd.read_csv(filepath, sep="\t", dtype=str)
        if target not in df_orth.columns or species not in df_orth.columns:
            print(f"    [ERROR] expected columns Orthogroup, {target}, {species}; "
                  f"got {list(df_orth.columns)}")
            continue

        matched = skipped = 0
        for _, row in df_orth.iterrows():
            og = str(row.get("Orthogroup", ""))
            pp_clean = [clean_pp_id(x) for x in split_ids(row[target])]
            ref_clean = [clean_protein_id(x) for x in split_ids(row[species])]
            pp_count, ref_count = len(pp_clean), len(ref_clean)
            if pp_count == 0 or ref_count == 0:
                skipped += 1
                continue

            orth_type = classify_orthology(pp_count, ref_count)
            ref_immune_count = sum(1 for rid in ref_clean if rid in immune_set)
            all_immune = (ref_immune_count == ref_count) and (ref_count > 0)

            for pp_id in pp_clean:
                if pp_id not in candidate_pp_ids:
                    continue
                matched += 1
                all_records.setdefault(pp_id, {}).setdefault(species, []).append({
                    "orthogroup": og,
                    "orth_type": orth_type,
                    "pp_count": pp_count,
                    "ref_total": ref_count,
                    "ref_immune": ref_immune_count,
                    "all_immune": all_immune,
                    "ref_ids_str": ", ".join(ref_clean),
                })
        print(f"    skipped empty rows: {skipped}; matched candidate records: {matched}")

    print(">>> Summarizing orthology to transcript level...")
    summary_rows = []
    for pp_id in candidate_pp_ids:
        rec = all_records.get(pp_id, {})
        row_out = {target_id_col: pp_id}
        oneto1_immune_species = []
        onetoN_allimmune_species = []

        for species in reference_species:
            prefix = f"Orth_{species}"
            sp_rows = rec.get(species, [])
            if not sp_rows:
                row_out[f"{prefix}_type"] = np.nan
                row_out[f"{prefix}_ref_total"] = np.nan
                row_out[f"{prefix}_ref_immune"] = np.nan
                row_out[f"{prefix}_all_immune"] = np.nan
                row_out[f"{prefix}_rows"] = 0
            else:
                row_out[f"{prefix}_rows"] = len(sp_rows)
                best = pick_best_row(sp_rows)
                row_out[f"{prefix}_type"] = best["orth_type"]
                row_out[f"{prefix}_ref_total"] = best["ref_total"]
                row_out[f"{prefix}_ref_immune"] = best["ref_immune"]
                row_out[f"{prefix}_all_immune"] = best["all_immune"]
                if best["orth_type"] == "1:1" and best["ref_immune"] >= 1:
                    oneto1_immune_species.append(species)
                if best["orth_type"] == "1:N" and best["all_immune"]:
                    onetoN_allimmune_species.append(species)

        row_out["Orth_1to1_immune_count"] = len(oneto1_immune_species)
        row_out["Orth_1to1_immune_species"] = "; ".join(oneto1_immune_species) or np.nan
        row_out["Orth_1toN_allimmune_count"] = len(onetoN_allimmune_species)
        row_out["Orth_1toN_allimmune_species"] = "; ".join(onetoN_allimmune_species) or np.nan
        summary_rows.append(row_out)

    df_summary = pd.DataFrame(summary_rows)
    print(f"    summary rows: {len(df_summary)}")

    df_merged = df_filtered.merge(df_summary, on=target_id_col, how="left")
    df_merged.to_csv(output_csv, index=False, encoding="utf-8-sig", quoting=csv.QUOTE_NONNUMERIC)
    print(f">>> Saved: {output_csv}  ({len(df_merged)} rows)")


if __name__ == "__main__":
    main()
