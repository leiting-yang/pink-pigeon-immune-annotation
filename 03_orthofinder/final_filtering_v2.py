#!/usr/bin/env python3
"""
final_filtering_v2.py
=====================
"Exist-to-keep" filtering: for every orthogroup that contains at least one
reference-species immune protein, keep all target-species proteins in that group
and attach the reference annotation.

Species-agnostic: reference species, the target species and the per-species
output columns all come from `species:` in the config.

Inputs : Orthogroups.tsv, master_lookup_table.csv
Output : PinkPigeon_Final_Filtered_List.csv (or the target's filtered list)
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_common import load_config, get_species, species_field_map, ref_annotation_columns


def clean_protein_id(pid):
    """Remove the Ensembl version suffix (e.g. ENSMUSP...9 -> ENSMUSP...)."""
    if pd.isna(pid):
        return ""
    pid = str(pid).strip()
    return pid.split(".")[0] if "." in pid else pid


def clean_text_field(text):
    if pd.isna(text):
        return ""
    text = str(text)
    text = text.replace("\n", " ").replace("\r", " ")
    text = text.replace('"', "'")
    text = text.replace("’", "'").replace("‘", "'")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def split_cell(row, key):
    raw = row.get(key)
    if pd.isna(raw):
        return []
    return [p.strip() for p in str(raw).split(",") if p.strip()]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--orthogroups", help="Orthogroups.tsv")
    p.add_argument("--master", help="master_lookup_table.csv")
    p.add_argument("--output", help="Filtered list CSV")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    sp = get_species(cfg)
    target = sp["target"]
    ref_names = sp["ref_names"]
    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    og_file = args.orthogroups or wd(orth.get("orthogroups_tsv", "Orthogroups.tsv"))
    master_file = args.master or wd(orth.get("master_lookup", "master_lookup_table.csv"))
    output_file = args.output or wd(orth.get("filtered_list", "PinkPigeon_Final_Filtered_List.csv"))

    # Per-species: which master fields to aggregate and into which output column.
    field_map = species_field_map(sp)               # (species, master_field, out_col)
    per_species_fields = {name: [] for name in ref_names}
    for name, field, col in field_map:
        per_species_fields[name].append((field, col))

    print("Loading reference immune-protein lookup...")
    try:
        df_master = pd.read_csv(master_file, encoding="utf-8-sig")
    except UnicodeDecodeError:
        df_master = pd.read_csv(master_file, encoding="utf-8")
    df_master["ProteinID_Clean"] = df_master["ProteinID"].apply(clean_protein_id)
    master_dict = df_master.set_index("ProteinID_Clean").to_dict("index")
    immune_set = set(df_master["ProteinID_Clean"])
    print(f"    reference immune IDs: {len(immune_set)}")

    print("Reading Orthogroups.tsv...")
    df_og = pd.read_csv(og_file, sep="\t")
    missing = [c for c in ref_names + [target] if c not in df_og.columns]
    if missing:
        raise SystemExit(
            f"ERROR: Orthogroups.tsv has no column(s) {missing}. The column names "
            f"come from the input FASTA basenames; they must match the config "
            f"species names ({[target] + ref_names}).")

    def aggregate_species(name, id_list):
        """Aggregate one species' annotation into its output columns."""
        buckets = {col: set() for (_f, col) in per_species_fields[name]}
        for pid in id_list:
            data = master_dict.get(clean_protein_id(pid))
            if not data:
                continue
            for field, col in per_species_fields[name]:
                val = clean_text_field(data.get(field))
                if val:
                    buckets[col].add(val)
        return {col: "; ".join(sorted(vals)) for col, vals in buckets.items()}

    print("Applying exist-to-keep filtering...")
    final_rows = []
    og_count = 0
    out_cols = [f"{target}_ProteinID", "Original_ID", "Orthogroup",
                "Total_Score", "Filter_Reason"] + ref_annotation_columns(sp)

    for _, row in df_og.iterrows():
        og_id = row["Orthogroup"]
        raw_ids = {name: split_cell(row, name) for name in ref_names}
        hits = {name: [p for p in (clean_protein_id(x) for x in raw_ids[name]) if p in immune_set]
                for name in ref_names}
        if not any(hits[name] for name in ref_names):
            continue

        og_count += 1
        reasons = [f"{name}_Hit" for name in ref_names if hits[name]]
        annotation = {}
        for name in ref_names:
            annotation.update(aggregate_species(name, raw_ids[name]))

        tgt_raw = row.get(target)
        if pd.isna(tgt_raw):
            continue
        for tgt_gene in [p.strip() for p in str(tgt_raw).split(",") if p.strip()]:
            clean_tgt = tgt_gene.replace("transcript_", "", 1) if tgt_gene.startswith("transcript_") else tgt_gene
            record = {
                f"{target}_ProteinID": clean_tgt,
                "Original_ID": tgt_gene,
                "Orthogroup": og_id,
                "Total_Score": len(reasons),
                "Filter_Reason": ", ".join(reasons),
            }
            record.update(annotation)
            final_rows.append(record)

    result_df = pd.DataFrame(final_rows).reindex(columns=out_cols)
    result_df.to_csv(output_file, index=False, encoding="utf-8-sig",
                     quoting=csv.QUOTE_NONNUMERIC)
    print("Done.")
    print(f"    immune-related orthogroups: {og_count}")
    print(f"    {target} protein rows:      {len(result_df)}")
    print(f"    saved to: {output_file}")


if __name__ == "__main__":
    main()
