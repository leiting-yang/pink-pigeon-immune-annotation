#!/usr/bin/env python3
"""
final_filtering_v2.py
=====================
"Exist-to-keep" filtering: for every orthogroup that contains at least one
reference-species immune protein, keep all Pink Pigeon proteins in that group
and attach the reference annotation.

Inputs : Orthogroups.tsv, master_lookup_table.csv
Output : PinkPigeon_Final_Filtered_List.csv
"""

import argparse
import csv
import os
import re

import pandas as pd

COL_MAP = {"Mouse": "Mouse", "Chicken": "Chicken",
           "ZebraFinch": "ZebraFinch", "PinkPigeon": "PinkPigeon"}

OUTPUT_COLS = [
    "PinkPigeon_ProteinID", "Original_ID", "Orthogroup", "Total_Score", "Filter_Reason",
    "Mouse_GeneSymbols", "Mouse_Category1", "Mouse_Subcategory", "Mouse_Description",
    "Mouse_UniProt_Function", "Mouse_Immune_Source",
    "Chicken_GeneSymbols", "Chicken_Description", "Chicken_GO_Terms",
    "ZebraFinch_GeneSymbols", "ZebraFinch_Description", "ZebraFinch_GO_Terms",
]


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


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


def aggregate_metadata(id_list, lookup_dict, fields):
    aggregated = {field: set() for field in fields}
    for pid in id_list:
        clean_pid = clean_protein_id(pid)
        data = lookup_dict.get(clean_pid)
        if not data:
            continue
        for field in fields:
            val = clean_text_field(data.get(field))
            if val:
                aggregated[field].add(val)
    return {k: "; ".join(sorted(v)) for k, v in aggregated.items()}


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
    orth = cfg.get("orthofinder", {})
    work = orth.get("work_dir", "")

    def wd(name):
        return os.path.join(work, name) if work else name

    og_file = args.orthogroups or wd(orth.get("orthogroups_tsv", "Orthogroups.tsv"))
    master_file = args.master or wd(orth.get("master_lookup", "master_lookup_table.csv"))
    output_file = args.output or wd(orth.get("filtered_list", "PinkPigeon_Final_Filtered_List.csv"))

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

    print("Applying exist-to-keep filtering...")
    final_rows = []
    og_count = 0

    for _, row in df_og.iterrows():
        og_id = row["Orthogroup"]
        mouse_raw = split_cell(row, COL_MAP["Mouse"])
        chk_raw = split_cell(row, COL_MAP["Chicken"])
        zf_raw = split_cell(row, COL_MAP["ZebraFinch"])

        found_mouse = [p for p in (clean_protein_id(x) for x in mouse_raw) if p in immune_set]
        found_chk = [p for p in (clean_protein_id(x) for x in chk_raw) if p in immune_set]
        found_zf = [p for p in (clean_protein_id(x) for x in zf_raw) if p in immune_set]

        if not (found_mouse or found_chk or found_zf):
            continue

        og_count += 1
        reasons = []
        if found_mouse:
            reasons.append("Mouse_Hit")
        if found_chk:
            reasons.append("Chicken_Hit")
        if found_zf:
            reasons.append("ZebraFinch_Hit")
        pseudo_score = len(reasons)

        mouse_info = aggregate_metadata(
            mouse_raw, master_dict,
            ["GeneSymbol", "Description", "Category1", "Subcategory",
             "UniProt_function", "Immune_Source"])
        chk_info = aggregate_metadata(chk_raw, master_dict, ["GeneSymbol", "Description", "GO_Terms"])
        zf_info = aggregate_metadata(zf_raw, master_dict, ["GeneSymbol", "Description", "GO_Terms"])

        pp_raw = row.get(COL_MAP["PinkPigeon"])
        if pd.isna(pp_raw):
            continue
        pp_genes = [p.strip() for p in str(pp_raw).split(",") if p.strip()]

        for pp_gene in pp_genes:
            clean_pp_id = pp_gene.replace("transcript_", "", 1) if pp_gene.startswith("transcript_") else pp_gene
            final_rows.append({
                "PinkPigeon_ProteinID": clean_pp_id,
                "Original_ID": pp_gene,
                "Orthogroup": og_id,
                "Total_Score": pseudo_score,
                "Filter_Reason": ", ".join(reasons),
                "Mouse_GeneSymbols": mouse_info["GeneSymbol"],
                "Mouse_Category1": mouse_info["Category1"],
                "Mouse_Subcategory": mouse_info["Subcategory"],
                "Mouse_Description": mouse_info["Description"],
                "Mouse_UniProt_Function": mouse_info["UniProt_function"],
                "Mouse_Immune_Source": mouse_info["Immune_Source"],
                "Chicken_GeneSymbols": chk_info["GeneSymbol"],
                "Chicken_Description": chk_info["Description"],
                "Chicken_GO_Terms": chk_info["GO_Terms"],
                "ZebraFinch_GeneSymbols": zf_info["GeneSymbol"],
                "ZebraFinch_Description": zf_info["Description"],
                "ZebraFinch_GO_Terms": zf_info["GO_Terms"],
            })

    result_df = pd.DataFrame(final_rows)[OUTPUT_COLS]
    result_df.to_csv(output_file, index=False, encoding="utf-8-sig",
                     quoting=csv.QUOTE_NONNUMERIC)

    print("Done.")
    print(f"    immune-related orthogroups: {og_count}")
    print(f"    Pink Pigeon protein rows:   {len(result_df)}")
    print(f"    saved to: {output_file}")


if __name__ == "__main__":
    main()
