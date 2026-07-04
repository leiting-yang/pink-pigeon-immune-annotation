#!/usr/bin/env python3
"""
filter_kofam.py
===============
Three-step processing of the raw KofamScan detail-tsv output:

  Step 1  parse the detail format (strip the leading '*' marker, split columns)
  Step 2  keep only rows where score > threshold          -> step1_filtered_score.txt
  Step 3  keep only KOs that belong to a KEGG Immune-system pathway
                                                            -> step2_immune_only.txt
  Step 4  deduplicate per transcript, keeping the lowest E-value
                                                            -> step3_final_deduplicated.txt

Immune KOs are collected from the KEGG hierarchy:
  ko00001 -> 09150 Organismal Systems -> 09151 Immune system -> pathways -> KOs
"""

import argparse
import json
import os
import sys

import pandas as pd


def load_config(path):
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--input", help="Raw KofamScan detail-tsv file")
    p.add_argument("--kegg-json", help="KEGG ko00001.json")
    p.add_argument("--outdir", help="Directory for the step outputs")
    return p.parse_args()


def parse_kofam_detail(filepath):
    """Parse KofamScan detail format, handling the leading '*' and ragged spacing."""
    print(f"Step 1: reading {filepath} ...")
    rows = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Drop the leading '*' significance marker so it does not shift columns.
            clean = line.lstrip("*").strip()
            # Columns: gene, KO, thrshld, score, E-value, definition (definition may hold spaces)
            parts = clean.split(maxsplit=5)
            if len(parts) < 5:
                continue
            if len(parts) < 6:
                parts.append("")
            rows.append({
                "gene_name": parts[0],
                "KO": parts[1],
                "thrshld": parts[2],
                "score": parts[3],
                "E_value": parts[4],
                "KO_definition": parts[5],
            })

    df = pd.DataFrame(rows)
    for col in ("thrshld", "score", "E_value"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"  -> parsed rows: {len(df)}")
    return df


def build_immune_ko_map(kegg_data):
    """Return {KO_id: 'pathway1; pathway2'} for KOs under the Immune system node."""
    ko_to_pathway = {}
    top = kegg_data.get("children", [])
    organismal = next((n for n in top if n["name"] == "09150 Organismal Systems"), None)
    if not organismal:
        return ko_to_pathway
    immune = next((n for n in organismal["children"] if n["name"] == "09151 Immune system"), None)
    if not immune:
        return ko_to_pathway
    for pathway in immune["children"]:
        pathway_name = pathway["name"]
        for ko_item in pathway.get("children", []):
            ko_id = ko_item["name"].split()[0]
            ko_to_pathway.setdefault(ko_id, []).append(pathway_name)
    return {ko: "; ".join(paths) for ko, paths in ko_to_pathway.items()}


def main():
    args = parse_args()
    cfg = load_config(args.config)

    kof = cfg.get("kofamscan", {})
    proc_dir = args.outdir or kof.get("processing_dir", ".")
    work_dir = kof.get("work_dir", "")

    input_file = args.input or (os.path.join(work_dir, kof.get("kofam_detail", "result_kofam_detail.txt"))
                                if work_dir else kof.get("kofam_detail", "result_kofam_detail.txt"))
    kegg_json = args.kegg_json or os.path.join(proc_dir, kof.get("kegg_json", "ko00001.json"))

    os.makedirs(proc_dir, exist_ok=True)
    out_score = os.path.join(proc_dir, kof.get("step_filtered_score", "step1_filtered_score.txt"))
    out_immune = os.path.join(proc_dir, kof.get("step_immune_only", "step2_immune_only.txt"))
    out_final = os.path.join(proc_dir, kof.get("step_deduplicated", "step3_final_deduplicated.txt"))

    if not os.path.exists(input_file):
        print(f"ERROR: input not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    df = parse_kofam_detail(input_file)

    # Step 2: score > threshold ('-' thresholds become NaN and are dropped).
    df_score = df[df["score"] > df["thrshld"]].copy()
    df_score.to_csv(out_score, sep="\t", index=False)
    print(f"Step 2: score > threshold -> {len(df_score)} rows ({out_score})")

    # Step 3: keep KOs in the Immune-system hierarchy.
    if not os.path.exists(kegg_json):
        print(f"ERROR: KEGG JSON not found: {kegg_json} (run prepare_db.py first)", file=sys.stderr)
        sys.exit(1)
    print("Step 3: extracting Immune-system KOs from KEGG JSON ...")
    with open(kegg_json) as fh:
        kegg_data = json.load(fh)
    ko_map = build_immune_ko_map(kegg_data)

    df_score["Immune_Pathway"] = df_score["KO"].map(ko_map)
    df_immune = df_score.dropna(subset=["Immune_Pathway"]).copy()
    df_immune.to_csv(out_immune, sep="\t", index=False)
    print(f"Step 3: immune-only -> {len(df_immune)} rows ({out_immune})")

    # Step 4: deduplicate per transcript, keep lowest E-value.
    df_sorted = df_immune.sort_values(by=["gene_name", "E_value"], ascending=[True, True])
    df_final = df_sorted.drop_duplicates(subset=["gene_name"], keep="first")
    df_final.to_csv(out_final, sep="\t", index=False)
    print(f"Step 4: deduplicated -> {len(df_final)} rows ({out_final})")


if __name__ == "__main__":
    main()
