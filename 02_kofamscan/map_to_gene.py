#!/usr/bin/env python3
"""
map_to_gene.py
==============
Lift the transcript-level KofamScan results to the gene level using the
transcript-to-gene map, then deduplicate again per gene (keep the lowest
E-value).

Inputs : step3_final_deduplicated.txt, transcript_to_gene.txt
Output : final_kofam_gene_level.tsv
"""

import argparse
import os
import sys

import pandas as pd

DESIRED_COLS = ["GeneID", "KO", "Immune_Pathway", "score", "E_value", "KO_definition", "gene_name"]


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
    p.add_argument("--input", help="Transcript-level deduplicated KofamScan TSV")
    p.add_argument("--mapping", help="transcript_to_gene TSV")
    p.add_argument("--output", help="Gene-level output TSV")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    kof = cfg.get("kofamscan", {})
    proc_dir = kof.get("processing_dir", ".")

    input_file = args.input or os.path.join(proc_dir, kof.get("step_deduplicated", "step3_final_deduplicated.txt"))
    mapping_file = args.mapping or os.path.join(proc_dir, kof.get("transcript_to_gene", "transcript_to_gene.txt"))
    output_file = args.output or os.path.join(proc_dir, kof.get("gene_level", "final_kofam_gene_level.tsv"))

    for path in (input_file, mapping_file):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    print("Mapping transcript-level KO annotations to gene level...")
    df_kofam = pd.read_csv(input_file, sep="\t")
    print(f"  -> KofamScan rows: {len(df_kofam)}")

    # Strip 'transcript:' so IDs match the mapping table keys.
    df_kofam["clean_tid"] = df_kofam["gene_name"].astype(str).str.replace("transcript:", "", regex=False)

    df_map = pd.read_csv(mapping_file, sep="\t")
    df_map.columns = ["clean_tid", "GeneID"]
    df_map["clean_tid"] = df_map["clean_tid"].astype(str).str.strip()
    df_map["GeneID"] = df_map["GeneID"].astype(str).str.strip()
    print(f"  -> mapping rows: {len(df_map)}")

    df_merged = pd.merge(df_kofam, df_map, on="clean_tid", how="left")

    missing = df_merged["GeneID"].isna().sum()
    if missing:
        print(f"  -> {missing} entries had no GeneID; filled with transcript ID.")
        df_merged["GeneID"] = df_merged["GeneID"].fillna(df_merged["clean_tid"])
    else:
        print("  -> all entries mapped to a GeneID.")

    # Gene-level deduplication: keep the lowest E-value per gene.
    df_merged["E_value"] = pd.to_numeric(df_merged["E_value"], errors="coerce")
    df_merged = df_merged.sort_values(by=["GeneID", "E_value"], ascending=[True, True])
    before = len(df_merged)
    df_final = df_merged.drop_duplicates(subset=["GeneID"], keep="first")
    print(f"  -> gene-level dedup: {before} -> {len(df_final)} rows")

    final_cols = [c for c in DESIRED_COLS if c in df_final.columns]
    df_final = df_final[final_cols]
    df_final.to_csv(output_file, sep="\t", index=False)
    print(f"Saved gene-level KofamScan table: {output_file}")


if __name__ == "__main__":
    main()
