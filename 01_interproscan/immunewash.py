#!/usr/bin/env python3
"""
immunewash.py
=============
Filter the merged InterProScan output down to proteins that carry at least one
immune-related GO term.

Inputs:
  - merged InterProScan TSV (headerless, 15 columns; produced by interproscan_run.sh)
  - immune GO-term whitelist (comma-separated; AmiGO export, see README)

Output:
  - CSV with the rows whose GO annotations hit the whitelist, plus a
    `immune_go_term` column listing the matched terms.

A protein is kept if ANY of its GO annotations contains ANY whitelist term.
"""

import argparse
import os
import sys

import pandas as pd

# Standard InterProScan (-goterms) TSV columns. The trailing empty name matches
# the extra tab some InterProScan versions emit at the end of each line.
COLUMN_NAMES = [
    "ID", "sequence_MD5_digest", "len", "Analysis", "signature_accession",
    "signature_description", "start", "stop", "evalue", "status", "date",
    "interpro_annotations_accession", "interpro_annotations_description",
    "go_annotations", "",
]


def load_config(path):
    """Load a YAML config; return {} if no path given."""
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def resolve(cli_value, cfg, *keys, default=None):
    """CLI value wins; otherwise walk nested config keys; else default."""
    if cli_value is not None:
        return cli_value
    node = cfg
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            return default
        node = node[key]
    return node


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="Path to config.yaml")
    p.add_argument("--input", help="Merged InterProScan TSV")
    p.add_argument("--go-terms", help="Immune GO-term whitelist file")
    p.add_argument("--output", help="Output CSV of immune-matched proteins")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    work = resolve(None, cfg, "interproscan", "work_dir", default="")
    merged_rel = resolve(None, cfg, "interproscan", "merged_tsv",
                         default="ipr_out/interproscan_merged.tsv")
    immune_rel = resolve(None, cfg, "interproscan", "immune_results",
                         default="interproscan_immune_results.csv")
    go_default = resolve(None, cfg, "external_data", "go_terms_immune",
                         default="go_terms_immune_system_process.txt")

    input_tsv = resolve(args.input, {}, default=os.path.join(work, merged_rel) if work else merged_rel)
    go_file = resolve(args.go_terms, {}, default=os.path.join(work, go_default) if work else go_default)
    output_csv = resolve(args.output, {}, default=os.path.join(work, immune_rel) if work else immune_rel)

    print("Filtering InterProScan output for immune-related GO terms...")

    # --- Read merged InterProScan TSV (headerless) ---
    if not os.path.exists(input_tsv):
        print(f"ERROR: input not found: {input_tsv}", file=sys.stderr)
        sys.exit(1)
    print(f"Reading {input_tsv} ...")
    df = pd.read_csv(input_tsv, sep="\t", header=None, names=COLUMN_NAMES,
                     dtype=str, engine="python")
    print(f"Read {len(df)} rows")

    # --- Normalize the GO-annotation column ---
    df["go_annotations"] = df["go_annotations"].fillna("")
    for token in ("(InterPro,PANTHER)", "(InterPro)", "(PANTHER)"):
        df["go_annotations"] = df["go_annotations"].str.replace(token, "", regex=False)
    df["go_annotations"] = df["go_annotations"].str.replace("|", ",", regex=False)

    # --- Load immune GO-term whitelist ---
    if not os.path.exists(go_file):
        print(f"ERROR: GO-term whitelist not found: {go_file}", file=sys.stderr)
        sys.exit(1)
    with open(go_file) as fh:
        immune_terms = [t.strip() for t in fh.read().strip().split(",") if t.strip()]
    print(f"Loaded {len(immune_terms)} immune GO terms")

    # --- Filter: keep rows whose GO annotations hit any whitelist term ---
    df["keep"] = False
    df["immune_go_term"] = ""
    for i, row in df.iterrows():
        go_str = str(row["go_annotations"])
        matched = {term for term in immune_terms if term in go_str}
        if matched:
            df.at[i, "keep"] = True
            df.at[i, "immune_go_term"] = " ".join(sorted(matched))

    immune_results = df[df["keep"]]

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    immune_results.to_csv(output_csv, index=False)

    total = len(df)
    kept = len(immune_results)
    print(f"Saved {kept} immune-matched rows to {output_csv}")
    print("Summary:")
    print(f"   total records:  {total}")
    print(f"   immune records: {kept}")
    if total:
        print(f"   match rate:     {kept / total * 100:.2f}%")


if __name__ == "__main__":
    main()
