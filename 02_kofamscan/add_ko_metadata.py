#!/usr/bin/env python3
"""
add_ko_metadata.py
==================
Add per-KO gene symbol and gene name (parsed from the KEGG hierarchy) to the
gene-level KofamScan table.

Inputs : final_kofam_gene_level.tsv, ko00001.json
Output : final_kofam_annotated.tsv  (adds KO_Gene_Symbol, KO_Gene_Name)
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
    p.add_argument("--input", help="Gene-level KofamScan TSV")
    p.add_argument("--kegg-json", help="KEGG ko00001.json")
    p.add_argument("--output", help="Annotated output TSV")
    return p.parse_args()


def parse_kegg_symbols(json_path):
    """Walk the KEGG tree; return {KO: {'symbol': ..., 'name': ...}}.

    Leaf names look like: 'K00001  E1.1.1.1, adh; alcohol dehydrogenase [EC:...]'
    Symbol is the part before ';', name is the part after.
    """
    print(f"Parsing KEGG hierarchy: {json_path} ...")
    with open(json_path) as fh:
        data = json.load(fh)

    ko_info = {}

    def traverse(node):
        if "children" not in node:
            name_str = node.get("name", "")
            if not name_str.startswith("K"):
                return
            parts = name_str.split(maxsplit=1)
            if len(parts) < 2:
                return
            ko_id, desc = parts[0], parts[1]
            if ";" in desc:
                symbol, name = desc.split(";", 1)
                ko_info[ko_id] = {"symbol": symbol.strip(), "name": name.strip()}
            else:
                ko_info[ko_id] = {"symbol": "-", "name": desc.strip()}
        else:
            for child in node["children"]:
                traverse(child)

    traverse(data)
    print(f"  -> parsed {len(ko_info)} KO annotations")
    return ko_info


def main():
    args = parse_args()
    cfg = load_config(args.config)
    kof = cfg.get("kofamscan", {})
    proc_dir = kof.get("processing_dir", ".")

    input_tsv = args.input or os.path.join(proc_dir, kof.get("gene_level", "final_kofam_gene_level.tsv"))
    kegg_json = args.kegg_json or os.path.join(proc_dir, kof.get("kegg_json", "ko00001.json"))
    output_tsv = args.output or os.path.join(proc_dir, kof.get("annotated", "final_kofam_annotated.tsv"))

    for path in (input_tsv, kegg_json):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    df = pd.read_csv(input_tsv, sep="\t")
    print(f"Loaded gene-level table: {len(df)} rows")

    ko_dict = parse_kegg_symbols(kegg_json)
    df["KO_Gene_Symbol"] = df["KO"].map(lambda ko: ko_dict.get(ko, {}).get("symbol", "-"))
    df["KO_Gene_Name"] = df["KO"].map(lambda ko: ko_dict.get(ko, {}).get("name", "-"))

    # Place the two new columns right after GeneID.
    cols = [c for c in df.columns if c not in ("KO_Gene_Symbol", "KO_Gene_Name")]
    idx = cols.index("GeneID") + 1 if "GeneID" in cols else 0
    cols[idx:idx] = ["KO_Gene_Symbol", "KO_Gene_Name"]
    df = df[cols]

    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Saved annotated KofamScan table: {output_tsv}")


if __name__ == "__main__":
    main()
