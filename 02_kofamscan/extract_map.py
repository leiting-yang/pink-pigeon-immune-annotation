#!/usr/bin/env python3
"""
extract_map.py
==============
Extract a transcript-to-gene mapping from the GFF3 (from mRNA / transcript
feature lines), used later to lift transcript-level KO annotations to the gene
level.

Output: 2-column TSV (TranscriptID <TAB> GeneID). IDs are stripped of the
common `transcript:` and `gene:` prefixes.
"""

import argparse
import os

import pandas as pd  # noqa: F401  (kept for a consistent env; not strictly required)


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
    p.add_argument("--gff", help="Input GFF3 annotation")
    p.add_argument("--output", help="Output transcript-to-gene TSV")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    gff_file = args.gff or cfg.get("reference", {}).get("gff3_raw")
    if args.output:
        out_path = args.output
    else:
        proc_dir = cfg.get("kofamscan", {}).get("processing_dir", ".")
        name = cfg.get("kofamscan", {}).get("transcript_to_gene", "transcript_to_gene.txt")
        out_path = os.path.join(proc_dir, name)

    if not gff_file:
        raise SystemExit("ERROR: no GFF3 given (use --gff or config reference.gff3_raw)")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    print(f"Extracting transcript -> gene mapping from {gff_file} ...")
    n = 0
    with open(gff_file) as f_in, open(out_path, "w") as f_out:
        f_out.write("TranscriptID\tGeneID\n")
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in ("mRNA", "transcript"):
                continue
            attrs = dict(item.split("=", 1) for item in parts[8].split(";") if "=" in item)
            t_id = attrs.get("ID", "").replace("transcript:", "")
            g_id = attrs.get("Parent", "").replace("gene:", "")
            if t_id and g_id:
                f_out.write(f"{t_id}\t{g_id}\n")
                n += 1

    print(f"Wrote {n} transcript -> gene rows to: {out_path}")


if __name__ == "__main__":
    main()
