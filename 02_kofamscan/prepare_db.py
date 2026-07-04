#!/usr/bin/env python3
"""
prepare_db.py
=============
Download the KEGG KO hierarchy (ko00001) as JSON for offline use by the
KofamScan filtering steps.

NOTE (manual fallback): HPC compute nodes are often offline. If the download
fails, fetch the JSON on a machine with internet access and copy it to the
target path instead.
"""

import argparse
import json
import os
import sys

import requests

KEGG_URL = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir="


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
    p.add_argument("--output", help="Destination path for ko00001.json")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)

    if args.output:
        out_path = args.output
    else:
        proc_dir = cfg.get("kofamscan", {}).get("processing_dir", ".")
        json_name = cfg.get("kofamscan", {}).get("kegg_json", "ko00001.json")
        out_path = os.path.join(proc_dir, json_name)

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    print(f"Downloading KEGG hierarchy from: {KEGG_URL}")
    try:
        response = requests.get(KEGG_URL, timeout=30)
        response.raise_for_status()
        data = response.json()
        with open(out_path, "w") as fh:
            json.dump(data, fh, indent=2)
        print(f"Saved KEGG hierarchy to: {out_path}")
    except Exception as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        print("If this node has no internet, download the JSON elsewhere and "
              f"copy it to: {out_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
