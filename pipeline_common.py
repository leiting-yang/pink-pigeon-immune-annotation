#!/usr/bin/env python3
"""
pipeline_common.py
==================
Shared helpers so every stage resolves species names and column names from the
same single source of truth (`species:` in config.yaml). This is what makes the
pipeline species-agnostic: change the config, and all reference-species columns
and logic follow.

Scripts in the stage sub-folders import this by adding the repo root to sys.path:

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from pipeline_common import load_config, get_species
"""


def load_config(path):
    """Load a YAML config; return {} if no path is given."""
    if not path:
        return {}
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def get_species(cfg):
    """Resolve the species configuration into a plain dict.

    Returns keys:
      target               focal species name (e.g. 'PinkPigeon')
      ref_names            list of reference species names, in config order
      ref_map              {name: {'name':.., 'biomart':.., 'curated_list':..}}
      curated              reference names that carry a curated functional list
      consensus_group      species whose intersection is high-confidence ([] to disable)
      prediction_priority  single-species fallback order for symbol prediction
      gene_info_enrichment bool, whether to run the gene_info.csv step
    """
    sp = cfg.get("species", {}) or {}
    target = sp.get("target", "PinkPigeon")

    ref_entries = []
    for r in sp.get("reference", []) or []:
        ref_entries.append(r if isinstance(r, dict) else {"name": r})

    ref_names = [r["name"] for r in ref_entries]
    ref_map = {r["name"]: r for r in ref_entries}
    curated = [r["name"] for r in ref_entries if r.get("curated_list")]

    return {
        "target": target,
        "ref_names": ref_names,
        "ref_map": ref_map,
        "curated": curated,
        "consensus_group": sp.get("consensus_group", []) or [],
        "prediction_priority": sp.get("prediction_priority", ref_names) or ref_names,
        "gene_info_enrichment": sp.get("gene_info_enrichment", True),
    }


# master_lookup field -> output-column suffix. Every reference species gets the
# base fields; species that carry a curated list additionally get the curated
# fields. Producer (final_filtering) and consumers (merge_immune_annotations,
# process_symbols, ...) all build their column names from here so they agree.
BASE_FIELD_SUFFIX = [
    ("GeneSymbol", "GeneSymbols"),
    ("Description", "Description"),
    ("GO_Terms", "GO_Terms"),
]
CURATED_FIELD_SUFFIX = [
    ("Category1", "Category1"),
    ("Subcategory", "Subcategory"),
    ("UniProt_function", "UniProt_Function"),
    ("Immune_Source", "Immune_Source"),
]


def species_field_map(sp):
    """List of (species, master_field, output_column) for all reference species."""
    out = []
    for name in sp["ref_names"]:
        for field, suffix in BASE_FIELD_SUFFIX:
            out.append((name, field, f"{name}_{suffix}"))
        if name in sp["curated"]:
            for field, suffix in CURATED_FIELD_SUFFIX:
                out.append((name, field, f"{name}_{suffix}"))
    return out


def ref_annotation_columns(sp):
    """Ordered per-species annotation column names (no target/meta columns)."""
    return [col for (_sp, _field, col) in species_field_map(sp)]


def symbol_columns(sp):
    """Per-reference-species gene-symbol columns, in config order."""
    return [f"{name}_GeneSymbols" for name in sp["ref_names"]]
