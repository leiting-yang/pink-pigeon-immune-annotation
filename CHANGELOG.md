# Changelog

## Unreleased — portability fixes and cleanup after a second-species test run

Validated by running the full pipeline end to end on a new target species
(Mauritius Kestrel, *Falco punctatus*, NCBI RefSeq GCF_963210335.1). The run
surfaced several bugs and rough edges that were not specific to that species;
all are fixed below.

### Bug fixes (would affect any dataset)

1. **`01_interproscan/immunewash.py` — path resolution returned a dict.**
   `resolve()` was called with no config keys, so it returned the empty `{}`
   instead of the computed default, and `os.path.exists()` raised
   `TypeError: stat: path should be string ... not dict`. The three input/output
   paths are now resolved directly (CLI value else computed default).

2. **`01_interproscan/immunewash.py` — InterProScan TSV parsing crashed on `"`.**
   Description fields can contain a literal double quote, which the CSV reader
   treated as a quote character (`ParserError: '\t' expected after '"'`). The
   merged TSV is now read with `quoting=csv.QUOTE_NONE`.

3. **`04_tiering/add_gene_info_3.py` — empty native-symbol file crashed.**
   An empty `native_symbols` file (explicitly allowed) raised
   `EmptyDataError: No columns to parse`. `load_native_symbols()` now catches it
   and returns an empty map.

4. **`04_tiering/merge_og_stats_to_final_gene_table.py` — pandas incompatibility.**
   Newer pandas excludes the grouping column from each group in
   `groupby(...).apply()`, so `primary["GeneID"]` raised `KeyError: 'GeneID'`.
   The gene id is now read from the group key (`group.name`).

### Robustness / usability improvements

5. **Stage 00 cleans gffread output.** `00_preprocess/extract_proteins.sh` now
   replaces `.` (untranslatable codon) with `X` and removes `*` (stop) from the
   protein sequences. InterProScan, diamond (OrthoFinder) and HMMER (KofamScan)
   all reject those characters, so cleaning once here prevents mid-run failures
   in all three downstream tools.

6. **Deterministic OrthoFinder output folder.** `03_orthofinder/run_orthofind.sh`
   now runs OrthoFinder with `-n "$ORTHO_RUN_NAME"`, producing a fixed
   `Results_<name>` folder instead of the date-stamped `Results_<date>`, so the
   OrthoFinder paths in `config.yaml` no longer need editing between runs.

7. **KofamScan is now in the conda env.** `envs/environment.yml` adds
   `kofamscan`; the README documents downloading its HMM database (profiles +
   `ko_list`) and pointing `KOFAM_DB_CONFIG` at a `config.yml`.

8. **Tidy per-stage outputs.** Config inputs use absolute paths; each stage
   writes to `<workspace>/outputs/0X_<stage>/`. `config/cluster.sh` gains
   `WORKSPACE`, `PROTEIN_FASTA_PATH`, `IPR_WORK_DIR`, `KOFAM_WORK_DIR` and
   `ORTHO_RUN_NAME`; `interproscan_run.sh` and `run_kofam.sh` read those (with
   fallbacks to the old behaviour).

9. **Native symbols from the GFF (optional, ON by default).**
   `00_preprocess/extract_native_symbols.sh` builds the `Symbol <TAB> GeneID`
   table from the annotation GFF3 (`gene=`/`Name=` on gene features). Skip it if
   you already have a native symbol table.

10. **Bundled reference data + citations.** `reference_data/` holds the small
    BioMart/immune-gene tables and the GO-term list, with
    `reference_data/CITATIONS.md` documenting provenance and how to download the
    large reference proteome FASTAs (kept out of git).

### Notes for NCBI RefSeq inputs

The transcript/gene IDs in a RefSeq annotation carry `rna-` / `gene-` prefixes
and a `.N` version. The pipeline handles these unchanged: the GFF `tx->gene` map
and the InterProScan/KofamScan/OrthoFinder IDs all keep the same prefixes and
version consistently, so the join keys line up without editing the ID-cleaning
`.replace(...)` calls. Verify on a first run by checking the "failed to map" /
"mapped X/Y" counts printed by the stage-02 and stage-04 merges are near zero.
