# Reference data — provenance and citations

This folder holds the small reference inputs bundled with the pipeline. The large
reference **proteome FASTAs are NOT stored here** (see "Reference proteomes" below
for how to download them). Point `config/config.yaml` at these files (or copy them
into your workspace) when running the pipeline.

> Note on redistribution: these are processed/derived tables. Before publishing
> the repository, confirm you have the right to redistribute the files derived
> from third-party sources (Nandakumar et al. and AvianImmunome, below). If in
> doubt, replace the file with a download/derivation script and keep only the
> citation.

---

## Files in this folder

| File | What it is | Source / citation |
|------|-----------|-------------------|
| `processed_mouse_biomart.csv` | Mouse immune-gene BioMart table (ProteinID, GeneID, GeneSymbol, Description, GO_Terms) | Nandakumar et al. 2025 (see below) |
| `ImmuneGeneFunction_20240520.csv` | Curated mouse immune-gene functional list (Category1 / Subcategory / UniProt_function) | Nandakumar et al. 2025 (see below) |
| `processed_chicken_combined.csv` | Chicken immune-gene table; `Immune_Source = AvianImmunome` | AvianImmunome database (see below) |
| `processed_zebrafinch_combined.csv` | Zebra finch immune-gene table; `Immune_Source = AvianImmunome` | AvianImmunome database (see below) |
| `go_terms_immune_system_process.txt` | Immune-related GO terms (children of GO:0002376 plus "immune" MF/CC hits), comma-separated | Gene Ontology / AmiGO export (see below) |

---

## Citations

**Mouse BioMart table and curated immune-gene function list**
`processed_mouse_biomart.csv`, `ImmuneGeneFunction_20240520.csv`

> Nandakumar M, Lundberg M, Carlsson F, et al. Positive selection on mammalian
> immune genes — effects of gene function and selective constraint. *Molecular
> Biology and Evolution*. 2025;42(1):msaf016.

**Chicken and zebra finch immune-gene lists**
`processed_chicken_combined.csv`, `processed_zebrafinch_combined.csv`
Derived from the **AvianImmunome** database (see the `Immune_Source` column).
Please cite the AvianImmunome resource per its own citation guidance.
<!-- TODO: add the exact AvianImmunome reference / URL you used. -->

**Immune GO-term list**
`go_terms_immune_system_process.txt`
Exported from the Gene Ontology (AmiGO), immune system process GO:0002376 and its
descendants.

> The Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023.
> *Genetics*. 2023;224(1):iyad031.

---

## Reference proteomes (NOT bundled — download from Ensembl)

The three reference-species protein FASTAs are standard Ensembl proteomes and are
too large to store in git. Download the `pep` (protein) FASTA for each and name it
`<SpeciesName>.fa`/`.faa` to match the species names in `config/config.yaml`.

| Species | Assembly | Ensembl file (pep, all) |
|---------|----------|--------------------------|
| Mouse (`Mouse.faa`) | GRCm39 | `Mus_musculus.GRCm39.pep.all.fa.gz` |
| Chicken (`Chicken.faa`) | bGalGal1.mat.broiler.GRCg7b | `Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz` |
| Zebra finch (`ZebraFinch.faa`) | bTaeGut1_v1.p | `Taeniopygia_guttata.bTaeGut1_v1.p.pep.all.fa.gz` |

Download from the Ensembl FTP (record the release you used, e.g. release-112):

```bash
REL=112   # set to the Ensembl release you use
base="https://ftp.ensembl.org/pub/release-${REL}/fasta"
wget "${base}/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz"
wget "${base}/gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz"
wget "${base}/taeniopygia_guttata/pep/Taeniopygia_guttata.bTaeGut1_v1.p.pep.all.fa.gz"
gunzip *.pep.all.fa.gz
# rename to match config species names, e.g.:
mv Mus_musculus.GRCm39.pep.all.fa Mouse.faa
```

> Cite Ensembl for the proteomes, e.g. Harrison PW, et al. Ensembl 2024.
> *Nucleic Acids Research*. 2024;52(D1):D891–D899 — and record the release number.
