# finch-go-reference-sets

Zebra finch (*Taeniopygia guttata*) gene annotation, strict 1:1 orthologs (human, mouse, chicken) and Gene Ontology gene sets for enrichment analysis with **ShinyGO** and **mitch**.

The build tool, `scripts/build_reference_sets.py`, combines public Ensembl, NCBI RefSeq and Gene Ontology files into ready-to-use tables. It covers both zebra finch genome versions:

| Genome | Accession | Annotation | Gene IDs look like |
|---|---|---|---|
| **bTaeGut7.mat** (new, maternal haplotype, Jarvis/Formenti 2026; autosomes + Z + W) | GenBank GCA_048771995.1 / RefSeq GCF_048771995.1 | Ensembl `ensemblorganisms` 2026_05; NCBI RefSeq RS_2025_04 | `ENSTGUG00025xxxxxx` |
| **bTaeGut1_v1** (old) | GCA_003957565.2 | Ensembl release 104 | `ENSTGUG00000xxxxxx` |

Use the bTaeGut1_v1 files for data counted against Ensembl 104. Use the bTaeGut7 files for data aligned to the new assembly. `crosswalk_bTaeGut1v1_to_bTaeGut7.tsv` carries gene IDs from old to new.

Developed in the Julie E. Miller Lab, University of Arizona, by Aleco Dominguez.

## Quick start: use the prebuilt files

The `reference_sets/` folder holds the outputs of the last build (October 2026). The GMT files are gzipped to keep them under GitHub's file-size limits. Unzip them before use:

```bash
gunzip reference_sets/*.gmt.gz
```

### ShinyGO (https://bioinformatics.sdstate.edu/go/)

ShinyGO takes a pasted gene list and a pasted background. It does not accept custom GMT files.

1. Convert your finch gene IDs to human symbols with the annotation table for your genome version (`human_symbol` column for bTaeGut1_v1, `human_ortholog_symbol_consensus` for bTaeGut7).
2. In ShinyGO, set species to **Human** and paste the symbols.
3. Click **Background (recommended)** and paste your *expressed* genes, converted to human symbols the same way. If you have no expressed-gene list, `ShinyGO_background_*_all_human_symbols.txt` is a genome-wide fallback, but it inflates significance. Never leave the background empty.

### mitch (R)

```r
genesets <- gmt_import("GO_finch_bTaeGut1v1_ENS104_EnsemblID_combined.gmt")
m <- mitch_import(x, DEtype = "deseq2")   # profile already keyed by old finch Ensembl IDs
res <- mitch_calc(m, genesets, priority = "effect")
```

To use human-symbol sets instead, pass `geneTable = read.table("mitch_geneTable_bTaeGut1v1_to_human_symbol.tsv")` together with `GO_human_symbol_goa.gmt`.

## Files in `reference_sets/`

| File | What it is | Use it for |
|---|---|---|
| `finch_bTaeGut7_gene_annotation.tsv` | One row per bTaeGut7 Ensembl gene (28,197). Symbol, biotype, description, coordinates; NCBI GeneID/symbol/description matched by coordinate overlap on the same assembly; 1:1 human, mouse and chicken orthologs; consensus human ortholog symbol; `display_symbol`; `ensembl_ncbi_symbol_conflict` flag. | Master annotation for bTaeGut7 data |
| `finch_bTaeGut1v1_ENS104_gene_annotation.tsv` | One row per Ensembl 104 gene (22,150). Strict Compara `ortholog_one2one` human (Compara 104) and mouse (Compara 116) orthologs, plus each gene's crosswalk to bTaeGut7. | Annotating Ensembl 104 data |
| `crosswalk_bTaeGut1v1_to_bTaeGut7.tsv` | Old → new gene ID map with evidence and confidence. | Moving results to the new genome |
| `orthologs_finch_bTaeGut1v1_to_mouse_1to1_Compara116_highconf.tsv` (10,120 pairs), `..._Compara116_all.tsv` (11,871) | Old-ID finch–mouse 1:1 orthologs from Ensembl Compara 116, with mouse symbols and NCBI GeneIDs. | Finch–mouse comparisons (high-confidence first) |
| `orthologs_finch_bTaeGut7_to_mouse_1to1.tsv` | New-ID finch–mouse 1:1 orthologs (BLAST reciprocal best hits). | Finch–mouse comparisons on bTaeGut7 |
| `GO_finch_bTaeGut1v1_ENS104_EnsemblID_combined.gmt.gz` | GO sets keyed by old Ensembl IDs. | mitch on Ensembl 104 data |
| `GO_finch_bTaeGut7_EnsemblID_combined.gmt.gz` | GO sets keyed by new Ensembl IDs: finch-native GO plus human GO projected through 1:1 orthologs. | mitch on bTaeGut7 data |
| `GO_finch_bTaeGut7_EnsemblID_finchNativeOnly.gmt.gz` | Finch-native GO only (Ensembl xref + NCBI GAF). | Sensitivity check |
| `GO_finch_bTaeGut7_NCBIGeneID_RefSeqGAF.gmt.gz` | NCBI RefSeq GO, keyed by NCBI GeneID. | Data counted against RefSeq |
| `GO_human_symbol_goa.gmt.gz` | Human GOA, keyed by HGNC symbol. | mitch after converting to human symbols |
| `mitch_geneTable_*_to_human_symbol.tsv` | Two columns, no header: finch ID → human symbol. | `geneTable=` in `mitch_import()` |
| `ShinyGO_background_*_all_human_symbols.txt` | All finch genes with a 1:1 human ortholog, as human symbols. | Genome-wide ShinyGO background (fallback only) |
| `build_log.txt` | Counts from the build. | Checking a rebuild |

GO sets are propagated up the GO hierarchy through `is_a` and `part_of`. Each set has 5–1000 genes. Annotations with the NOT qualifier or ND evidence are dropped. Set names follow `GO:0005739 mitochondrion [CC]`.

## How the 1:1 ("strict") orthologs are defined

- **bTaeGut1_v1 IDs:** Ensembl Compara gene-tree orthologs restricted to `ortholog_one2one`. Human from Compara 104 (11,631; 11,033 high-confidence). Mouse from Compara 116 (11,871; 10,120 high-confidence) when that file is present, otherwise Compara 104. Filter on `*_high_confidence` for the strictest set.
- **bTaeGut7 IDs:** the 2026_05 `ensemblorganisms` release has no Compara gene trees yet, only BLAST reciprocal best hits (`homolog_rbbh`). Non-1:1 RBBHs are dropped (9,956 human, 9,809 mouse, 10,470 chicken). NCBI `gene_orthologs` calls are added, and a **consensus** human symbol is recorded (13,330 genes). Symbols used by more than one finch gene are removed so the consensus stays 1:1.

## How the old → new crosswalk is built

Evidence per old gene: identical protein sequence; DIAMOND reciprocal best protein hit (≥80% identity, ≥50% coverage); same gene symbol (only when unique in both annotations); same 1:1 human ortholog.

- **High:** ≥2 lines of evidence, or an identical protein (11,200 genes)
- **Medium:** one line of evidence (2,392)
- **Ambiguous:** a tie between candidates (93)
- **Unmapped:** 2,934 protein-coding and 5,531 non-coding genes (non-coding genes would need coordinate liftover)

## Rebuild from source

Requirements: Python 3.10, pandas 2.3.3. About 20 seconds.

```bash
pip install pandas==2.3.3
python3 scripts/build_reference_sets.py      # run from the repository root
```

The script reads `reference_sources/` and writes `reference_sets/` (uncompressed) plus `build_log.txt`. Raw downloads are not stored in this repository (`.gitignore`). Download them and save them under the exact file names below.

| File name in `reference_sources/` | Source | Required? |
|---|---|---|
| `ENS2026_bTaeGut7mat_genes.gff3.gz`, `ENS2026_bTaeGut7mat_xref.tsv.gz`, `ENS2026_bTaeGut7mat_pep.fa.gz`, `ENS2026_bTaeGut7mat_homology.tsv.gz` | Ensembl `ensemblorganisms` (ftp.ebi.ac.uk), GCA_048771995.1, 2026_05 gene set | Yes |
| `ENS104_bTaeGut1v1_genes.gff3.gz`, `ENS104_bTaeGut1v1_xref.tsv.gz`, `ENS104_bTaeGut1v1_pep.fa.gz` | Ensembl `ensemblorganisms`, GCA_003957565.2 (2019_12 = release 104) | Yes |
| `ENS104_Compara_homo_sapiens_vs_tguttata.tsv` | Ensembl Compara 104 protein homologies (ftp.ensembl.org), human file filtered to rows with zebra finch | Optional (no old-ID human orthologs without it) |
| `ENS116_Compara_mus_musculus_vs_tguttata.tsv` (or `ENS104_Compara_mus_musculus_vs_tguttata.tsv`) | Ensembl Compara 116 (or 104) mouse homologies, filtered to rows with zebra finch | Optional |
| `NCBI_GCF_048771995.1_bTaeGut7.mat_genomic.gff.gz`, `NCBI_GCF_048771995.1_bTaeGut7.mat_assembly_report.txt`, `NCBI_GCF_048771995.1-RS_2025_04_gene_ontology.gaf.gz` | NCBI RefSeq GCF_048771995.1 | Yes |
| `NCBI_gene_orthologs_taxid59729_subset.tsv` | NCBI `gene_orthologs`, filtered to taxon 59729 | Yes |
| `NCBI_Homo_sapiens.gene_info.gz` | NCBI `Homo_sapiens.gene_info` | Optional |
| `NCBI_mouse_gene_annotations_merged_FMcCarthy.tsv` | NCBI GRCm39 mouse genes with Ensembl IDs | Optional (improves mouse symbols) |
| `GO_go-basic.obo` | geneontology.org | Yes |
| `GO_goa_human.gaf.gz` | EBI GOA `goa_human.gaf` | Optional (no human-projected sets without it) |
| `crosswalk_protein_besthits_diamond.tsv` | Included here. DIAMOND v2.1.10 blastp, both directions between the two proteomes, `--fast`, `--evalue 1e-10` | Optional |

Columns the script expects in the non-standard files:

- Compara files: `gene_stable_id`, `homology_gene_stable_id`, `homology_species`, `homology_type`, `homology_identity`, `is_high_confidence`
- `ENS2026_bTaeGut7mat_homology.tsv.gz`: `ref_species`, `homology_type`, `query_gene_stable_id`, `ref_gene_stable_id`, `ref_gene_name`, `query_perc_id`
- NCBI orthologs: `#tax_id`, `GeneID`, `Other_tax_id`, `Other_GeneID`
- Mouse annotation: `Ensembl_GeneID`, `Symbol`, `GeneID`
- DIAMOND table: `old_gene`, `new_gene`, `pident`, `qcov`, `scov`, `bits`, `reciprocal`

**Check a rebuild.** With the October 2026 sources, `build_log.txt` must match the copy in this repository. Two spot checks: PAX6 (ENSTGUG00000004824) maps to ENSTGUG00025009380 with all evidence types, and SNCA is ENSTGUG00025000678 = NCBI GeneID 751582.

## Notes

- **Ensembl and NCBI disagree on 662 gene symbols.** These are flagged in `ensembl_ncbi_symbol_conflict`. Where one of the two symbols matches the 1:1 human ortholog, `display_symbol` uses it. Check this flag before naming a gene in a figure.
- **One-to-many genes.** Some genes have two old finch paralogs that both map to one new gene. Because they are not 1:1, they have no 1:1 human ortholog.
- **Compara 116 vs. 104.** Compara 116 still uses bTaeGut1_v1 IDs. Of the one-to-one finch–mouse pairs, 11,540 are identical between releases; 331 are new in 116 and 169 are gone. The bTaeGut7 BLAST orthologs agree with Compara 116 for 99.7% of the testable high-confidence pairs.
- BioMart is being retired, so orthology here is rebuilt from the Ensembl FTP homology files.

## Citing

If you use these files, please cite the underlying resources (Ensembl, NCBI RefSeq, Gene Ontology / GOA) and, for ShinyGO, Ge SX, Jung D, Yao R (2020) *Bioinformatics* 36(8):2628–2629, doi:10.1093/bioinformatics/btz931.
