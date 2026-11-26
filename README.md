# Scripts for the shotgun metagenomic analysis of a partial nitritation/anammox structured-bed reactor (SBRIA)

Scripts used to generate the numbers, tables and figures of the metagenomic sections (Methods 2.5, Results and Discussion 3.5, Supplementary Methods S1, Supplementary Results S2, Supplementary Tables and Figures) of the manuscript on nitrogen removal under sensor-guided aeration in a structured-bed reactor. Sample NGS1140_6S (day 595, end of Phase VI).

The repository contains code only. Raw reads and the 25 bins are deposited in the European Nucleotide Archive (ENA) under BioProject accession PRJEB130422. Intermediate files (assembly, BAM files, Kaiju output) are not included because of their size.

## Layout

| Folder | Content |
|---|---|
| `results/review_checks/scripts/` | Read-pair linkage of *hzs*/*hdh* contigs to bins (`link_hzs_contigs_to_bins.py`), attribution, depth and GC of candidate contigs (`attribution_18.py`), taxonomic classification of candidate contigs (`classify_anammox_contigs.py`), comparison of the targeted search with NCycDB (`ncycdb_vs_tableS2.py`), rRNA/tRNA detection and MIMAG categories (`rrna_trna_16.*`), supporting scripts (`base.py`, `hits.py`, `eutils.py`) and the figure scripts of the composite attribution figure (`fig_attribution_v10_27.R`) |
| `results/review_checks/qc_update/scripts/` | CheckM, CheckM2 and GUNC comparison (Table S3) |
| `results/review_checks/supp_tables/`, `results/tables_text/` | Construction of Table S1, Table S2 and figure values |
| `results/review_checks/comammox/` | Contig-level checks of the *amo* genes |
| `results/contig_check/work/` | Coverage of contigs by primary alignments (identity ≥97%) |
| `results/16S_biodiversity/` | 16S rRNA gene biodiversity calculations (R) |
| `verificacao_itens_1_2/` | Recalculation of the numbers reported for bin.009 (skani/FastANI against *Ca.* Brocadia genomes, Kaiju profile, read counts, read-pair linkage control) |
| root `*.R` | Final versions of the figure scripts (R) |

## Software versions

fastp v1.3.6; metaSPAdes v3.15.3; MetaBAT2, MaxBin2 and CONCOCT with DAS Tool v1.1.2; CheckM v1.0.18; CheckM2 v1.1.0; GUNC v1.1.1; GTDB-Tk v2.3.2 (R214); DRAM v0.1.2; Kaiju v1.9.0 (RefSeq genomes without eukaryotes, 2024-08-14); Bowtie2 v2.5.5; samtools v1.24; skani v0.4.0; FastANI v1.34; NCBI Datasets v18.39.0; NCycDB (NCyc_95). Analyses were run on macOS (Apple Silicon, ARM64) and in R and Python.

## Notes on use

- The scripts were written for a fixed directory layout (`~/sbria_anammox_shotgun/...`). Paths must be adjusted to reproduce the analyses.
- Several scripts skip a step when its output file already exists; delete the output file to recompute it.
- Comments inside the scripts are partly in Portuguese.

## Citation

Archived at Zenodo: https://doi.org/10.5281/zenodo.24547209 (this DOI always resolves to the latest version).

## License

MIT (see `LICENSE`).
