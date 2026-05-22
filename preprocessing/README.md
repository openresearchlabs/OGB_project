# Metagenomic preprocessing

Raw paired-end metagenomic reads were processed on [CSC Puhti](https://www.csc.fi/en/puhti) before taxonomic and functional profiling. Scripts are in [`raw_reads/`](raw_reads/).

## Raw read processing (taxprofiler)

Pipeline: [nf-core/taxprofiler](https://nf-co.re/taxprofiler) v1.1.5 (Nextflow).

Steps enabled in `raw_reads/taxprofiler_runs.sh`:

1. **QC** — FastQC
2. **Trimming** — fastp (adapters and low-quality bases)
3. **Complexity filter** — BBduk
4. **Host removal** — Bowtie2 against human reference T2T-CHM13v2
5. **Taxonomy** — MetaPhlAn4 (`merge_metaphlan_tables.py` for merged profiles)

Run example: `sbatch raw_reads/taxprofiler_runs.sh` (edit samplesheet and config paths for your Puhti project).

More setup notes: [workflow_metagenome](https://github.com/erawijantari/workflow_metagenome).

## Functional profiling (optional)

HUMAnN3 via Snakemake (`raw_reads/Snakefile_humann.txt`, submit with `raw_reads/run_humann_workflow.sh`) on host-depleted FASTQs (`*.unmapped_1.fastq.gz`). Post-processing: `raw_reads/humann_modif.sh`.

## Downstream (this repo)

MetaPhlAn abundance tables used in the Quarto pipeline are cleaned with:

- `remove_columns_metaphlan_db_meta4_combined_reports.R` — drop sample columns and write `latest_metaphlan_db_meta4_combined_reports.txt`

## Note on PCR deduplication

Read deduplication was **not** included in this workflow.
