# Nanopore tumor-vs-normal pipeline

Snakemake workflow for long-read Oxford Nanopore tumor analysis with a matched normal reference BAM.

## Current workflow

This pipeline currently performs:

1. Merge/pass all FASTQ files from two input datasets
2. Align reads with `minimap2`
3. Filter BAM with `samtools` (`MAPQ >= 20`, excluding supplementary/secondary flags)
4. Generate QC report with `NanoPlot`
5. Call structural variants with `cuteSV`
6. Call structural variants with `Sniffles`
7. Generate `.seqz` input for copy-number / HRD analysis using a matched normal BAM
8. Run HRD scoring and alternative copy-number solution plotting scripts

## Files

- `Snakefile`
- `config/config.example.yaml`
- `.gitignore`

## Expected config

```yaml
dataset1: /path/to/fastq_pass_run1
dataset2: /path/to/fastq_pass_run2
output_dir: /path/to/output
barcode: SAMPLE_BARCODE
bam_normal: /path/to/normal_sample.bam
```

## Example run

```bash
snakemake --cores 32 --use-conda --configfile config/config.example.yaml
```

## Important notes before publishing

This workflow currently contains environment-specific absolute paths in the original project, including:

- reference genome path
- GC wiggle file path
- R and Python helper scripts for HRD plotting
- local conda environment names

Before publishing, replace hard-coded local paths with configuration variables or document them clearly.

## Suggested improvements for the public repository

- Move reference genome and helper script paths into `config.yaml`
- Add `envs/*.yaml` conda environment files
- Add a workflow diagram
- Add test data or a minimal dry-run example
- Add license and citation information
