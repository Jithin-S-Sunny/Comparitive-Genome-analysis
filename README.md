# Comparitive-Genome-analysis

## Variant Analysis Pipeline

### Overview

This pipeline performs genome variant analysis using Snakemake. It aligns a query genome to a reference genome, detects small variants and structural variations, and annotates them using SnpEff. The pipeline is automated using Snakemake and requires Conda for dependency management.

### Prerequisites:
Install Conda (Miniconda or Anaconda)
Install Snakemake:
conda install -c bioconda snakemake
Ensure required dependencies are installed using the provided Conda environment file.

### Setup

Clone or download this repository.
Create a Conda environment:
conda env create -f envs/variant_analysis.yaml
Activate the environment:
conda activate variant_analysis

### Running the Pipeline

To execute the Snakemake workflow, run the following command:

snakemake --use-conda -j 4

--use-conda ensures that Snakemake creates and uses the specified environment.
-j 4 specifies the number of parallel jobs (adjust as needed).

### Pipeline Steps

Index Reference Genome: Creates an index file for minimap2.
Align Query Genome: Aligns the query genome to the reference.
Sort & Index BAM: Converts SAM to sorted BAM and indexes it.
Mapping Statistics: Generates read alignment statistics.
Variant Calling: Calls SNPs and small indels using bcftools.
Filter Variants: Removes low-quality variants.
Structural Variants Analysis: Detects large genomic variations using MUMmer.
Identify Structural Variations: Extracts specific SV types.
Genome Coverage Analysis: Computes genome coverage using Bedtools.
Variant Annotation: Annotates filtered variants using SnpEff.
Generate Annotation Summary: Creates an annotation report.
Plot Alignment: Visualizes genome alignment using MUMmer.

### Outputs

filtered_variants.vcf: Filtered high-confidence variant calls.
annotated_snpeff.vcf: Annotated variants.
alignment_summary.txt: Structural variation summary.
mummer_plot.png: Genome alignment visualization.
uncovered_regions.bed: Regions with no coverage.

### Troubleshooting

Ensure all dependencies are installed by running:
conda list
If a job fails, rerun with:
snakemake --use-conda --rerun-incomplete
Check logs for error details in .snakemake/log/.
