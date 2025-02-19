conda: "envs/variant_analysis.yaml"

REFERENCE = "pseudo_ref.fna"
QUERY = "F27_B1.fna"
PREFIX = "genome_comparison"

rule all:
    input:
        f"{PREFIX}.sorted.bam",
        f"{PREFIX}.sorted.bam.bai",
        "filtered_variants.vcf",
        "annotated_snpeff.vcf",
        "alignment_summary.txt",
        "inversions.txt",
        "deletions.txt",
        "low_identity_regions.txt",
        "uncovered_regions.bed",
        "mummer_plot.png"

# Step 1: Index reference genome
rule index_reference:
    input: REFERENCE
    output: REFERENCE + ".mmi"
    shell: "minimap2 -d {output} {input}"

# Step 2: Align query genome to reference
rule align_genome:
    input:
        ref_index = REFERENCE + ".mmi",
        query = QUERY
    output: "strain_aligned.bam"
    shell: "minimap2 -ax asm5 {input.ref_index} {input.query} | samtools view -Sb -o {output}"

# Step 3: Sort and index BAM file
rule sort_index_bam:
    input: "strain_aligned.bam"
    output:
        sorted_bam = f"{PREFIX}.sorted.bam",
        index = f"{PREFIX}.sorted.bam.bai"
    shell:
        "samtools sort -o {output.sorted_bam} {input} && "
        "samtools index {output.sorted_bam}"

# Step 4: Generate mapping statistics
rule mapping_stats:
    input: f"{PREFIX}.sorted.bam"
    output:
        flagstat = "flagstat.txt",
        depth = "depth.txt",
        idxstats = "idxstats.txt"
    shell:
        "samtools flagstat {input} > {output.flagstat} && "
        "samtools depth {input} | awk '{sum+=$3} END {print \"Average depth:\", sum/NR}' > {output.depth} && "
        "samtools idxstats {input} > {output.idxstats}"

# Step 5: Variant calling
rule call_variants:
    input:
        ref = REFERENCE,
        bam = f"{PREFIX}.sorted.bam"
    output: "raw_variants.vcf"
    shell: "bcftools mpileup -Ou -f {input.ref} {input.bam} | bcftools call -mv -Ov -o {output}"

# Step 6: Filter low-quality variants
rule filter_variants:
    input: "raw_variants.vcf"
    output: "filtered_variants.vcf"
    shell: "bcftools filter -i 'QUAL>20 && DP>10' {input} -o {output}"

# Step 7: Structural Variant Analysis with MUMmer
rule structural_variants:
    input:
        ref = REFERENCE,
        query = QUERY
    output:
        delta = f"{PREFIX}.delta",
        snps = "genome_variants.txt"
    shell:
        "nucmer --maxmatch -c 100 {input.ref} {input.query} -p {PREFIX} && "
        "show-snps -Clr {output.delta} > {output.snps}"

# Step 8: Identify Structural Variations
rule identify_sv:
    input: f"{PREFIX}.delta"
    output:
        alignment_summary = "alignment_summary.txt",
        inversions = "inversions.txt",
        deletions = "deletions.txt",
        low_identity = "low_identity_regions.txt"
    shell:
        "show-coords -rcl {input} > {output.alignment_summary} && "
        "grep -E '\\s-\\s' {output.alignment_summary} > {output.inversions} && "
        "awk '$3-$4 > 5000' {output.alignment_summary} > {output.deletions} && "
        "awk '$6 < 95' {output.alignment_summary} > {output.low_identity}"

# Step 9: Genome Coverage Analysis
rule genome_coverage:
    input: f"{PREFIX}.sorted.bam"
    output: "uncovered_regions.bed"
    shell: "bedtools genomecov -ibam {input} -bga > {output}"

# Step 10: Variant Annotation with SnpEff
rule annotate_variants:
    input: "filtered_variants.vcf"
    output: "annotated_snpeff.vcf"
    shell:
        "java -Xmx4G -jar snpEff.jar Pseudomonas_aeruginosa {input} > {output}"

# Step 11: Generate Annotation Summary
rule annotation_summary:
    input: "annotated_snpeff.vcf"
    output: "snpEff_summary.html"
    shell: "java -jar snpEff.jar -stats {output} {input}"

# Step 12: Plot alignment
rule plot_alignment:
    input: f"{PREFIX}.delta"
    output: "mummer_plot.png"
    shell: "mummer.py {input} -o {output}"
