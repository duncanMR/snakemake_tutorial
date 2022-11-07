configfile: "config.json"

import pandas as pd
samples = pd.read_table(config["samples"])
ref_path = "data/reference/" + config["reference_name"] + ".fa.gz"
ref_path_unzipped = "data/reference/" + config["reference_name"] + ".fa"

#Run this to load data for troubleshooting
# import os
# os.chdir("/home/duncan/classes/ngs_pipeline")
# samples = pd.read_table("samples.tsv")
# config = {
#     "reference_name": "Pf3D7_v3",
#     "fastq_filename_template": "results/subsampled/{ENA}_{read}.fastq.gz",
#     "samples": "samples.tsv"
# }

include: "functions.snakefile"

rule all:
    input:
        html1 = expand("results/qc/{ID}_1_fastqc.html", ID = samples.Sample),
        html2 = expand("results/qc/{ID}_2_fastqc.html", ID = samples.Sample),
        multiqc = "results/qc/multiqc_report.html",
        aligned = expand("results/aligned/{ID}.bam", ID = samples.Sample),
        bedgraph = expand("results/coverage/{ID}.coverage.bedgraph", ID = samples.Sample),
        vcf = "results/joint_call.vcf.gz"

#------------------------------------------------------------------------------------------------
# QC
#------------------------------------------------------------------------------------------------

rule rename_reads:
    input:
        fq1 = lambda w: get_input_filename_from_ID( w.ID, read = 1 ),
        fq2 = lambda w: get_input_filename_from_ID( w.ID, read = 2 )
    output:
        fq1 = temp("results/tmp/{ID}_1.fastq.gz"),
        fq2 = temp("results/tmp/{ID}_2.fastq.gz")
    shell:
        """
        cp {input.fq1} {output.fq1}
        cp {input.fq2} {output.fq2}
        """

rule fastqc_reads:
    input:
        fq1 = rules.rename_reads.output.fq1,
        fq2 = rules.rename_reads.output.fq2
    output:
        html1 = "results/qc/{ID}_1_fastqc.html",
        html2 = "results/qc/{ID}_2_fastqc.html",
    params:
        outputdir = "results/qc"
    shell:
        """
        fastqc -q -o {params.outputdir} {input.fq1} {input.fq2}
        """

rule multiqc_reads:
    input:
        qc_dir = "results/qc"
    output:
        multiqc = "results/qc/multiqc_report.html"
    shell:
        """
        multiqc -f {input.qc_dir} -o {input.qc_dir}
        """

#------------------------------------------------------------------------------------------------
# ALIGNMENT
#------------------------------------------------------------------------------------------------

rule index_ref:
    input:
        ref = ref_path
    output:
        idx = multiext(ref_path, ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell:
        """
        bwa index {input.ref}
        """

rule align_reads:
    input:
        fq1 = rules.rename_reads.output.fq1,
        fq2 = rules.rename_reads.output.fq2,
        ref = ref_path
    output:
        sam = temp("results/tmp/{ID}.sam")
    params:
        read_group_spec = lambda w: get_read_group_line( w.ID )
    shell:
        """
        bwa mem -R {params.read_group_spec} \
        -t 2 -o {output.sam} {input.ref} {input.fq1} {input.fq2}
        """

rule convert_to_bam:
    input:
        sam = rules.align_reads.output.sam
    output:
        bam = temp("results/tmp/{ID}.bam")
    shell:
        """
        samtools view -b -o {output.bam} {input.sam}
        """

rule fix_mate_pairs:
    input:
        bam = rules.convert_to_bam.output.bam
    output:
        bam = temp("results/tmp/{ID}_fixmate.bam")
    shell:
        """
        samtools fixmate -m {input.bam} {output.bam}
        """

rule sort_reads:
    input:
        bam = rules.fix_mate_pairs.output.bam
    output:
        bam = temp("results/tmp/{ID}_sorted.bam")
    shell:
        """
        samtools sort -T "tmp" -o {output.bam} {input.bam}
        """

rule mark_dupes_and_index:
    input:
        bam = rules.sort_reads.output.bam
    output:
        bam = "results/aligned/{ID}.bam"
    shell:
        """
        samtools markdup -s {input.bam} {output.bam}
        samtools index {output.bam}
        """

rule compute_coverage:
    input:
        bam = rules.mark_dupes_and_index.output.bam
    output:
        bg = "results/coverage/{ID}.coverage.bedgraph"
    shell:
        """
        bedtools genomecov -bg -ibam {input.bam} > {output.bg}
        """

#------------------------------------------------------------------------------------------------
# VARIANT CALLING
#------------------------------------------------------------------------------------------------

rule unzip_and_index_ref:
    input:
        ref = ref_path
    output:
        ref = ref_path_unzipped,
        fai = ref_path_unzipped + ".fai"
    shell:
        """
        gunzip -c {input.ref} > {output.ref}
        samtools faidx {output.ref}
        """

rule call_variants:
    input:
        fa = rules.unzip_and_index_ref.output.ref,
        bam = expand("results/aligned/{ID}.bam", ID = samples.Sample)
    output:
        vcf = "results/joint_call.vcf.gz"
    params:
        regions = config["regions"]
    shell:
        """
        octopus -R {input.fa} -I {input.bam} --regions {params.regions} \
            -o {output.vcf}

        """
