configfile: "config.json"

import pandas as pd
samples = pd.read_table(config["samples"])#.set_index("Sample", drop=False)

rule all:
    input:
        fq1 = expand("results/subsampled_reads/{ena}_1.fastq.gz", ena = samples.ENA),
        fq2 = expand("results/subsampled_reads/{ena}_2.fastq.gz", ena = samples.ENA)

rule subsample_reads:
    input:
        fq1 = "data/reads/{ena}_1.fastq.gz",
        fq2 = "data/reads/{ena}_2.fastq.gz"
    output:
        fq1 = "results/subsampled_reads/{ena}_1.fastq.gz",
        fq2 = "results/subsampled_reads/{ena}_2.fastq.gz"
    shell:
        """
        set +o pipefail;\
        gunzip -c {input.fq1} | head -n 4000 | gzip -c > {output.fq1}
        gunzip -c {input.fq2} | head -n 4000 | gzip -c > {output.fq2}
        """
