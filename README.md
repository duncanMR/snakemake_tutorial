# Snakemake pipeline for Plasmodium NGS data

Written as part of a course on next-generation sequencing bioinformatics at the Wellcome Centre for Human Genetics, Oxford University. To install:
```
$ git clone https://github.com/duncanMR/snakemake_tutorial snakemake_tutorial
```
Put the data to be analysed in `snakemake_tutorial/data/reads` and update `config.json` with your parameters of interest. You can then run the pipeline as follows:
```
$ cd snakemake_tutorial
$ snakemake -j1 -s pipelines/analysis.snakefile
```

`
