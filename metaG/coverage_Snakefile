"""
Author: Susheel Bhanu BUSI
Affiliation: ESB group LCSB UniLU
Date: [2021-01-31]
Run: snakemake -s coverage_Snakefile --use-conda --cores 45 -rp
Latest modification:
"""

import os, fnmatch
import glob
import pandas as pd

configfile:"coverage_config.yaml"
DATA_DIR=config['data_dir']
RESULTS_DIR=config['results_dir']
FASTQ_DIR=config['fastq_dir']
ASSEMBLY_DIR=config['assembly_dir']
SAMPLES=[line.strip() for line in open("sample_list", 'r')]    # if using a sample list instead of putting them in a config file

###########
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "coverage/{sample}_covstats.txt"), sample=SAMPLES)

################################
# rules for files and analyses #
################################
rule coverage:
    input:
        fa=os.path.join(ASSEMBLY_DIR, "{sample}_prokka.fna"),
        r1=os.path.join(FASTQ_DIR, "{sample}/{sample}_pass_1.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}/{sample}_pass_2.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "coverage/{sample}_covstats.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.coverage.log")
#    params:
#        tmp_dir=config["coverm"]["tmp_dir"]
    conda:
        os.path.join("envs/coverm.yaml")
    message:
        "calculating coverage for for {wildcards.sample}"
    shell:
        "(date && coverm contig -1 {input.r1} -2 {input.r2} --reference {input.fa} --output-file {output} -t {threads} && date) &> {log}"
