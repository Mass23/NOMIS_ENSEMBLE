"""
Author: Susheel Bhanu BUSI
Affiliation: ESB group LCSB UniLU
Date: [2019-06-02]
Modification on: [2021-01-11]
Run: snakemake -s Snakefile --use-conda --cores 10 -rp
Latest modification:
"""

import os
import glob
# import pandas 

configfile:"config.yaml"
DATA_DIR=config['data_dir']
RESULTS_DIR=config['results_dir']
SAMPLES=[line.strip() for line in open("sample_list", 'r')]    # if using a sample list instead of putting them in a config file
ENV_DIR=config['env_dir']
DEDUP_DIR=config['dedup_dir']

# print(f"number of cores: {workflow.cores}")

###########
rule all: 
    input: 
#        expand(os.path.join(DEDUP_DIR, "{sample}_R{reads}.fastq.gz"), reads=["1", "2"], sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "preproc/concat/merged_R{reads}.fastp.fastq.gz"), reads=["1", "2"]),
#        expand(os.path.join(RESULTS_DIR, "kraken2/{sample}.kraken.summary.out"), sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "headers/{sample}.headers.txt"), sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "extracted/{sample}_{read}.fastq"), sample=SAMPLES, read=["R1", "R2"]),
        expand(os.path.join(RESULTS_DIR, "assembly/ASSEMBLY.fasta"))
#        expand(os.path.join(RESULTS_DIR, "preproc/{sample}_R1.fastp.fastq.gz"), sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "qc/{sample}_{rid}.fastp_fastqc.html"), rid=["R1", "R2"], sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "nonpareil/{sample}/{sample}.fasta"), sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "metaxa2/{sample}_output/LEVEL-7_{sample}_rarefaction_out"), sample=SAMPLES)


#################
# Preprocessing #
#################
rule fastp_sr:
    input:
        r1=os.path.join(DATA_DIR, "coassembly/{sample}_Pool1_R1.fastq.gz"),
        r2=os.path.join(DATA_DIR, "coassembly/{sample}_Pool1_R2.fastq.gz")
    output:
        o1=os.path.join(RESULTS_DIR, "preproc/{sample}_R1.fastp.fastq.gz"),
        o2=os.path.join(RESULTS_DIR, "preproc/{sample}_R2.fastp.fastq.gz"),
        html=os.path.join(RESULTS_DIR, "preproc/{sample}_fastp.html"),
        json=os.path.join(RESULTS_DIR, "preproc/{sample}_fastp.json")
    log:
        out="logs/fastp.{sample}.out.log",
        err="logs/fastp.{sample}.err.log"
#    wildcard_constraints:
#        mtype="metag|metat"
    threads:
        config["fastp"]["threads"]
    conda:
        os.path.join(ENV_DIR, "fastp.yaml")
    message:
        "Preprocessing short reads: FastP"
    shell:
        "(date && fastp -l {config[fastp][min_length]} -i {input.r1} -I {input.r2} -o {output.o1} -O {output.o2} -h {output.html} -j {output.json} -w {threads} && date) 2> {log.err} > {log.out}"

rule fastqc_fastp_sr:
    input:
        os.path.join(RESULTS_DIR, "preproc/{sample}_{rid}.fastp.fastq.gz")
    output:
        html=os.path.join(RESULTS_DIR, "qc/{sample}_{rid}.fastp_fastqc.html"),
        zip=os.path.join(RESULTS_DIR, "qc/{sample}_{rid}.fastp_fastqc.zip")
    log:
        out="logs/fastqc.{sample}_{rid}.out.log",
        err="logs/fastqc.{sample}_{rid}.err.log"
    wildcard_constraints:
#        mtype="metag|metat",
        rid="|".join(["R1", "R2"])
    threads:
        config["fastqc"]["threads"]
    conda:
        os.path.join(ENV_DIR, "fastqc.yaml")
    message:
        "Preprocessing short reads: FastQC"
    shell:
        "(date && fastqc {config[fastqc][params]} -t {threads} -o $(dirname {output.html}) {input} && date) 2> {log.err} > {log.out}"


#################
# Preprocessing #
#################
rule deduplicate:
    input:
        r1=os.path.join(RESULTS_DIR, "preproc/{sample}_R1.fastp.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "preproc/{sample}_R2.fastp.fastq.gz")
    output:
        odup1=os.path.join(DEDUP_DIR, "{sample}_R1.fastq.gz"), 
        odup2=os.path.join(DEDUP_DIR, "{sample}_R2.fastq.gz")
    log:
        out="logs/dedup.{sample}.out.log",
        err="logs/dedup.{sample}.err.log"
    threads:
        config["clumpify"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bbmap.yaml")
    message:
        "Removing duplicate reads for easier downstream assembly"
    shell:
        "(date && clumpify.sh in={input.r1} in2={input.r2} out={output.odup1} out2={output.odup2} dupedist={config[clumpify][dupedist]} dedupe=t optical=t threads={threads} groups={config[clumpify][groups]} -Xmx{config[clumpify][memory]} && date) 2> {log.err} > {log.out}"

###################
# Read extraction #
###################
# Extracting only those reads classified as Eukaryotes
# Taxonomy with Kraken2
rule kraken2:
    input:
        dedup1=os.path.join(DEDUP_DIR, "{sample}_R1.fastq.gz"),
        dedup2=os.path.join(DEDUP_DIR, "{sample}_R2.fastq.gz"),
        database=config['kraken2']['db']
    output:
        rep=os.path.join(RESULTS_DIR, "kraken2/{sample}.kraken.report.txt"),
        summary=os.path.join(RESULTS_DIR, "kraken2/{sample}.kraken.summary.out")
    log:
        out="logs/{sample}.kraken2.out.log",
        err="logs/{sample}.kraken2.err.log"
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    message:
        "Kraken2 taxonomy"
    shell:
        "(date && kraken2 --threads {threads} --db {input.database} --use-names --confidence 0.5 --paired {input.dedup1} {input.dedup2} --gzip-compressed --output {output.summary} --report {output.rep} && date) 2> {log.err} > {log.out}"

rule headers:
    input:
        os.path.join(RESULTS_DIR, "kraken2/{sample}.kraken.summary.out")
    output:
        os.path.join(RESULTS_DIR, "headers/{sample}.headers.txt")
    log:
        out="logs/{sample}.headers.out.log",
        err="logs/{sample}.headers.err.log"
    message:
        "Extracting 'unclassified' read headers"
    shell:
        "(date && awk '{{if ($3 ~ /unclassified/ || $3 ~ /Eukaryota/) print $2}}' {input} > {output} && date) 2> {log.err} > {log.out}"

rule extract_reads:
    input:
        ids=rules.headers.output,
        dedup1=os.path.join(DEDUP_DIR, "{sample}_R1.fastq.gz"),
        dedup2=os.path.join(DEDUP_DIR, "{sample}_R2.fastq.gz")
    output:
        ex1=os.path.join(RESULTS_DIR, "extracted/{sample}_R1.fastq"),
        ex2=os.path.join(RESULTS_DIR, "extracted/{sample}_R2.fastq")
    log:
        out="logs/{sample}.reads.out.log",
        err="logs/{sample}.reads.err.log"
    conda:
        os.path.join(ENV_DIR, "seqtk.yaml")
    message:
        "Extracting reads from {wildcards.sample}"
    shell:
        "(date && seqtk subseq {input.dedup1} {input.ids} > {output.ex1} && seqtk subseq {input.dedup2} {input.ids} > {output.ex2} && date) 2> {log.err} > {log.out}"

# TODO
rule concatenate:
    input:
        read1=expand(os.path.join(RESULTS_DIR, "extracted/{sample}_R1.fastq"), sample=SAMPLES),
        read2=expand(os.path.join(RESULTS_DIR, "extracted/{sample}_R2.fastq"), sample=SAMPLES)
    output:
        or1=os.path.join(RESULTS_DIR, "preproc/concat/merged_R1.preprocessed.fastq"),
        or2=os.path.join(RESULTS_DIR, "preproc/concat/merged_R2.preprocessed.fastq")
    log:
        out="logs/concat.out.log",
        err="logs/concat.err.log"
    message:
        "Concatenating the reads for co-assembly"
    shell:
        "(date && cat {input.read1} > {output.or1} && cat {input.read2} > {output.or2} && date) 2> {log.err} > {log.out}"

rule compress:
    input:
        comp1=rules.concatenate.output.read1, 
        comp2=rules.concatenate.output.read2
    output:
        out1=os.path.join(RESULTS_DIR, "preproc/concat/merged_R1.preprocessed.fastq.gz"),
        out2=os.path.join(RESULTS_DIR, "preproc/concat/merged_R2.preprocessed.fastq.gz")
    message:
        "Compressing the fastq files"
    shell:
        "(date && gzip {input.comp1} > {output.out1} && gzip {input.comp2} > {output.out2} && date)"

############
# Assembly #
############
rule assembly_sr_megahit:
    input:
        sr1=rules.compress.output.out1,
        sr2=rules.compress.output.out2
    output:
        os.path.join(RESULTS_DIR, "assembly/ASSEMBLY.fasta")
    log:
        out="logs/megahit.out.log",
        err="logs/megahit.err.log"
    threads:
        config["megahit"]["threads"]
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    message:
        "Assembly: short reads: MEGAHIT"
    shell:
#        "ifiles1=$(echo \"{input.r1}\" | sed 's/ /,/g') && "
#        "ifiles2=$(echo \"{input.r2}\" | sed 's/ /,/g') && "
#        "echo ${{ifiles1}} && echo ${{ifiles2}} && "
##        "R1s=$(ls /scratch/users/sbusi/metaG_JULY_2020/coassembly/results/preproc/GL_*R1.fastp.fastq.gz | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])') && "
##        "R2s=$(ls /scratch/users/sbusi/metaG_JULY_2020/coassembly/results/preproc/GL_*R2.fastp.fastq.gz | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])') && "
#        "(date && megahit -1 $ifiles1 -2 $ifiles2 --kmin-1pass -m 60e+10 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t {threads} -o $(dirname {output})/tmp && "
        "(date && megahit -1 {input.sr1} -2 {input.sr2} --kmin-1pass -m 60e+10 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t {threads} -o $(dirname {output})/tmp && "
        "cd $(dirname {output}) && "
        "rsync -avP tmp/ . && "
        "ln -sf final.contigs.fa $(basename {output}) && "
        "rm -rf tmp/ && "
        "date) 2> {log.err} > {log.out}"
            

############
# Analyses #
############
rule Eukrep_fasta:
    input:
        rules.assembly_sr_megahit.output
    output:
        os.path.join(RESULTS_DIR, "eukrep/eukrep_contigs.fa")
    log:
        out="logs/eukrep.out.log",
        err="logs/eukrep.err.log"
    conda:
        os.path.join(ENV_DIR, "eukrep.yaml")
    threads:config["threads"]
    shell:
        "(date && EukRep -i {input} -o {output} --min 2000 -m strict && date) 2> {log.err} > {log.out}"


#####################################
##### NONPAREIL & METAXA2 ###########
#####################################                      
# TO DO: gunzip the fastq files for nonpareil and metaxa2  
rule nonpareil:
    input:
        os.path.join(DEDUP_DIR, "{sample}_R1.fastp.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "nonpareil/{sample}/{sample}.fasta")
    log:
        out="logs/nonpareil.{sample}.out.log",
        err="logs/nonpareil.{sample}.err.log"
    conda:
        os.path.join(ENV_DIR, "nonpareil.yaml")
    threads:config["threads"]
    shell:
        """
        (date && 
        zcat {input} | paste - - - - | awk 'BEGIN{{{{FS="\t"}}{{}}{{print ">"substr($1,2)"\n"$2}}}}' > {output} &&
        nonpareil -s {output} -T kmer -f fasta -b {wildcards.sample} &&
        date) 2> {log.err} > {log.out}
        """

rule metaxa2:
    input:
        mt1=os.path.join(DEDUP_DIR, "{sample}_R1.fastp.fastq.gz"),
        mt2=os.path.join(DEDUP_DIR, "{sample}_R2.fastp.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "metaxa2/{sample}_output/{sample}_metaxa_output.taxonomy.txt")
    log:
        out="logs/metaxa2.{sample}.out.log",
        err="logs/metaxa2.{sample}.err.log"
    conda:
        os.path.join(ENV_DIR, "metaxa2.yaml")
    threads:config["threads"]
    shell:
        "(date && "
        "metaxa2 --threads {threads} -f fastq -1 {input.mt1} -2 {input.mt2} -o $(dirname {output}) && "
        "date) 2> {log.err} > {log.out}"

rule rarefaction_metaxa2:
    input:
        rules.metaxa2.output
    output:
        os.path.join(RESULTS_DIR, "metaxa2/{sample}_output/LEVEL-7_{sample}_rarefaction_out")
    log:
        out="logs/rarefaction_metaxa2.{sample}.out.log",
        err="logs/rarefaction_metaxa2.{sample}.err.log"
    threads:config["threads"]
    conda:
        os.path.join(ENV_DIR, "metaxa2.yaml")
    shell:
        "(date && "
        "metaxa2_rf --threads {threads} -n 7 --resamples 10000 --scale 1000000 -i {input} -o {output} && "
        "date) 2> {log.err} > {log.out}"


