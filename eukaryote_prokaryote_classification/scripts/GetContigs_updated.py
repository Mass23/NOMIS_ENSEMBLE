#!/usr/bin/env python

"""
# Usage python3 GetContigs_updated.py -n 500000 -k 5 -s 5000
"""

##################################################
# IMPORT
##################################################
import os
import glob
from Bio import SeqIO
import argparse
import pandas as pd
import random
import numpy as np
import itertools
import gzip
import logging

##################################################
# LOGGING
##################################################
# logger
logging.basicConfig(
    filename=str(snakemake.log),
    filemode="w",
    level=logging.DEBUG,
    format='[%(asctime)s] %(name)s %(levelname)s: %(message)s'
)
logger = logging.getLogger(__file__)

##################################################
# FUNCTIONS
##################################################
def GetContig(genome_file, contig_size):
    with gzip.open(genome_file, "rt") as handle:
        parsed_fasta = list(SeqIO.parse(handle, "fasta"))
        genome_size = False
        while genome_size == False:
            random_genome_scaffold = random.sample(parsed_fasta, 1)[0]
            random_genome_sequence = str(random_genome_scaffold.seq)
            if len(random_genome_sequence) > contig_size:
                genome_size = True

        start_index = np.random.randint(0,len(random_genome_sequence)-contig_size)
        return(str(random_genome_sequence[start_index:start_index+contig_size]))

def main():
    refseq_path = snakemake.params.refseq_path
    sample_size = snakemake.params.sample_size
    contig_size = snakemake.params.contig_size
    k = snakemake.params.kmer

    pro_dirs = ['archaea','bacteria']
    euk_dirs  = ['fungi','protozoa','invertebrate','vertebrate_mammalian','plant','vertebrate_other']

    # flattened list of all genomes in a group
    pro_genomes = [item for sublist in [glob.glob(refseq_path + folder + '/*/*.fna.gz') for folder in pro_dirs] for item in sublist]
    euk_genomes = [item for sublist in [glob.glob(refseq_path + folder + '/*/*.fna.gz') for folder in euk_dirs] for item in sublist]
    print('Prokaryotic genomes:')
    print(pro_genomes)
    print('Eukaryotic genomes:')
    print(euk_genomes)

    count = 0
    for i in range(sample_size):
        count += 1
        print('Sample number: ' + str(count))

        euk_genome_file = random.sample(euk_genomes, 1)[0]
        pro_genome_file = random.sample(pro_genomes, 1)[0]

        euk_contig = GetContig(euk_genome_file, contig_size)
        pro_contig = GetContig(pro_genome_file, contig_size)
        
        os.makedirs(str("results/pro"), exist_ok=True)
        os.makedirs(str("results/euk"), exist_ok=True)

        with open("results/pro/pro_" + str(count) + ".fasta", "w") as pro_file:
            pro_file.write(">pro_" + str(count) + "\n" + str(pro_contig) + "\n")
        
        with open("results/euk/euk" + str(count) + ".fasta", "w") as euk_file:
            euk_file.write(">euk_" + str(count) + "\n" + str(euk_contig) + "\n")


if __name__ == "__main__":
    main()
