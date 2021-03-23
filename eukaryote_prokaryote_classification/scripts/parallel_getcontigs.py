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
from multiprocessing import Pool

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
def GetContig(arg_list):
    workerid, contig_size, genomes, n_contigs, category = arg_list
    with open("results/" + str(category) + "/" + str(category) + "_"  + str(workerid) + ".fasta", "w") as f:
        for i in range(0, n_contigs):
            genome_file = random.sample(genomes, 1)[0]
            with gzip.open(genome_file, "rt") as handle:
                parsed_fasta = list(SeqIO.parse(handle, "fasta"))
                genome_size = False
                while genome_size == False:
                    random_genome_scaffold = random.sample(parsed_fasta, 1)[0]
                    random_genome_sequence = str(random_genome_scaffold.seq)
                    if len(random_genome_sequence) > contig_size:
                        genome_size = True

                start_index = np.random.randint(0,len(random_genome_sequence)-contig_size)
                contigseq=str(random_genome_sequence[start_index:start_index+contig_size])
                f.write(">" + str(category) + "_workerid=" + str(workerid) + "_" + str(i) + "\n" + str(contigseq) + "\n")

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

#    count = 0
#    for i in range(sample_size):
#        count += 1
#        print('Sample number: ' + str(count))

#        euk_genome_file = random.sample(euk_genomes, 1)[0]
#        pro_genome_file = random.sample(pro_genomes, 1)[0]

#        euk_contig = GetContig(euk_genome_file, contig_size)
#        pro_contig = GetContig(pro_genome_file, contig_size)
       
    n_contigs = 1000
    pro_workers_list = [[i, contig_size, pro_genomes, n_contigs, 'pro'] for i in range(0,sample_size//n_contigs)]
    euk_workers_list = [[i, contig_size, euk_genomes, n_contigs, 'euk'] for i in range(0,sample_size//n_contigs)]
    
    n_cores = snakemake.params.cores
    
    with Pool(n_cores) as p:
        p.map(GetContig, pro_workers_list)
        p.map(GetContig, euk_workers_list)

if __name__ == "__main__":
    main()
