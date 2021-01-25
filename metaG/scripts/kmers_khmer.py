#!/usr/bin/env python

"""
Exact k-mer counting w/ khmer

Based on: https://github.com/dib-lab/khmer/blob/fe0ce116456b296c522ba24294a0cabce3b2648b/examples/python-api/exact-counting.py
"""

##################################################
# IMPORT
##################################################
import khmer
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
def define_canonical_kmers(cg, nkmers):
    """
    Define canonical k-mers, i.e. exclude palindromic and rev. compl. k-mers

    Parameters
    ----------
    cg : khmer.Countgraph
        a k-mer countgraph
    nkmers : int
        number of all possible k-mers

    Returns
    -------
    set
        a set of canonical k-mers
    """
    canonical_kmers = set() # TODO Consider a sorting step to guarantee order of canonical kmers/kmer-hashes
    for i in range(nkmers):
        kmer = cg.reverse_hash(i)
        kmer_rev_comp = khmer.reverse_complement(kmer)
        # Store only the lexicographically *smaller* kmer or the palindromic kmer
        if kmer < kmer_rev_comp or kmer == kmer_rev_comp:
            canonical_kmers.add(i)
        else:
            continue
    return canonical_kmers

def build_countgraph(read, ksize, tabsize):
    """
    Build a k-mer countgraph from given nucleotide sequence, i.e. its k-mer profile

    Parameters
    ----------
    read : khmer read object
        nucleotide sequence (read, contig) to be used to build the count graph
    ksize : int
        k-mer size
    tabsize : int
        ???

    Returns
    -------
    khmer.Countgraph
        a k-mer countgraph
    """
    cg = khmer.Countgraph(ksize, tabsize, 1)
    cg.consume(read.sequence)
    return cg

def build_kmer_header(cg, canonical_kmers):
    """
    Create a k-mer profile header, i.e. list of k-mers

    Parameters
    ----------
    cg : khmer.Countgraph
        a k-mer countgraph
    canonical_kmers: set
        set of canonical k-mers

    Returns
    -------
    list
        list of k-mers
    """
    kmer_header = [cg.reverse_hash(kmer) for kmer in canonical_kmers]
    return kmer_header
        
def extract_kmer_profile(cg, canonical_kmers):
    """
    Extract k-mer counts, i.e. k-mer profile

    Parameters
    ----------
    cg : khmer.Countgraph
        a k-mer countgraph
    canonical_kmers: set
        set of canonical k-mers

    Returns
    -------
    list
        list of k-mer counts
    """
    kmer_profile = []
    for i in canonical_kmers:
        count = cg.get(i)
        if cg.get(i):
            kmer_profile.append(count)
        else:
            kmer_profile.append(0)
    return kmer_profile

##################################################
# VARS
##################################################
# Settings for exact k-mer counting
# Following https://github.com/dib-lab/khmer/blob/fe0ce116456b296c522ba24294a0cabce3b2648b/examples/python-api/exact-counting.py
KSIZE   = int(snakemake.wildcards.kmer)
NKMERS  = 4**KSIZE
TABSIZE = NKMERS + 10

# Initialize countgraph
CGRAPH = khmer.Countgraph(KSIZE, TABSIZE, 1)
# Get set of canonical kmers
CKMERS = define_canonical_kmers(CGRAPH, NKMERS)
logger.info("Number of canonical kmers for k = %i: %i" % (KSIZE, len(CKMERS)))

##################################################
# ANALYSIS
##################################################
with open(snakemake.output.kmers, "w") as ofile:
    # Header
    kmer_header = build_kmer_header(CGRAPH, CKMERS)
    ofile.write("ID\t%s\n" % "\t".join(kmer_header))
    # Counts
    for read in khmer.ReadParser(snakemake.input.fasta):
        # Create new Countgraph (start counting kmers from 0) and get the profile
        read_cgraph  = build_countgraph(read, KSIZE, TABSIZE)
        read_profile = extract_kmer_profile(read_cgraph, CKMERS)
        ofile.write("{}\t{}\n".format(read.name.split(" ")[0], "\t".join(map(str, read_profile))))

