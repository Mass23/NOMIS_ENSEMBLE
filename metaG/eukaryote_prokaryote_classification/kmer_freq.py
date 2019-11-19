import pandas as pd
from Bio import SeqIO
import argparse
from multiprocessing import Pool
import numpy as np
import gzip
import itertools

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fastafile'    , help='file to process, contigs fasta'          , type=str, action = 'store', required = True)
parser.add_argument('-k', '--kmersize'     , help='Kmer size to use for analysis'           , type=int, action = 'store', required = True)
parser.add_argument('-t', '--ThreadNumber' , help='Number of thread for parrallel computing', type=int, action = 'store', required = True)

args = parser.parse_args()
fasta_file = args.fastafile
k = args.kmersize
t = args.ThreadNumber

# Get Kmer spectrum of a sequence
def KmersFreqSeq(sequence, k):
    kmer_dict = dict()
    combs = [''.join(i) for i in itertools.combinations_with_replacement(['A','C','G','T'], 4)]
    for char in range(0,len(sequence) - k + 1):
        kmer = sequence[char:char + k]
        if kmer in combs:
            if kmer in kmer_dict.keys():
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    count = sum([kmer_dict[i] for i in kmer_dict.keys()])
    for kmer in kmer_dict.keys():
        kmer_dict[kmer] = kmer_dict[kmer] / count
    return(kmer_dict)

def GetGC(sequence):
    return((sequence.count('G') + sequence.count('C')) / (sequence.count('A') + sequence.count('C') + sequence.count('G') + sequence.count('T')))

def Process_file(record):
    dataset = pd.DataFrame()
    for seq in record:
        # Calculate tnfs
        record_output = dict()
        record_output.update(KmersFreqSeq(seq.seq, k))
        # Calculate sequence length and gc
        record_output.update({'contig' : str(seq.id), 'length' : len(str(seq.seq)),'gc_content' : GetGC(str(seq.seq))})
        dataset = dataset.append(record_output, ignore_index=True)
    return(dataset)

def Main():
    print('Parsing fasta file...')
    main_table = pd.DataFrame()
    if fasta_file.endswith('.gz'):
        handle = gzip.open(fasta_file, 'rt')
    else:
        handle = fasta_file
    data = list(SeqIO.parse(handle, 'fasta'))

    print('Calculating TNFs...')
    if t > 1:
        record_list = np.array_split(data,t)
        p = Pool(t)
        output = p.map(Process_file, record_list)
        dataset = pd.concat(output, sort=False)
    else:
        dataset = Process_file(data)

    print('Merging results...')
    dataset.to_csv(fasta_file + '_' + str(k) + 'mers.csv', sep='\t')

    print('Done!')

if __name__ == '__main__':
    Main()
