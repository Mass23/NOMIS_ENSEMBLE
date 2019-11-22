import pandas as pd
from Bio import SeqIO
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastafile' , help='file to process, contigs fasta'          , type=str, action = 'store', required = True)
args = parser.parse_args()
fasta_file = args.fastafile

def Main():
    print('parsing fasta file...')
    data = list(SeqIO.parse(fasta_file, 'fasta'))
    gene_set = {}
    os.mkdir(fasta_file.replace('.ffn',''))

    print('creating the gene set dictionary...')
    count = 0
    for gene_annotation in data:
        count += 1
        print(str(count) + '/' + str(len(data)), end="\r")
        gene_name = ' '.join(str(gene_annotation.description).split(' ')[1:])
        if '(partial)' in gene_name:
            continue
        elif gene_name in gene_set.keys():
            gene_set['_'.join(gene_name.split(' '))].append(str(gene_annotation.id))
        else:
            gene_set['_'.join(gene_name.split(' '))] = [str(gene_annotation.id)]

    print('writing fasta files and compiling stats...')
    stats = pd.DataFrame(columns=['gene','copy_number'])
    count = 0
    for gene in gene_set.keys():
        count += 1
        print(str(count) + '/' + str(len(gene_set.keys())), end="\r")

        annotation_ids = gene_set[gene]
        stats = stats.append({'gene': gene, 'copy_number': len(annotation_ids)}, ignore_index = True)
        gene_data = [i for i in data if i.id in annotation_ids]
        for copy in gene_data:
            copy.id = fasta_file.replace('.ffn','') + '_' + str(copy.description).split(' ')[0]
            copy.description = ''
        SeqIO.write(gene_data, fasta_file.replace('.ffn','') + '/' + str(gene).replace('/','-') + '.fasta', 'fasta')

    stats.to_csv(fasta_file.replace('.ffn','') + '/' + fasta_file.replace('.ffn','') + '_stats.csv',index=False)

if __name__ == '__main__':
    Main()
