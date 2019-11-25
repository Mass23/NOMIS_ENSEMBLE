import subprocess
import argparse

# python3 q2_create_visualise.py -i NAME
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--id', help='ID of the dataset', type=str, action = 'store', required = True)
parser.add_argument('-n', '--n_threads', help='Number of threads to use', type=int, action = 'store', required = True)
args = parser.parse_args()
dataset_id = args.id
n_threads = args.n_threads

# IMPORTATION ##################################################################
import_args = ['qiime', 'tools', 'import', '--type', "'SampleData[PairedEndSequencesWithQuality]'", '--input-path', dataset_id + '_manifest.csv', '--output-path', dataset_id + '_raw.qza',  '--input-format', 'PairedEndFastqManifestPhred33']
subprocess.call(' '.join(import_args), shell=True)

visualise_args = ['qiime', 'demux', 'summarize', '--i-data', dataset_id + '_raw.qza', '--o-visualization', dataset_id + '_raw.qzv']
subprocess.call(' '.join(visualise_args), shell = True)

trim_args = ['qiime', 'quality-filter', 'q-score',
  '--i-demux', coremicrobiome1'_raw.qza',
  '--p-min-quality', '15',
  '--p-quality-window', '4',
  '--p-min-length-fraction', '0.7',
  '--o-filtered-sequences', dataset_id + '_trimmed.qza',
  '--o-filter-stats', dataset_id + '_trimming_stats.qza']
subprocess.call(' '.join(trim_args), shell = True)

trim_stats_args = ['qiime', 'metadata', 'tabulate',
  '--m-input-file', dataset_id + '_trimming_stats.qza',
  '--o-visualization', dataset_id + '_trimming_stats.qzv']
subprocess.call(' '.join(trim_stats_args), shell = True)

# DENOISING ####################################################################
denoise_args = ['qiime','dada2', 'denoise-paired',
  '--i-demultiplexed-seqs', dataset_id + '_trimmed.qza',
  '--p-n-threads', str(n_threads),
  '--p-trim-left-f', '17',
  '--p-trunc-len-f', '300',
  '--p-trim-left-r', '21',
  '--p-trunc-len-r', '300',
  '--o-table', dataset_id + '_dada2_table.qza',
  '--o-representative-sequences', dataset_id + '_dada2_seqs.qza',
  '--o-denoising-stats', dataset_id + '_dada2_stats.qza']
subprocess.call(' '.join(denoise_args), shell = True)

dada2_stats = ['qiime', 'metadata', 'tabulate', '--m-input-file', dataset_id + '_dada2_stats.qza', '--o-visualization', dataset_id + '_dada2_stats.qzv']
subprocess.call(' '.join(dada2_stats), shell = True)

dada2_seqs = ['qiime', 'feature-table', 'tabulate-seqs', '--i-data', dataset_id + '_dada2_seqs.qza', '--o-visualization', dataset_id + '_dada2_seqs.qzv']
subprocess.call(' '.join(dada2_seqs), shell = True)

dada2_table = ['qiime', 'feature-table', 'summarize', '--i-table', dataset_id + '_dada2_table.qza', --o-visualization dataset_id + '_dada2_table.qzv']
subprocess.call(' '.join(dada2_table), shell = True)

dada2_table_export = ['qiime', 'tools', 'export',
  '--input-path', dataset_id + '_dada2_table.qza']
  '--output-path', dataset_id + '_dada2_table.biom']
subprocess.call(' '.join(dada2_table_export), shell = True)

# OTU CLUSTERING ###############################################################
vsearch_args = ['qiime', 'vsearch', 'cluster-features-de-novo',
  '--i-sequences', dataset_id + '_dada2_seqs.qza',
  '--i-table', dataset_id + '_dada2_table.qza',
  '--p-perc-identity', '0.97',
  '--p-threads', n_threads,
  '--o-clustered-table', dataset_id + '_vsearch_OTU97_table.qza',
  '--o-clustered-sequences', dataset_id + '_vsearch_OTU97_seqs.qza']
subprocess.call(' '.join(vsearch_args), shell = True)

vsearch_table = ['qiime', 'feature-table', 'tabulate-seqs',
  '--i-data', dataset_id + '_vsearch_OTU97_seqs.qza',
  '--o-visualization', dataset_id + '_vsearch_OTU97_seqs.qzv']
subprocess.call(' '.join(vsearch_table), shell = True)

vsearch_seqs = ['qiime', 'feature-table', 'summarize',
  '--i-table', dataset_id + '_vsearch_OTU97_table.qza'
  '--o-visualization', dataset_id + '_vsearch_OTU97_table.qzv']
subprocess.call(' '.join(vsearch_seqs), shell = True)

vsearch_table_export = ['qiime', 'tools', 'export',
  '--input-path', dataset_id + '_vsearch_OTU97_table.qza']
  '--output-path', dataset_id + '_vsearch_OTU97_table.biom']
subprocess.call(' '.join(vsearch_table_export), shell = True)
