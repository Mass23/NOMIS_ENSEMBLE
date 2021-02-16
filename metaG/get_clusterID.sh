#!/usr/bin/env python

import os
import pandas as pd

df=pd.read_csv("unassigned_clusters_EUCI.tsv", sep='\t')
df.drop(df.columns[0], axis=1, inplace=True)    # remove index column
df1 = (df.assign(Sequences = df['Sequences'].str.split(',')).explode('Sequences').reset_index(drop=True)) # separate single row comma-separated values into columns as new rows

gp = df1.groupby('ClusterID')
gp_edited = gp[['Sequences']]
gp_edited.apply(lambda x: x.to_csv(str(x.name) + '.txt', sep='\t', header=False, index=False)) # output multiple files based on common values in column (https://stackoverflow.com/questions/37216230/output-multiple-files-based-on-column-value-python-pandas)

# saving old cluster files before creating new ones
# mkdir orig_cluster_files && mv *.txt orig_cluster_files

df=pd.read_csv("Unassigned_only_min_7_samples.tsv", sep='\t')
df.drop(df.columns[0], axis=1, inplace=True)
df.drop(df.columns[1], axis=1, inplace=True)
new_df = df[~df.Sequence.str.contains("K", na=False)]

gp = new_df.groupby('ClusterID')
gp_edited = gp[['Sequence']]
gp_edited.apply(lambda x: x.to_csv(str(x.name) + '.txt', sep='\t', header=False, index=False)) # output multiple files based on common values in column (https://stackoverflow.com/questions/37216230/output-multiple-files-based-on-column-value-python-pandas)

# for those with min-9 sequnces and at least in 2 samples
df=pd.read_csv("clusters_min_9_seq_2_samp.tsv", sep='\t')
df.drop(df.columns[0], axis=1, inplace=True)
gp=df.groupby('ClusterID')
gp_edited=gp[['Sequence']]
gp_edited.apply(lambda x: x.to_csv(str(x.name) + '.txt', sep='\t', header=False, index=False))

# Did the following in BASH to get the file clusters. Commented here due to python script
# for the min_7 and clusters_min_9_seq_2 samples, have to rename all the files to have an underscore (_)
# renaming part of filename
# rename "Cluster" "Cluster_" *.txt
# ls -1 *.txt | sed 's/Cluster_//g' | sed 's/.txt//g' > cluster_list
