#!/usr/bin/env python
# Script used to extract individual cluster sequenceIDs into separate files

import os
import pandas as pd

df=pd.read_csv("unassigned_clusters_EUCI.tsv", sep='\t')
df.drop(df.columns[0], axis=1, inplace=True)    # remove index column
df1 = (df.assign(Sequences = df['Sequences'].str.split(',')).explode('Sequences').reset_index(drop=True)) # separate single row comma-separated values into columns as new rows

gp = df1.groupby('ClusterID')
gp_edited = gp[['Sequences']]
gp_edited.apply(lambda x: x.to_csv(str(x.name) + '.txt', sep='\t', header=False, index=False)) # output multiple files based on common values in column (https://stackoverflow.com/questions/37216230/output-multiple-files-based-on-column-value-python-pandas)
