#!/usr/bin/env python

import collections
from sklearn.cluster import DBSCAN
import time
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder

labels_arr = np.genfromtxt(open(str(snakemake.input.subset), "rb"), delimiter="\t", skip_header=1, usecols=0, dtype='str')
le = LabelEncoder()
labels = le.fit(labels_arr).transform(labels_arr)
unassigned_label = labels[np.where(labels_arr == 'Unassigned')[0]][0]
# unassigned_label = labels[np.where(labels_arr == Unassigned)[0]]

chunk_size=25000
data_array = np.empty([0, 8192])
count=0
for chunk in pd.read_csv(str(snakemake.input.subset), sep='\t', header=0, chunksize=chunk_size):
    count += 1
    print('Chunk ' + str(count))
    chunk_array = chunk[[i for i in chunk.columns if i not in ['Unassigned','ID','Sample']]].to_numpy()
    chunk_array = chunk_array / chunk_array.sum(axis=1, keepdims=True)
    data_array = np.append(data_array, chunk_array, axis=0)

def CalculateScore(labels, prediction,param):
    score_df = pd.DataFrame(columns=['Label','N_cluster','param', 'Sensitivity', 'Specificity'])
    for label in list(set(labels)):
        try:
            label_idxs = np.where(labels == label)[0]
            other_idxs = np.where(labels != label)[0]
            label_counts = collections.Counter(prediction[label_idxs])
            del label_counts[-1]
            other_counts = collections.Counter(prediction[other_idxs])
            del other_counts[-1]
            max_count = max(zip(label_counts.values(), label_counts.keys()))
            score_df = score_df.append({'Label':label,
                                        'Sensitivity':label_counts[max_count[1]] / len(label_idxs),
                                        'Specificity':(1-(other_counts[max_count[1]] / len(other_idxs))),
                                        'N_cluster':len(set(prediction[label_idxs])), 'param':param}, ignore_index=True)
        except:
            score_df = score_df.append({'Label':label,
                                        'Sensitivity':np.nan,
                                        'Specificity':np.nan,
                                        'N_cluster':np.nan, 'param':param}, ignore_index=True)        
    return(score_df)

#def CalculateScore(labels, prediction,param):
#    score_df = pd.DataFrame(columns=['Label','N_cluster','param', 'Sensitivity', 'Specificity'])
#
#    for label in list(set(labels)):
#        label_idxs = np.where(labels == label)[0]
#        other_idxs = np.where(labels != label)[0]
#        label_counts = collections.Counter(prediction[label_idxs])
#        del label_counts[-1]
#        other_counts = collections.Counter(prediction[other_idxs])
#        del other_counts[-1]
#        max_count = max(zip(label_counts.values(), label_counts.keys()))
#        score_df = score_df.append({'Label':label,
#                                    'Sensitivity':label_counts[max_count[1]] / len(label_idxs),
#                                    'Specificity':(1-(other_counts[max_count[1]] / len(other_idxs))),
#                                    'N_cluster':len(set(prediction[label_idxs])), 'param':param}, ignore_index=True)
#
#    return(score_df)

full_score_df = pd.DataFrame(columns=['Label','Sensitivity','Specificity','N_cluster','param'])

t0 = time.time()
for eps in np.logspace(-10,0,20,endpoint=True):
    print('eps = ' + str(eps))
    for i in range(1,6):
        print('  - iteration: ' + str(i))
        indexes = np.random.randint(data_array.shape[0], size=5000)
        data_subset = data_array[indexes,:]
        labels_subset = labels[indexes]
        clustering = DBSCAN(eps=eps, min_samples=2,n_jobs=10).fit(data_subset)
        full_score_df = full_score_df.append(CalculateScore([i for i in labels_subset if i != unassigned_label], clustering.labels_,eps))
    print('Done!')
    print(time.time() - t0)

full_score_df.to_csv(str(snakemake.output.score), sep = '\t', index = False)
