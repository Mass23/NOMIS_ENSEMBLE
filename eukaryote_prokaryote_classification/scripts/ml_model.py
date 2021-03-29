#!/usr/bin/env python

"""
# Usage: python3 ml_model.py <file names>
"""

##################################################
# IMPORT
##################################################
import os
import glob
import logging
import pandas as pd
import numpy as np
import itertools
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split

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


##################################################
# SCRIPT
##################################################
data = pd.read_csv(snakemake.input.kmers, sep = '\t')

eukaryotic_contigs = pd.read_csv(snakemake.input.euk, sep = '\t', header = None)
eukaryotic_contigs = [i.split(' ')[0] for i in eukaryotic_contigs[0]]

prokaryotic_contigs = pd.read_csv(snakemake.input.pro, sep = '\t', header = None)
prokaryotic_contigs = [r[0] for i, r in prokaryotic_contigs.iterrows()]

data['domain'] = 'unknown'
for idx, row in data.iterrows():
    if (row.contig in prokaryotic_contigs) & (row.contig not in eukaryotic_contigs):
        data.loc[idx,'domain'] = 'prokaryota'
    if (row.contig not in prokaryotic_contigs) & (row.contig in eukaryotic_contigs):
        data.loc[idx,'domain'] = 'eukaryota'
    if (row.contig in prokaryotic_contigs) & (row.contig in eukaryotic_contigs):
        data.loc[idx,'domain'] = 'ambiguous'

print('Prokaryota: ' + str(len(data[data['domain'] == 'prokaryota'])))
print('Eukaryota: '  + str(len(data[data['domain'] == 'eukaryota'])))
print('Ambiguous: '  + str(len(data[data['domain'] == 'ambiguous'])))
print('Unknown: '    + str(len(data[data['domain'] == 'unknown'])))

full_set_pro = data[data['domain'] == 'prokaryota'].sample(n = 1500).fillna(0)
full_set_euk = data[data['domain'] == 'eukaryota' ].sample(n = 1500).fillna(0)

full_set = pd.concat([full_set_pro, full_set_euk], ignore_index=True)

print('Training neural network...')
x = full_set[[''.join(i) for i in itertools.product(['A', 'C', 'G', 'T'], repeat=4)]]
y = full_set['domain']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.33)

parameter_space = {
# https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
'hidden_layer_sizes': [(128,64,), (256,128,), (256,256,), (128,128), (512,256,), (256,256,256,), (128,128,128,)],
'activation': ['relu'],
'solver': ['adam'],
'alpha': [0.0003, 0.0004, 0.0005, 0.0006, 0.0007],
'learning_rate': ['invscaling']}

mlp = MLPClassifier(max_iter=100)

clf = GridSearchCV(mlp, parameter_space, cv=3,n_jobs=12)
clf.fit(x_train, y_train)

# Best paramete set : https://datascience.stackexchange.com/questions/36049/how-to-adjust-the-hyperparameters-of-mlp-classifier-to-get-more-perfect-performa
print('Best parameters found:\n', clf.best_params_)

# All results
means = clf.cv_results_['mean_test_score']
stds = clf.cv_results_['std_test_score']
for mean, std, params in zip(means, stds, clf.cv_results_['params']):
    print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))

y_true, y_pred = y_test , clf.predict(x_test)

print('Results on the test set:')
print(classification_report(y_true, y_pred))
