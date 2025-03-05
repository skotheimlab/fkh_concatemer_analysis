#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 16:25:17 2025

@author: xies
"""

import numpy as np
from random import sample
import pandas as pd
from os import path
import os
import matplotlib.pyplot as plt
import seaborn as sb
from glob import glob
os.environ['PYSNEMBL_CACHE_DIR'] = '/Users/xies/Desktop/Ensembl cache'

from pyfaidx import Fasta, Faidx

dirname = '/Users/xies/Library/CloudStorage/OneDrive-Stanford/Bioinformatics/Whi5/'

entrezIDs = pd.read_csv(path.join(dirname,'blastP_genes.csv'),index_col=0)

#%%

entry = entrezIDs.iloc[9]

assembly = entry['Assembly']
genome_fasta = glob(path.join(dirname,f'Genomes/{assembly}/*.fa'))[0]

genome = Fasta(genome_fasta)
chrmID = entry['Chr'] + '.1'

if entry['Strand'] == 'plus':
    start = entry['Start']
    end = entry['End']
    seq = genome[chrmID][start:end+1]
    promoter = genome[chrmID][start-1000:start]
    promoter_rc = genome[chrmID][start-1000:start].reverse.complement
    promoter_uppercase = promoter.seq.upper()
    promoter_rc_uppercase = promoter_rc.seq.upper()
    rand_promoter_uppercase = ''.join(sample(promoter_uppercase,1000))
    
    
if entry['Strand'] == 'minus':
    start = entry['Start']
    end = entry['End']
    seq = genome[chrmID][start:end+1]
    promoter = genome[chrmID][end:end+1000].reverse.complement
    promoter_rc = promoter.reverse.complement
    promoter_uppercase = promoter.seq.upper()
    promoter_rc_uppercase = promoter_rc.seq.upper()
    rand_promoter_uppercase = ''.join(sample(promoter_uppercase,1000))

Fkh = pd.read_csv(path.join(dirname,'FKH1_pwm.txt'),sep='\t',index_col=0,header=None).T

entropy = -np.sum(Fkh*np.log(Fkh),axis=1)
high_entropy_positions = [1,2,10,11,12]
Fkh = Fkh.drop(labels=high_entropy_positions).reset_index().drop(columns='index')
L = len(Fkh)

scores = np.zeros(1000-L)
rc_scores = np.zeros(1000-L)
random_scores = np.zeros(1000-L)
for i in range(len(promoter_uppercase) - L):
    
    probes = {l:s for l in range(L) for s in list(promoter_uppercase[i:i+L])}
    rc_probes = {l:s for l in range(L) for s in list(promoter_rc_uppercase[i:i+L])}
    random_probes = {l:s for l in range(L) for s in list(rand_promoter_uppercase[i:i+L])}
    scores[i] = np.sum([Fkh.loc[k,v] for k,v in probes.items()])
    rc_scores[i] = np.sum([Fkh.loc[k,v] for k,v in rc_probes.items()])
    random_scores[i] = np.sum([Fkh.loc[k,v] for k,v in random_probes.items()])

plt.plot(scores)
plt.plot(rc_scores)
plt.plot(random_scores)


