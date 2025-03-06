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

Fkh = pd.read_csv(path.join(dirname,'FKH1_pwm.txt'),sep='\t',index_col=0,header=None).T
entropy = -np.sum(Fkh*np.log(Fkh),axis=1)
high_entropy_positions = [1,2,10,11,12]
Fkh = Fkh.drop(labels=high_entropy_positions).reset_index().drop(columns='index')
L = len(Fkh)

#%% Extract promoter sequence and calculate motif scores

def extract_promoter_seq(entry,genome):
    
    # Figure out the exact key
    chrmID = pd.DataFrame(genome.keys())[pd.DataFrame(genome.keys())[0].str.startswith(entry['Chr'])].iloc[0].values[0]
    
    start = entry['Start']
    end = entry['End']
    # seq = genome[chrmID][start:end+1]
    if entry['Strand'] == 'plus':
        promoter = genome[chrmID][start - 1000:start]
    elif entry['Strand'] == 'minus':
        promoter = genome[chrmID][end:end+100].reverse.complement
        
    promoter_ = promoter.seq.upper()
    promoter_rc = promoter.reverse.complement
    promoter_rc_ = promoter_rc.seq.upper()
    
    return promoter_, promoter_rc_

entry = entrezIDs.iloc[9]
assembly = entry['Assembly']
genome_fasta = glob(path.join(dirname,f'Genomes/{assembly}/*.fa'))[0]

genome = Fasta(genome_fasta)
chrmID = entry['Chr']
    
promoter,promoter_rc = extract_promoter_seq(entry, genome)    


def get_pswm_score(mat,seq):
    L = len(mat)
    N = len(seq)
    return np.array([ np.sum([mat.loc[k,v] for k,v in enumerate( list(seq[i:i+L]) )]) for i in range(N-L) ])

scores = get_pswm_score(Fkh,promoter)
scores_rc = get_pswm_score(Fkh,promoter_rc)

scores_random = np.zeros((100,len(promoter) - len(Fkh)))
for i in range(100):
    rand_seq = ''.join(sample(promoter, len(promoter)) )
    scores_random[i,:] = get_pswm_score(Fkh, rand_seq)
mean_random_scores = np.mean(scores_random,axis=0)

plt.plot(scores)
plt.plot(scores_rc)
plt.plot(mean_random_scores)
plt.title('')





