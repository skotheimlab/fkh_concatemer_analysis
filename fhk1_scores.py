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
entrezIDs = entrezIDs.rename(columns={'Protein.1':'Protein'})

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
    seq = genome[chrmID][start:end+1]
    
    # Check that the ORF starts with ATG, if so then extract 1kb upstream as promoter resion
    if entry['Strand'] == 'plus':
        if not seq[0:3].seq.upper() == 'ATG':
            print(f'{entry} does not begin with ATG')
        promoter = genome[chrmID][start - 1000:start]
    elif entry['Strand'] == 'minus':
        if not seq.reverse.complement[0:3].seq.upper() == 'ATG':
            print(f'{entry} does not begin with ATG')
        promoter = genome[chrmID][end:end+1000].reverse.complement
        
    promoter_ = promoter.seq.upper()
    promoter_rc = promoter.reverse.complement
    promoter_rc_ = promoter_rc.seq.upper()
    
    return promoter_, promoter_rc_

def get_pswm_score(mat,seq):
    L = len(mat)
    N = len(seq)
    return np.array([ np.sum([mat.loc[k,v] for k,v in enumerate( list(seq[i:i+L]) )]) for i in range(N-L) ])

#%%

scores_by_org = {}
scores_rc_by_org = {}
random_by_org = {}
for _,entry in entrezIDs.iterrows():
    
    species_name = entry['Organism']
    gene_entrez = entry['Gene']
    print(f'{species_name}')
    assembly = entry['Assembly']
    genome_fasta = glob(path.join(dirname,f'Genomes/{assembly}/*.fa'))[0]
    
    genome = Fasta(genome_fasta)
    chrmID = entry['Chr']
        
    promoter,promoter_rc = extract_promoter_seq(entry, genome)    
    
    
    scores = get_pswm_score(Fkh,promoter)
    scores_rc = get_pswm_score(Fkh,promoter_rc)
    
    scores_random = np.zeros((100,len(promoter) - len(Fkh)))
    for i in range(100):
        rand_seq = ''.join(sample(promoter, len(promoter)) )
        scores_random[i,:] = get_pswm_score(Fkh, rand_seq)
    mean_random_scores = np.mean(scores_random,axis=0)
    
    scores_by_org[entry['Protein']] = scores
    scores_rc_by_org[entry['Protein']] = scores_rc
    random_by_org[entry['Protein']] = mean_random_scores

#%% Generate plots

for org in scores_by_org.keys():
    
    entry = entrezIDs.loc[org]
    plt.plot(scores_by_org[org])
    plt.plot(scores_rc_by_org[org])
    plt.plot(random_by_org[org])
    
    species_name = entry['Organism']
    plt.title(species_name + '_' + org)
    plt.legend(['Promoter score','RevComp','Mean of permuted'])
    
    plt.savefig(path.join(path.join(dirname,f'Fkh1 scores/{species_name}_{org}.png')), dpi=199)
    plt.close()

#%% Stats

for org in entrezIDs.index:
    
    entrezIDs.loc[org,'Mean Fkh1'] = scores_by_org[org].mean()
    entrezIDs.loc[org,'Max Fkh1'] = np.max(scores_by_org[org])
    
    entrezIDs.loc[org,'Mean Fkh1 RC'] = scores_rc_by_org[org].mean()
    entrezIDs.loc[org,'Max Fkh1 RC'] = np.max(scores_rc_by_org[org])
    
    entrezIDs.loc[org,'Random Fkh1'] = random_by_org[org].mean()

entrezIDs['Max Fkh1 RC bgsub'] = entrezIDs['Mean Fkh1 RC'] - entrezIDs['Random Fkh1']
sb.catplot(entrezIDs,y='Max Fkh1 RC bgsub',x='Organism')
plt.tight_layout()
plt.xticks(rotation=90)
plt.show()



