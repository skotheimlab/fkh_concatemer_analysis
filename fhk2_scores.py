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
from re import findall
import matplotlib.pyplot as plt
import seaborn as sb
from glob import glob

from Bio import Seq, SeqIO
from pyfaidx import Fasta

dirname = '/Users/xies/Library/CloudStorage/OneDrive-Stanford/Bioinformatics/Whi5/'

entrezIDs = pd.read_csv(path.join(dirname,'blastP_genes.csv'),index_col=0)
entrezIDs = entrezIDs.rename(columns={'Protein.1':'Protein'})

Fkh = pd.read_csv(path.join(dirname,'FKH2_pwm.txt'),sep='\t',index_col=0,header=None).T
entropy = -np.sum(Fkh*np.log(Fkh),axis=1)
high_entropy_positions = [1,2,10,11,12]
Fkh = Fkh.drop(labels=high_entropy_positions).reset_index().drop(columns='index')
L = len(Fkh)

# Drop suppressed entries
entrezIDs = entrezIDs.drop('EUN24805')

#%% Extract promoter sequence and calculate motif scores

def extract_promoter_seq(entry,genome):
    
    # Figure out the exact key
    chrmID = pd.DataFrame(genome.keys())[pd.DataFrame(genome.keys())[0].str.startswith(entry['Chr'])].iloc[0].values[0]
    
    start = int(entry['Start'])
    end = int(entry['End'])
    
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
promoters = {}
for _,entry in entrezIDs.iterrows():
    
    species_name = entry['Organism']
    gene_entrez = entry['Gene']
    print(f'{species_name}')
    assembly = entry['Assembly']
    genome_fasta = glob(path.join(dirname,f'Genomes/{assembly}/*.fa'))[0]
    
    genome = Fasta(genome_fasta)
    chrmID = entry['Chr']
        
    promoter,promoter_rc = extract_promoter_seq(entry, genome)
    promoters[entry['Protein']] = promoter
    
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

import pickle as pkl
with open(path.join(dirname,'promoters.pkl'),'wb') as f:
    pkl.dump(promoters,f)
with open(path.join(dirname,'scores_fkh2.pkl'),'wb') as f:
    pkl.dump(scores_by_org,f)
with open(path.join(dirname,'scores_rc_fkh2.pkl'),'wb') as f:
    pkl.dump(scores_rc_by_org,f)
with open(path.join(dirname,'scores_random_fkh2.pkl'),'wb') as f:
    pkl.dump(random_by_org,f)
    
SeqIO.write([SeqIO.SeqRecord(id=key, seq = Seq.Seq(seq)) for key,seq in promoters.items()],
            path.join(dirname,'promoters.fa'),format='fasta')

#%% Read in data

import pickle as pkl
with open(path.join(dirname,'promoters.pkl'),'rb') as f:
    promoters = pkl.load(f)
with open(path.join(dirname,'scores_fkh2.pkl'),'rb') as f:
    scores_by_org = pkl.load(f)
with open(path.join(dirname,'scores_rc_fkh2.pkl'),'rb') as f:
    scores_rc_by_org = pkl.load(f)
with open(path.join(dirname,'scores_random_fkh2.pkl'),'rb') as f:
    random_by_org = pkl.load(f)

#%%

from scipy import interpolate

for org in scores_by_org.keys():
    
    entry = entrezIDs.loc[org]
    score = scores_by_org[org]
    f = interpolate.make_smoothing_spline(range(len(score)),score)
    plt.plot(f(range(len(score))),alpha = 0.5)
    score_rc = scores_rc_by_org[org]
    f = interpolate.make_smoothing_spline(range(len(score_rc)),score_rc)
    plt.plot(f(range(len(score))),alpha = 0.5)
    plt.hlines(y=random_by_org[org].mean(),xmin=0,xmax=1000, alpha = 0.5)
    
    species_name = entry['Organism']
    plt.title(species_name + '_' + org)
    plt.legend(['Promoter score','RevComp','Mean of permuted'])
    
    plt.savefig(path.join(path.join(dirname,f'Fkh2 scores/{species_name}_{org}.png')), dpi=199)
    plt.close()



#%% Stats

for org in entrezIDs.index:
    
    score = scores_by_org[org]
    f = interpolate.make_smoothing_spline(range(len(score)),score)
    score_rc = scores_rc_by_org[org]
    f_rc = interpolate.make_smoothing_spline(range(len(score_rc)),score_rc)
    entrezIDs.loc[org,'Fkh2'] = max(f(range(len(score))).max() - random_by_org[org].mean(),
                                    f_rc(range(len(score))).max() - random_by_org[org].mean())
    
    entrezIDs.loc[org,'Fkh2'] = f(range(len(score))).max() - random_by_org[org].mean()
        
    # entrezIDs.loc[org,'Random Fkh1'] = random_by_org[org].mean()

sb.catplot(entrezIDs,y='Fkh2',x='Organism')
plt.tight_layout()
plt.xticks(rotation=90)
plt.show()


#%% Do some string op

for protID in entrezIDs.index:
    
    prom = promoters[protID]
    entrezIDs.loc[protID,'AACAA'] = len(findall('AACAA',prom))
    entrezIDs.loc[protID,'AACAAAACAA'] = len(findall('AACAAAACAA',prom))
    
    entrezIDs.loc[protID,'TTGTT'] = len(findall('TTGTT',prom))
    entrezIDs.loc[protID,'TTGTTTTGTT'] = len(findall('TTGTTTTGTT',prom))

entrezIDs.to_csv(path.join(dirname,'results_fkh2.csv'))

