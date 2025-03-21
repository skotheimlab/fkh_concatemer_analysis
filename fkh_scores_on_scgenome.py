#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 17:46:21 2025

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

#23288841 is Sc assemblyID
yeast_genome = glob(path.join(dirname,f'Genomes/23288841/GCF_000146045.2_R64_genomic.fa'))
yeast_genome = Fasta(yeast_genome[0])

# All promoters from SCG
promoters = glob(path.join(dirname,f'all_S228C_ORFs/orf_genomic_1000_all.fasta'))
promoters = Fasta(promoters[0])
                     
Fkh = pd.read_csv(path.join(dirname,'FKH2_pwm.txt'),sep='\t',index_col=0,header=None).T
entropy = -np.sum(Fkh*np.log(Fkh),axis=1)
high_entropy_positions = [1,2,10,11,12]
Fkh = Fkh.drop(labels=high_entropy_positions).reset_index().drop(columns='index')
L = len(Fkh)

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


#%% Calculate Fkh score on all promoter sequence and all chromosomal sequence

prom_scores = {gene: get_pswm_score(Fkh,seq[:].seq) for gene,seq in promoters.items() }

mean_scores = {gene: score.mean() for gene,score in prom_scores.items()}
max_scores = {gene: score.max() for gene,score in prom_scores.items()}
sum_scores = {gene: score.sum() for gene,score in prom_scores.items()}

df = pd.DataFrame([mean_scores,max_scores,sum_scores]).T
df.columns = ['Mean','Max','Sum']

df.sort_values('Mean',ascending=False).head(10)

#%% Scores of 2mers

Fkh_2mer = pd.concat((Fkh,Fkh),ignore_index=True)
prom_scores = {gene: get_pswm_score(Fkh_2mer,seq[:].seq) for gene,seq in promoters.items() }

mean_scores = {gene: score.mean() for gene,score in prom_scores.items()}
max_scores = {gene: score.max() for gene,score in prom_scores.items()}
sum_scores = {gene: score.sum() for gene,score in prom_scores.items()}

df_2mer = pd.DataFrame([mean_scores,max_scores,sum_scores]).T
df_2mer.columns = ['Mean','Max','Sum']

df_2mer.sort_values('Max',ascending=False).head(10)

#%% Generate 'background' via permutation and calculate scores

nt_counts = {}

for chrom in yeast_genome.keys():
    if chrom == 'NC_027264.1': # MT DNA
        continue
    
    nts = pd.DataFrame(list(yeast_genome[chrom][:].seq.upper()),columns=['Nucleotide'])
    nts['Count'] = 1
    nt_counts[chrom] = nts.groupby('Nucleotide')['Count'].sum()

nt_counts = pd.DataFrame(nt_counts).T
nt_counts.loc['Genome'] = nt_counts.sum(axis=0)

nt_freq = nt_counts.loc['Genome'] / nt_counts.loc['Genome'].sum()

# Assume i.i.d. nucleotides and generate a random 1000mer sequence

from random import choices

rand_scores = np.array(
    [ get_pswm_score(Fkh,choices(nt_freq.index.values, weights=nt_freq.values,k=len(Fkh)+1))[0]
    for i in range(100000)])

#%%

plt.hist(rand_scores,100, weights=np.ones(len(rand_scores))/len(rand_scores))

plt.vlines(x = df.loc['YOR083W','Max'],ymin=0,ymax=0.02,color='r')
plt.xlabel('Fkh1 motif 2mer match scores')
plt.ylabel('Frequency')
plt.legend(['Max score within Whi5 promoter','Random sequences (S.c. nt frequency)'])

p_val = (rand_scores>df.loc['YOR083W','Max']).sum() / 1e6
plt.title(f'P-value = {p_val}')





