#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 17:22:54 2025

@author: xies
"""

import pandas as pd
from os import path
import ncbi_genome_download as ngd

from Bio import Entrez

dirname = '/Users/xies/Library/CloudStorage/OneDrive-Stanford/Bioinformatics/Whi5 promoter/'

# Fungal BlastP
blastp = pd.read_csv(path.join(dirname,'Fungal BLASTP/blastp.txt'),sep='|',header=None)
# Subject is from GenBank accession number
blastp['GenBank'] = blastp[3].str.split('.',expand=True)[0]

#%%

# Fungal BlastN

blastn = pd.read_csv(path.join(dirname,'Fungal BLASTN/blastn.txt'),sep='|',header=None)
# Subject is from GenBank accession number
blastn['GenBank'] = blastn[3].str.split('.',expand=True)[0]
Entrez.email = 'xies@stanford.edu'
hits = []
for entry in blastn.iloc[0:10]['GenBank']:
    stream = Entrez.efetch(db='nucleotide',rettype='gb',id=entry,retmode='xml')
    rec = Entrez.read(stream)[0]
    hits.append(rec)

#%%

