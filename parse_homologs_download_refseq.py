#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 17:22:54 2025

@author: xies
"""

import pandas as pd
from os import path
import os
from re import findall

from Bio import Entrez

dirname = '/Users/xies/Library/CloudStorage/OneDrive-Stanford/Bioinformatics/Whi5/'

#%% Fungal BlastP

entrezIDs = pd.DataFrame(columns=['Protein','mRNA','Gene','TaxonID','Chr','Start','End','Strand','Assembly','Organism'])

blastp = pd.read_csv(path.join(dirname,'Fungal BLASTP/blastp.txt'),sep='|',header=None)
Nhead = len(blastp)

# Subject is from GenBank accession number
blastp['GenBank'] = blastp[3].str.split('.',expand=True)[0]
Entrez.email = 'xies@stanford.edu'

proteins = { entry['GenBank'] :
    Entrez.read(Entrez.efetch(db='protein',rettype='gb',retmode='xml',
                              id=entry['GenBank']))[0]
                for _,entry in blastp.iloc[0:Nhead].iterrows() }

entrezIDs['Protein'] = blastp.iloc[0:Nhead]['GenBank']
entrezIDs.index = entrezIDs['Protein']

#% Kick it up to the mRNA

mRNAs = {}
for protID,entry in proteins.items():
    # transcriptID = entry['GBSeq_source-db'].split(' ')[2].split('.')[0]
    matches = findall('REFSEQ: accession (\S+)',entry['GBSeq_source-db'])
    if len(matches) > 0:
        transcriptID = matches[0].split('.')[0]
        entrezIDs.loc[protID,'mRNA'] = transcriptID
        stream = Entrez.efetch(db='nucleotide',rettpe='gb',id=transcriptID,retmode='xml')
        rec = Entrez.read(stream)[0]
        mRNAs[protID] = rec
        print(transcriptID)

#% Find the gene records

genes = {}
for protID,entry in mRNAs.items():
    
    feature_table = entry['GBSeq_feature-table']
    for feature in entry['GBSeq_feature-table']:
        if feature['GBFeature_key'] == 'gene':
           for qualifier in feature['GBFeature_quals']:
               if qualifier['GBQualifier_name'] == 'db_xref':
                   # Find valud genes
                   geneID = qualifier['GBQualifier_value'].split(':')[1]
                   entrezIDs.loc[protID,'Gene'] = geneID
                   
                   # Retrieve the gene record
                   stream = Entrez.efetch(db='Gene',rettype='gb',retmode='xml',
                                          id=geneID)
                   rec = Entrez.read(stream)[0]
                   genes[protID] = rec
                   
                   # Retrieve organism name
                   orgname = rec['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_orgname']['OrgName']['OrgName_name']['OrgName_name_binomial']['BinomialOrgName']
                   gen = orgname['BinomialOrgName_genus']
                   spe = orgname['BinomialOrgName_species']
                   entrezIDs.loc[protID,'Organism'] = gen + ' ' + spe
                   
                   print(geneID)

# Parse gene records to grab coordinates

for protID,entry in genes.items():
    
    print(protID)
    # organisms
    taxaID = entry['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_db'][0]['Dbtag_tag']['Object-id']['Object-id_id']
    
    # Chromosome
    chrID = entry['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
    
    locus = entry['Entrezgene_locus'][0]
    # Find the correct commentary type: 'genomic'
    for locus in entry['Entrezgene_locus']:
        if locus['Gene-commentary_type'].attributes['value'] == 'genomic':
            if locus['Gene-commentary_accession'] == chrID:
            
                interval = locus['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
                start = interval['Seq-interval_from']
                end = interval['Seq-interval_to']
                strand = interval['Seq-interval_strand']['Na-strand'].attributes['value']
                        
                entrezIDs.loc[protID,'Chr'] = chrID
                entrezIDs.loc[protID,'Start'] = int(start)
                entrezIDs.loc[protID,'End'] = int(end)
                entrezIDs.loc[protID,'Strand'] = strand
                entrezIDs.loc[protID,'TaxonID'] = taxaID

# Calculate gene length

entrezIDs = entrezIDs.dropna(subset='Gene')
entrezIDs['Length'] = entrezIDs['End'] - entrezIDs['Start']

#%% Fetch genomes

import shlex

for protID,entry in entrezIDs.iloc[0:].iterrows():
    
    stream = Entrez.esearch(db='assembly',
        term=entry['TaxonID']+'[taxid]', rettype='docsum',retmode='xml')
    record = Entrez.read(stream)
    assemblyID = record['IdList'][0]
    entrezIDs.loc[protID,'Assembly'] = assemblyID
    
    assemb = Entrez.read(Entrez.efetch(id = assemblyID,
                                       db='assembly',rettype='docsum',retmode='xml'))
    
    for docsum in assemb["DocumentSummarySet"]["DocumentSummary"]:
        ftp_path = docsum["FtpPath_RefSeq"]
        if ftp_path:
            fasta_url = ftp_path + "/" + ftp_path.split("/")[-1] + "_genomic.fna.gz"
            
            # Download the genomes
            genome_dir = path.join(dirname,f'Genomes/{assemblyID}')
            os.makedirs(genome_dir, exist_ok=True)
            os.system(f"wget -c --read-timeout=5 --tries=0 {shlex.quote(fasta_url)} -P {genome_dir}")


entrezIDs.to_csv(path.join(dirname,'blastP_genes.csv'))










