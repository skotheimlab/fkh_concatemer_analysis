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
import shlex

from Bio import Entrez, Seq, SeqIO
Entrez.email = 'xies@stanford.edu'

dirname = '/Users/xies/Library/CloudStorage/OneDrive-Stanford/Bioinformatics/Whi5/'

#%% Fungal BlastP

entrezIDs = pd.DataFrame(columns=['Protein','mRNA','Gene','TaxonID','Chr','Start','End','Strand','Assembly','Organism'])

#Manually put in Whi5 and Whi7
# entrezIDs['DAA09241'] = {'Protein':'DAA09241','mRNA':


blastp = pd.read_csv(path.join(dirname,'Fungal BLASTP/blastp.txt'),sep='|',header=None)
Nhead = len(blastp)

# Subject is from GenBank accession number
blastp['GenBank'] = blastp[3].str.split('.',expand=True)[0]
entrezIDs['Protein'] = blastp.iloc[0:Nhead]['GenBank']
entrezIDs.index = entrezIDs['Protein']

proteins = { protID :
    Entrez.read(Entrez.efetch(db='protein',rettype='gb',
                              id=protID,retmode='xml'))[0]
    for protID,entry in entrezIDs.iterrows()
}

#%% DB hierarchy depends on project format

def extract_coordinates_refseq(entry, entrezIDs):
    
    matches = findall('REFSEQ: accession (\S+)',entry['GBSeq_source-db'])

    transcriptID = matches[0].split('.')[0]
    entrezIDs.loc[protID,'mRNA'] = transcriptID
    stream = Entrez.efetch(db='nucleotide',rettpe='gb',id=transcriptID,retmode='xml')
    mRNA = Entrez.read(stream)[0]

    for feature in mRNA['GBSeq_feature-table']:
        if feature['GBFeature_key'] == 'gene':
           for qualifier in feature['GBFeature_quals']:
               if qualifier['GBQualifier_name'] == 'db_xref':
                   # Find valud genes
                   geneID = qualifier['GBQualifier_value'].split(':')[1]
                   entrezIDs.loc[protID,'Gene'] = geneID
                   
                   # Retrieve the gene record
                   stream = Entrez.efetch(db='Gene',rettype='gb',retmode='xml',
                                          id=geneID)
                   gene = Entrez.read(stream)[0]
                   
                   # Retrieve organism name
                   orgname = gene['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_orgname']['OrgName']['OrgName_name']['OrgName_name_binomial']['BinomialOrgName']
                   gen = orgname['BinomialOrgName_genus']
                   spe = orgname['BinomialOrgName_species']
                   entrezIDs.loc[protID,'Organism'] = gen + ' ' + spe
                   
                   print(geneID)

                   taxaID = gene['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_db'][0]['Dbtag_tag']['Object-id']['Object-id_id']                    
                   # Chromosome
                   chrID = gene['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']                    
                   locus = gene['Entrezgene_locus'][0]
                   # Find the correct commentary type: 'genomic'
                   for locus in gene['Entrezgene_locus']:
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

    return entrezIDs

def extract_coordinates_nonrefseq(entry, entrezIDs):
    
    protID = entry['GBSeq_locus']
    matches = findall('accession (\S+)',entry['GBSeq_source-db'])
    chromID = matches[0].split('.')[0]
    entrezIDs.loc[protID,'Chr'] = chromID
    # entrezIDs.loc[protID,'Assembly'] = chromID
    entrezIDs.loc[protID,'Organism'] = ' '.join(entry['GBSeq_organism'].split(' ')[0:2])
    
    chrm = Entrez.read(Entrez.efetch(db='nucleotide',id=chromID,rettype='gb',retmode='xml'))[0]
    
    for feature in chrm['GBSeq_feature-table']:
        if feature['GBFeature_key'] == 'source':
            for qual in feature['GBFeature_quals']:
                if 'db_xref' in qual.values():
                    entrezIDs.loc[protID,'TaxonID'] = qual['GBQualifier_value'].split(':')[1]
        if feature['GBFeature_key'] == 'CDS':
            for qual in feature['GBFeature_quals']:
                if 'protein_id' in qual.values():
                    if qual['GBQualifier_value'].startswith(protID):
                        print(protID)
                        start = int(feature['GBFeature_intervals'][0]['GBInterval_from'])
                        end = int(feature['GBFeature_intervals'][0]['GBInterval_to'])
                        if end > start:
                            entrezIDs.loc[protID,'Start'] = start
                            entrezIDs.loc[protID,'End'] = end
                            entrezIDs.loc[protID,'Strand'] = 'plus'
                        elif start > end:
                            entrezIDs.loc[protID,'Start'] = end
                            entrezIDs.loc[protID,'End'] = start
                            entrezIDs.loc[protID,'Strand'] = 'minus'

    # # Fetch the sequence
    # seq = Seq.Seq(chrm['GBSeq_sequence'])
    # # Format and write FASTA
    # SeqIO.SeqRecord(seq).write(path.join)
    # genome_dir = path.join(dirname,f'Genomes/{assemblyID}')
    # os.makedirs(genome_dir, exist_ok=True)

    return entrezIDs

    

for protID,entry in proteins.items():
    
    _matches = findall('REFSEQ: accession (\S+)',entry['GBSeq_source-db'])
    if len(_matches) > 0:
        entrezIDs = extract_coordinates_refseq(entry,entrezIDs)
        continue

    _matches = findall('accession (\S+)',entry['GBSeq_source-db'])
    if len(_matches) > 0:
        entrezIDs = extract_coordinates_nonrefseq(entry,entrezIDs)        
        
entrezIDs = entrezIDs.dropna(subset='Chr')
entrezIDs['Length'] = entrezIDs['End'] - entrezIDs['Start']

#%% Fetch genomes

for protID,entry in entrezIDs.iterrows():
    
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

#%% REFSEQ
# #% If there
# mRNAs = {}
# for protID,entry in proteins.items():
#     # transcriptID = entry['GBSeq_source-db'].split(' ')[2].split('.')[0]
#     matches = findall('REFSEQ: accession (\S+)',entry['GBSeq_source-db'])
#     if len(matches) > 0:
#         transcriptID = matches[0].split('.')[0]
#         entrezIDs.loc[protID,'mRNA'] = transcriptID
#         stream = Entrez.efetch(db='nucleotide',rettpe='gb',id=transcriptID,retmode='xml')
#         rec = Entrez.read(stream)[0]
#         mRNAs[protID] = rec
#         print(transcriptID)

# #% Find the gene records\
# genes = {}
# for protID,entry in mRNAs.items():
    
#     feature_table = entry['GBSeq_feature-table']
#     for feature in entry['GBSeq_feature-table']:
#         if feature['GBFeature_key'] == 'gene':
#            for qualifier in feature['GBFeature_quals']:
#                if qualifier['GBQualifier_name'] == 'db_xref':
#                    # Find valud genes
#                    geneID = qualifier['GBQualifier_value'].split(':')[1]
#                    entrezIDs.loc[protID,'Gene'] = geneID
                   
#                    # Retrieve the gene record
#                    stream = Entrez.efetch(db='Gene',rettype='gb',retmode='xml',
#                                           id=geneID)
#                    rec = Entrez.read(stream)[0]
#                    genes[protID] = rec
                   
#                    # Retrieve organism name
#                    orgname = rec['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_orgname']['OrgName']['OrgName_name']['OrgName_name_binomial']['BinomialOrgName']
#                    gen = orgname['BinomialOrgName_genus']
#                    spe = orgname['BinomialOrgName_species']
#                    entrezIDs.loc[protID,'Organism'] = gen + ' ' + spe
                   
#                    print(geneID)

# # Parse gene records to grab coordinates

# for protID,entry in genes.items():
    
#     print(protID)
#     # organisms
#     taxaID = entry['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_db'][0]['Dbtag_tag']['Object-id']['Object-id_id']
    
#     # Chromosome
#     chrID = entry['Entrezgene_gene-source']['Gene-source']['Gene-source_src-str1']
    
#     locus = entry['Entrezgene_locus'][0]
#     # Find the correct commentary type: 'genomic'
#     for locus in entry['Entrezgene_locus']:
#         if locus['Gene-commentary_type'].attributes['value'] == 'genomic':
#             if locus['Gene-commentary_accession'] == chrID:
            
#                 interval = locus['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
#                 start = interval['Seq-interval_from']
#                 end = interval['Seq-interval_to']
#                 strand = interval['Seq-interval_strand']['Na-strand'].attributes['value']
                        
#                 entrezIDs.loc[protID,'Chr'] = chrID
#                 entrezIDs.loc[protID,'Start'] = int(start)
#                 entrezIDs.loc[protID,'End'] = int(end)
#                 entrezIDs.loc[protID,'Strand'] = strand
#                 entrezIDs.loc[protID,'TaxonID'] = taxaID

# Calculate gene length









