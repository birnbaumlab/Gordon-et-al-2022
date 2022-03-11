#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 09:02:33 2021

Maps barcode frequencies generated by Illumina sequencing onto long-read PacBio data

@author: khloegordon
"""

# -- import standard libraries -- #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import xlsxwriter
import editdistance
from datetime import datetime
        
#%% -- SETTINGS -- %%#

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'

dataset_filenames = ['Unselected.fastq.xlsx',
                     'EGFP+ sorted.fastq.xlsx',
                     '1st CD69+ sorted.fastq.xlsx',
                     '2nd CD69+ sorted.fastq.xlsx',
                     '3rd CD69+ sorted.fastq.xlsx',
                     '1st CD69+PD-1- sorted.fastq.xlsx',
                     '2nd CD69+PD-1- sorted.fastq.xlsx']
                     

settings = {'selection': ['CD69+','CD69+PD-1-'], 
            'round': [[0,1,2,3,4],[3,4]],
            'sheets': ['Unselected','EGFP+ sorted','CD69+ R1','CD69+ R2',
                       'CD69+ R3','CD69+PD-1- R1','CD69+PD-1- R2']}

#selections = ['Unselected','CD69+ R1','CD69+ R2','CD69+ R3','CD69+PD-1- R2']
# selections entries should match sheet names in dataset_filename spreadsheet
reference_filenames = ['reference_0009.xlsx',
                       'reference_0010.xlsx']
ICDfile = 'ICD lists and sequences updated 06.24.21.xlsx'
# ICDfile should have ICD names in column A and ICD family in column B
tolerance = 2 
# levenshtein edit distance tolerance for matching PacBio and Illumina barcodes

#%% -- function -- %%#

def domain_mapping(ILdf,PBdict,ICDs,tolerance,positions = 3):
    '''
    Maps ICD domains and metadata from PacBio reference onto Illumina barcodes 
    using levenshtein edit distance as tolerance for matching. positions = number of domains
    '''
                
    mapped = 0 # counter for matched sequences
    ILdf['Consensus Domains'] = ''   # create Consensus column in Illumina for PacBio Consensus domain combinations      

    bcs = ILdf['BC'].astype(str).values.tolist() # extract Illumina barcodes
    
    # compare Illumina and PacBio barcodes
    matches = ['']*len(ILdf.index)
    for i,bc in enumerate(bcs):
        
        if i % 100 == 0: 
            print('{} of {} barcodes mapped'.format(i,len(ILdf.index)))
            
        if bc in PBdict: 
            # if Illumina barcode is an exact match for a PacBio barcode
            matches[i] = PBdict[bc]
            #ILdf.loc[ILdf.index[i],'Consensus'] = PBdict[bc]
            # copy Consensus sequence associated with PacBio barcode
            mapped += 1 
            
        else:
            for k,v in PBdict.items(): # loop through PacBio barcodes
                if editdistance.eval(bc,k) <= tolerance: 
                    # compute levenshtein distance between barcodes
                    matches[i] = v
                    #ILdf.loc[ILdf.index[i],'Consensus'] = v 
                    # copy Consensus sequence associated with PacBio barcode if number of edits <= tolerance
                    mapped += 1
                    break
    
    ILdf['Consensus Domains'] = matches
    
    # split consensus domains into individual ICDs
    for i in range(positions):
        ILdf['ICD {}'.format(i+1)] = ILdf['Consensus Domains'].str.split(',',n = positions).str[i]
    
    ILdf.fillna('', inplace = True)
    #print(ILdf)
                
    # overlay ICD metadata
    for i in range(positions):
        ILdf['Fam {}'.format(i+1)] = '' # create empty column for each ICD position
        for j,(k,v) in enumerate(ICDdict.items()): # loop through ICD dictionary
            #meta_positions = ILdf.loc[ILdf['ICD {}'.format(i+1)] == k,['Fam {}'.format(i+1)]]
            #ILdf['Fam {}'.format(i+1)].loc[meta_positions] = v
            ILdf.loc[ILdf['ICD {}'.format(i+1)] == k,['Fam {}'.format(i+1)]] = v
            # assign ICD family metadata based on ICDs listed
                            
    #print(ILdf)
            
    return ILdf, mapped

#%% -- MAIN -- %%#

# import ICD metadata
ICDmeta = pd.read_excel(base + ICDfile, usecols = 'A:B')
ICDdict = {row['ICD']:row['Domain family'] for i,row in ICDmeta.iterrows()}
    
# import PacBio data
PBdf = pd.DataFrame()

for reference in reference_filenames: # combine multiple references from multiple PacBio runs
    print('Importing {} reference'.format(reference))
    PBdata = pd.read_excel(base + 'references/' + reference, 
                           sheet_name = 'DATA', 
                           usecols = 'A:AB')
    PBdf = PBdf.append(PBdata)
    
# filter barcodes
PBdf = PBdf[PBdf['% Freq.'] > 50] # remove anything with a consensus score of less than 50%
PBdf = PBdf[PBdf['Consensus Domains'] != ''] # filter out barcodes that are missing consensus sequences
PBdf['Barcode length'] = PBdf.Barcode.str.len()
PBdf = PBdf[PBdf['Barcode length'] < 24] # remove anything with a barcode of > 24 bp

# count number of occurrences that barcode was found in PacBio data
bc_instances = [0]*len(PBdf.index)
for i,row in PBdf.iterrows():
    if row['Consensus Domains'] != '':
        bc_instances[i] = sum(row != '') - 3

PBdf['Instances'] = bc_instances
    
# keep domain assignments based on sequencing read number and agreement
PBdf['Score'] = PBdf['% Freq.']*PBdf['Instances'] # calculate a consensus score
PBdf.sort_values(by = 'Score', ascending = False) # sort by consensus score
PBdf.drop_duplicates(subset = 'Barcode', keep = 'first', inplace = True) # drop duplicates that were found less often

# convert PacBio barcodes and consensus domains to dictionary    
PBdict = {row['Barcode']:row['Consensus Domains'] for i,row in PBdf.iterrows()}

# create writer object with new file path
for i in range(10000):
    fname = 'dataset_mapped_' + str(i).zfill(4) + '.xlsx'
    if os.path.isdir(os.path.join(base,'processed datasets/',fname)):
        continue
    else:
        break

writer = pd.ExcelWriter(base + 'processed datasets/' + fname, 
                        engine = 'xlsxwriter') # create new Excel spreadsheet

# map Illumina barcodes to PacBio data for each round of selection
idx = 0
for i,(sel,n) in enumerate(zip(settings['selection'],settings['round'])):
    for j,rnd in enumerate(n):
        
        # import Illumina dataset
        ILdf = pd.read_excel(base + 'datasets/CP output/' + dataset_filenames[idx], 
                             names = ['BC','Count'],
                             sheet_name = 'Sheet1') # import Illumina data
        
        print('{} R{}'.format(sel,rnd))
        print('Mapping {} dataset...'.format(settings['sheets'][idx]))
        # sort by count
        data = ILdf.sort_values(by = ['Count'], ascending = False)
        # map barcodes to domains
        ILmapped, mapped = domain_mapping(data,PBdict,ICDfile,tolerance = tolerance)
        ILmapped.to_excel(writer, sheet_name = settings['sheets'][idx]) # export to new Excel sheet
        print('{} dataset mapping finished! {}/{} barcodes mapped.'.format(settings['sheets'][idx],mapped,len(data.index)))
        
        idx += 1

# record settings used to generate processed data
today = datetime.strftime(datetime.now(), "%m/%d/%Y %H:%M:%S")    
    
settings = pd.DataFrame({'dataset_filename':        dataset_filenames,
                        'reference_filenames':      reference_filenames,
                        'ICD list':                 ICDfile,
                        'levenshtein threshold':    tolerance,
                        'date and time':            datetime.strftime(datetime.now(),
                                                    "%m/%d/%Y %H:%M:%S")},
                        orient = 'index')
    
settings.to_excel(writer, sheet_name = 'SETTINGS')

writer.save()

#%% -- test -- %%#

# import and map Illumina data

'''
writer = pd.ExcelWriter('Mapped Illumina test data.xlsx', engine = 'xlsxwriter') # create new Excel spreadsheet

ILdf = pd.read_excel('/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/processed datasets/test_dataset_unnested.xlsx',
                     sheet_name = 'test data')
ILmapped,mapped,gen = DomainSeqMap(ILdf,PBdict,ICDfile,threshold = threshold)
print('Test dataset mapping finished! {}/{} barcodes mapped.'.format(mapped,len(ILmapped.index)))
print('{} are 1st gen, {} are 2nd gen, {} are 3rd gen.'.format(G1,G2,G3))
ILmapped.to_excel(writer, sheet_name = 'test data mapped') # export to new Excel sheet

writer.save()    
writer.close()
'''