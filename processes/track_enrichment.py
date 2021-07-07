#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 13:55:47 2021

@author: khloegordon
"""

# -- import standard libraries -- #
import os
import pandas as pd
import matplotlib.pyplot as plt
import xlsxwriter

plt.rcParams['figure.dpi'] = 1200

#%% -- SETTINGS -- %%#

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/processed datasets/'
os.chdir(base)

dataset_filename = '{}dataset_processed_0004.xlsx'.format(base)
selections = ['Unselected','EGFP+ sorted','CD69+ R1','CD69+ R2','CD69+ R3','CD69+PD-1- R1','CD69+PD-1- R2']

# sequences to look for
seqs = {'Consensus Domains':['4-1BB,CD3z'],
         'CD69+ barcodes': ['GTAGGGAACCCGCCGGTG',
                            'CAGAAAATGAGCACCCGT',
                            'GACGGTCGCACACGATTT',
                            'CAACAAAGAGATGTTATC',
                            'CCCGTTGAGGTGATGCAC',
                            'ATGCACACCGCTTACCGT']}
        
# top N barcode frequencies to track enrichment for
top = [10,25,50,100,500]

#%% -- tools -- %%#

def findseq(df,seq,col):
                
    found = df.loc[df[col] == seq]
    
    return found

def findtopN(df,n):
    
    df = df.sort_values(by = 'Frequency', ascending = False)
    topN = df.loc[df.index[:n]]
    topNsum = sum(topN['Frequency'])
    
    return topN,topNsum

#%% -- MAIN -- %%#

tracking = {key:[0]*len(selections) for key in seqs['Consensus Domains'] + seqs['CD69+ barcodes']}
top_enriched = {key:[0]*len(selections) for key in top}

for i,sel in enumerate(selections):
    
    df = pd.read_excel(dataset_filename, sheet_name = sel)
    
    for v in seqs['Consensus Domains']:
        
        found = findseq(df,v,'Consensus Domains')
        tracking[v][i] = sum(found['Frequency'])
        
    '''
    for v in seqs['CD69+PD-1- barcodes']:
        
        found = findseq(df,v,'BC')
        tracking[v][i] = sum(found['Frequency'])
    '''
    
    for v in seqs['CD69+ barcodes']:
        
        found = findseq(df,v,'BC')
        tracking[v][i] = sum(found['Frequency'])
    
    for n in top:
        
        _,topNsum = findtopN(df,n)
        top_enriched[n][i] = topNsum
 
tracked = pd.DataFrame.from_dict(tracking, 
                                 orient = 'index',
                                 columns = selections)

top_df = pd.DataFrame.from_dict(top_enriched,
                                orient = 'index',
                                columns = selections)

# make new Excel file export path
counter = 1
path = base + 'barcode_tracking_'
while os.path.exists(path + str(counter).zfill(4) + '.xlsx'):
    counter += 1

writer = pd.ExcelWriter(path + str(counter).zfill(4) + '.xlsx', engine = 'xlsxwriter')

tracked.to_excel(writer, sheet_name = 'Barcode frequency tracking')
top_df.to_excel(writer, sheet_name = 'Top enriched frequency tracking')

writer.save()
