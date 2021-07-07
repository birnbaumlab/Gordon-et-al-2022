#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 15:45:29 2021

@author: khloegordon
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["figure.figsize"] = (15,8)

#%% -- SETTINGS -- %%#

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'
os.chdir(base)

ICDfile = 'ICD lists and sequences updated 06.24.21.xlsx'
ref1 = 'reference_0009.xlsx'
ref2 = 'reference_0010.xlsx'

dataset_filename = '{}/processed datasets/dataset_processed_0004.xlsx'.format(base)
selections = ['EGFP+ sorted','CD69+ R1','CD69+ R2','CD69+ R3','CD69+PD-1- R1','CD69+PD-1- R2']

#%% -- functions -- %%#

def PB_intersection(PB1,PB2,ref1,ref2,on = 'Barcode'):
    '''calculates the number of unique barcodes mapped in two PacBio runs,
    and the number of barcodes found in both'''
    
    PB1 = PB1.loc[PB1['Consensus Domains'].notnull()] # remove barcodes without consensus domains 
    PB1 = PB1.drop_duplicates(subset = 'Barcode', keep = 'first') # remove duplicates
    PB1_seqs = PB1.drop_duplicates(subset = 'Consensus Domains', keep = 'first')
    
    PB2 = PB2.loc[PB2['Consensus Domains'].notnull()] # remove sequences without consensus domains
    PB2 = PB2.drop_duplicates(subset = 'Barcode', keep = 'first') # remove duplicates
    PB2_seqs = PB2.drop_duplicates(subset = 'Consensus Domains', keep = 'first')
    
    # find common barcodes
    PB_merged = PB2.merge(PB1, how = 'inner', on = on,
                suffixes = ('_2','_1')).sort_values(by = on)
    
    # find common domain combinations
    PB_merged_seqs = PB2_seqs.merge(PB1_seqs, how = 'inner', on = 'Consensus Domains',
                                    suffixes = ('_2','_1')).sort_values(by = 'Consensus Domains')
    
    PB1_matched = len(PB1.index) # count unique barcodes mapped in PB1
    PB2_matched = len(PB2.index) # count unique barcodes mapped in PB2
    both_matched = len(PB_merged) # count unique barcodes mapped in both
    
    PB1_matched_seqs = len(PB1_seqs.index) # count unique domain combinations mapped in PB1
    PB2_matched_seqs = len(PB2_seqs.index) # count unique domain combinations mapped in PB2
    both_matched_seqs = len(PB_merged_seqs) # count unique domain combinations mapped in both

    
    print('{} unique mapped barcodes in {} PacBio data\n{} unique mapped barcodes in {} PacBio data\n{} barcodes in both'.format(
            PB1_matched,ref1,PB2_matched,ref2,both_matched))
    print('{} unique mapped domain combinations in {} PacBio data\n{} mapped domain combinations in {} PacBio data\n{} mapped domain combinations in both'.format(
            PB1_matched_seqs,ref1,PB2_matched_seqs,ref2,both_matched_seqs))
    
    return PB1_matched,PB2_matched,both_matched

def plot_lengths(ILdf,ICDdict,rnd,cols = ['ICD 1','ICD 2','ICD 3'],
                 ILdf0 = False,constant = 307):
    '''
    plot ICD and CAR length distributions unweighted and weighted by count
    '''
    
    fig, axs = plt.subplots(2,2)
    fig.suptitle(rnd) # label with round of selection

    ICDlengths = [] # initialize list of each length of ICD that appears in given round
    weighted_ICDlengths = []
    lengths = [0]*len(ILdf.index) # initialize list of each CAR length in a given round
    weighted_lengths = []
    
    # remove unmapped CARs and CARs with Count = 0
    ILdf = ILdf[ILdf['Consensus Domains'] != '']
    ILdf = ILdf[ILdf['Count'].notna() & ILdf['Count'] != 0]
        
    # loop through Illumina dataset barcodes
    for i,row in ILdf.iterrows():
                
        # loop through CAR ICD positions
        for col in cols:
            # if ICD found in key
            if row[col] in ICDdict.keys():
                
                ICDlengths.append(ICDdict[row[col]]) # add ICD length to ICDlengths
                lengths[i] += ICDdict[row[col]] # add ICD length this barcode's CAR length
                
                weighted_ICDlengths.extend([ICDdict[row[col]] for x in range(int(row['Count']))])
                '''
                for n in range(row['Count']):
                    
                    weighted_ICDlengths.append(ICDdict[row[col]])
                '''
        if lengths[i] != 0: 
            
            lengths[i] += constant
            weighted_lengths.extend([lengths[i] for x in range(int(row['Count']))])        
        '''
        for n in range(row['Count']):
            
            weighted_lengths.append(lengths[i])
        '''
        
    lengths = [i for i in lengths if i != 0]
    weighted_lengths = [i for i in weighted_lengths if i != 0]
                        
    # if unselected Illumina dataset provided, repeat loop for unselected population
    if isinstance(ILdf0, pd.DataFrame):
    
        ICDlengths0 = []
        weighted_ICDlengths0 = []
        lengths0 = [0]*len(ILdf0.index)
        weighted_lengths0 = []
        
        # remove unmapped CARs and CARs with a Count of 0
        ILdf0 = ILdf0[ILdf0['Consensus Domains'] != '']
        ILdf0 = ILdf0[ILdf0['Count'].notna() & ILdf0['Count'] != 0]
            
        for i,row in ILdf0.iterrows():
                    
            # loop through CAR ICD positions
            for col in cols:
                # if ICD found in key
                if row[col] in ICDdict.keys():
                    
                    ICDlengths0.append(ICDdict[row[col]]) # add ICD length to ICDlengths
                    lengths0[i] += ICDdict[row[col]] # add ICD length this barcode's CAR length
                    
                    weighted_ICDlengths0.extend([ICDdict[row[col]] for x in range(int(row['Count']))])
                    
                    '''
                    for n in range(row['Count']):
                        
                        weighted_ICDlengths0.append(ICDdict[row[col]])
                    '''
            
            if lengths0[i] != 0: 
                
                lengths0[i] += constant
                weighted_lengths0.extend(lengths0[i] for x in range(int(row['Count'])))

        lengths0 = [i for i in lengths0 if i != 0]
        weighted_lengths0 = [i for i in weighted_lengths0 if i != 0]
        
        '''
        for n in range(row['Count']):
            
            weighted_lengths0.append(lengths0[i])
        '''
        
        # set bins based on min and max lengths in unselected
        bins1 = np.linspace(min(ICDlengths0), max(ICDlengths0), 10)
        bins2 = np.linspace(min(lengths0), max(lengths0), 10)
        
        # plot ICD length and CAR length distributions
        axs[0,0].hist([ICDlengths0, ICDlengths], bins1, color = ['blue','orange'], 
                 alpha = 1, density = True, label = ['Unselected',rnd])
        axs[0,1].hist([lengths0, lengths], bins2, color = ['blue','orange'], 
             alpha = 1, density = True, label = ['Unselected',rnd])
        axs[1,0].hist([weighted_ICDlengths0, weighted_ICDlengths], bins1, color = ['blue','orange'], 
                 alpha = 1, density = True, label = ['Unselected',rnd])
        axs[1,1].hist([weighted_lengths0, weighted_lengths], bins2, color = ['blue','orange'], 
             alpha = 1, density = True, label = ['Unselected',rnd])
                
    else:
        
        # set bins based on min and max lengths in selected
        bins1 = np.linspace(min(ICDlengths), max(ICDlengths), 10)
        bins2 = np.linspace(min(lengths), max(lengths), 10)
        
        # plot ICD length and CAR length distributions
        axs[0,0].hist(ICDlengths, bins1, color = ['blue','orange'], 
                 alpha = 1, density = True, label = ['Unselected',rnd])
        axs[0,1].hist(lengths, bins2, color = ['blue','orange'], 
             alpha = 1, density = True, label = ['Unselected',rnd])
        plt.legend(loc = 'upper right')
        axs[1,0].hist(weighted_ICDlengths, bins1, color = ['blue','orange'], 
                 alpha = 1, density = True, label = ['Unselected',rnd])
        axs[1,1].hist(weighted_lengths, bins2, color = ['blue','orange'], 
             alpha = 1, density = True, label = ['Unselected',rnd])
    
    axs[0,0].set_xlabel('ICD lengths')
    axs[0,1].set_xlabel('CAR lengths')
    axs[1,0].set_xlabel('Weighted ICD lengths')
    axs[1,1].set_xlabel('Weighted CAR lengths')
    axs[0,0].set_ylabel('Frequency')
    axs[1,0].set_ylabel('Frequency')
    
    axs[0,1].legend(loc = 'upper right')
    
    plt.savefig('figures/'+rnd+' length distribution')
    plt.show()
    
    return lengths,weighted_lengths,lengths0,weighted_lengths0

#%% -- MAIN -- %%#

# import PacBio references
PB1 = pd.read_excel(base + 'references/' + ref1, 
                       sheet_name = 'DATA', 
                       usecols = 'A:D')

PB2 = pd.read_excel(base + 'references/' + ref2, 
                       sheet_name = 'DATA', 
                       usecols = 'A:D')

# import unselected Illumina data
ILdf0 = pd.read_excel(dataset_filename, sheet_name = 'Unselected')

# import ICD lengths
ICDmeta = pd.read_excel(ICDfile, usecols = 'A:C')
ICDs = ICDmeta['ICD'].values.tolist()
lengths = ICDmeta['Length (bp)'].values.tolist()
# make dictionary
ICDdict = {ICD:length for ICD,length in zip(ICDs,lengths)}

# make new Excel file export path
counter = 1
path = base + 'processed datasets/dataset_distributions_'
while os.path.exists(path + str(counter).zfill(4) + '.xlsx'):
    counter += 1
    
writer = pd.ExcelWriter(path + str(counter).zfill(4) + '.xlsx', engine = 'xlsxwriter')

# calculate ICD and CAR length distribution in each round of selection relative to unselected
for sel in selections:
    ILdf = pd.read_excel(dataset_filename, sheet_name = sel)
    lengths,weighted_lengths,lengths0,weighted_lengths0 = plot_lengths(ILdf,ICDdict,sel,ILdf0 = ILdf0)
    distributions = pd.DataFrame(list(zip(lengths,weighted_lengths,lengths0,weighted_lengths0)),
                 columns = ['ICD lengths (bp)','Weighted ICD lengths (bp)','Unselected ICD lengths (bp)','Weighted unselected ICD lengths (bp)'])
    distributions.to_excel(writer, sheet_name = sel)
    
# calculate how many things were uniquely identified in each PacBio run
PB_intersection(PB1,PB2,ref1,ref2) 
