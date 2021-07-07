#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:26:17 2021

@author: khloegordon
"""

# -- import standard libraries -- #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import xlsxwriter
import seaborn as sns
#from SeqStat import statistical_analysis

### SETTINGS ###

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'

dataset_filename = '{}/processed datasets/Mapped Illumina data unnested.xlsx'.format(base)
unselected = 'CD69+PD R0'
selections = ['CD69+PD R1','CD69+PD R2','CD69+PD R3','CD69+PD-1- R2']
ICDfile = 'ICD lists and sequences updated 06.24.21.xlsx'

def count2freq(ILdf):
    
    ILdf['Frequency'] = ''
    
    #total = ILdf['Count'].sum(axis = 0)
    ILdf['Frequency'] = ILdf['Count'].divide(ILdf['Count'].sum(axis = 0))
            
    return ILdf

def log2fc(df1,df2,on = 'BC',cols = ['BC','Frequency']):

    '''
    #df1.rename(columns = {'Frequency':'Unselected frequency'})
    df1 = df1.merge(df2['BC'], how = 'inner', on = 'BC', copy = False, 
                    suffixes = (None,'_selected')).sort_values(by = 'BC')
    print(df1)
    '''
    df2 = df2.merge(df1[cols], how = 'inner', on = on,
                    suffixes = ('','_unselected')).sort_values(by = on)
    #print(df2)
    df2['Fold change'] = np.log2(np.divide(df2['Frequency'],df2['Frequency_unselected']))
    
    #df2['Fold change'] = fc
    
    return df2

def feature_filter(df,features,colnames=['ICD 1','ICD 2','ICD 3']):
    '''
    Filters a dataframe for features and returns a list of the number 
    of total occurrences of that feature as well as an array of 
    occurrences of those features within specific columns
        
    Parameters
    ----------    
    df 
        A dataframe containing feature metadata and sequencing counts
    features
        A list of strings containing feature names to be filtered for
    colnames
        A list of column names to look for features in 
        
    Outputs
    ----------
    totals 
        A list of total occurences for each feature in a sequenced population 
        Note: if a feature occurs twice within a given entry, it is only counted once
    positions
        An array of counts for each feature in each column provided in colnames
        Has dimensions (number of ICDs) x (number of positions)
    '''
    
    
    # count instances of a feature in each specified column in a dataframe
    featuredict = {feat:[] for feat in features}
    '''positions = np.zeros((len(features),len(colnames)))'''
    
    for key in featuredict.keys():
        
        for j,col in enumerate(colnames):
            
            featuredict[key].append(df['Count'].loc[df[col] == key].sum(axis = 0))
            '''positions[i,j] = df.loc(df[col]==feature,['Count']).sum(axis = 0)'''
            
        featuredict[key].append(sum(featuredict[key][0:len(colnames)]))
    
    
    # if normalization to all ICD counts at all positions
    total = sum([featuredict[key][-1] for key in featuredict.keys()])
    
    
    '''
    # if normalization of each ICD position to ICD counts at each position
    total = [0]*3
    for j in range(len(colnames)):
        for k,v in featuredict.items():        
            total[j] += v[j]
    '''
    
    '''
    # if normalization of each ICD to each ICD's total counts across all positions
    total = {key:featuredict[key][-1] for key in featuredict.keys()}
    '''
    
    for key in featuredict.keys():
        
        for j in range(len(colnames)):
            
            featuredict[key][j] /= total # if normalization to all ICD counts at all positions
            #featuredict[key][j] /= total[j] # if normalization of each ICD position to ICD counts at each position
            #featuredict[key][j] /= total[key] # if normalization of each ICD to each ICD's total counts across all positions
            '''positions[i,j] = df.loc(df[col]==feature,['Count']).sum(axis = 0)'''
                        
    #print('Feature dictionary\n')
    #print(featuredict)
        
    '''        
    # count all instances of a feature in the dataframe
    totals = sum(np.array(positions),axis = 1)
    
    counts = list(len(ICDs))
    for i, ICD in enumerate(ICDs):
        TF = df.isin(ICD).all(axis = 'columns')
        counts[i] = df[TF]['Count'].sum(axis = 'columns')
        
    # count all instances of a feature in the dataframe
    totals = [df[df[colnames].isin(feature).any(axis = 1)]['Count'].sum(axis = 0) for feature in features]
        
    return positions
    '''
    
    return featuredict

def famfilter(feature_freq,feature_dict):
    
    fams = list(set([val for val in feature_dict.values()]))
        
    famdict = {fam:[] for fam in fams}
    #dict.fromkeys(fams)
    totals = {fam:[] for fam in fams}
                
    for fam in fams:
                
        for feat,freq in feature_freq.items():
            
            if feature_dict[feat] == fam:
                
                famdict[fam].append(freq[-1])
        
        #print(famdict[fam])
                                                                        
        totals[fam] = sum(famdict[fam])
        
    total_feat = sum(totals.values())
    
    for k,v in totals.items():
        
        totals[k] = v/total_feat
            
    #print(totals)
        
    return totals

def CARgen(df,positions = 3):
    
    mapped = df['Count'].loc[df['ICD 1'].notnull()].sum(axis = 0)
    gen = {'Gen {}'.format(i+1):0 for i in range(positions)}
    df['Gen'] = ''
    
    for i in range(positions-1):
        
        df['Gen'].loc[(df['ICD {}'.format(i+1)].notnull()) & (df['ICD {}'.format(i+2)].isnull())] = i+1
        gen['Gen {}'.format(i+1)] = df['Count'].loc[(df['ICD {}'.format(i+1)].notnull()) & 
            (df['ICD {}'.format(i+2)].isnull())].sum(0)/mapped
    
    df['Gen'].loc[df['ICD {}'.format(positions)].notnull()] = positions
    gen['Gen {}'.format(positions)] = df['Count'].loc[df['ICD {}'.format(positions)].notnull()].sum(0)/mapped

    '''
    for i,k in enumerate(gen.keys()):
        print(i)
        if i+1 == positions:
            gen[k] = df['Count'].loc[df['ICD {}'.format(i+1)] != ''].sum(axis = 0)
            df['Gen'].loc[df['ICD {}'.format(i+1)] != ''] = i+1
        else:
            gen[k] = df['Count'].loc[(df['ICD {}'.format(i+1)] != '') & 
               (df['ICD {}'.format(i+2)] == '')].sum(axis = 0)
            print('ICD {} counts: {}'.format(i+1,gen[k]))
            df['Gen'].loc[(df['ICD {}'.format(i+1)] != '') & (df['ICD {}'.format(i+2)] == '')] = i+1
    '''
    
    #print(gensum)
    
    '''
            
    df.loc(df['ICD 1'] != '' and df['ICD 2' == ''],['Gen']) = '1st'
    df.loc(df['ICD 2'] != '' and df['ICD 3' == ''],['Gen']) = '2nd'
    df.loc(df['ICD 3'] != '',['Gen']) = '3rd'
    
    # normalize frequency of each generation of CAR to the proportion that was mapped
    gen = np.divide(gen,mapped)
    
    '''
            
    return df,gen

def heatmap(ICDdf,feat_dict,scaletitle,cols = ['ICD 1','ICD 2','ICD 3'],
            remove = False,vmin = -3, vmax = 3):
    
    ICDs = ICDdf.index.values.tolist()
    ICDdf['Fam'] = ''
    
    for ICD in ICDs:
        
        ICDdf['Fam'] = feat_dict[ICD]
        
    if remove:
        
        ICDdf = ICDdf.loc[~(data==0).all(axis=1)]
    
    ICDdf.replace(np.nan,0) 
    
    ICDdf.sort_values(by = ['Fam'], ascending = [True], inplace = True)
    
    plt.figure(figsize = (4,20))
    plt.rcParams.update({'font.size': 8}) #adjust font size
    ax=sns.heatmap(ICDdf[cols],
                   cmap=plt.cm.get_cmap('Blues'),
                   cbar_kws={'label': scaletitle},
                   vmin=vmin, vmax=vmax) #plot
    ax.tick_params(left=False, bottom=False) #remove tick lines
    
    plt.xlabel('')
    plt.ylabel('')
    
    plt.show()   
    
    return plt

def write2excel(ILdf,ICDdf,Famdf,Gendf,sheet_name):
        
    ILdf.to_excel(writer, sheet_name = sheet_name)
        
    workbook = writer.book
    worksheet = workbook.add_worksheet('{} STATS'.format(sheet_name))
    writer.sheets['{} STATS'.format(sheet_name)] = worksheet
    
    ICDdf.name = 'ICD frequencies'
    Famdf.name = 'ICD family frequencies'
    Gendf.name = 'CAR generation frequencies'

    worksheet.write_string(0,0, ICDdf.name)    
    ICDdf.to_excel(writer,sheet_name = '{} STATS'.format(sheet_name),startrow = 1, 
                   startcol = 0)
    worksheet.write_string(0, ICDdf.shape[1] + 2, Famdf.name)
    Famdf.to_excel(writer,sheet_name = '{} STATS'.format(sheet_name),startrow = 1, 
                   startcol = ICDdf.shape[1] + 2)
    worksheet.write_string(0, ICDdf.shape[1] + Famdf.shape[1] + 4, Gendf.name)
    Gendf.to_excel(writer,sheet_name = '{} STATS'.format(sheet_name),startrow = 1, 
                   startcol = ICDdf.shape[1] + Famdf.shape[1] + 4)
    
    return


### MAIN ###
    
# import ICD metadata

ICDmeta = pd.read_excel(ICDfile, usecols = 'A:B')
ICDdict = {row['ICD']:row['Domain family'] for i,row in ICDmeta.iterrows()}
ICDs = [key for key in ICDdict.keys()]
positions = 3
    
# import and map Illumina data

ILdf0 = pd.read_excel(dataset_filename,sheet_name = unselected)
ILdf0,gen0 = CARgen(ILdf0)
ILdf0 = count2freq(ILdf0)
ICDfreq0 = feature_filter(ILdf0,ICDs)
Famfreq0 = famfilter(ICDfreq0,ICDdict)
# copy stats to dataframes
ICDdf0 = pd.DataFrame.from_dict(ICDfreq0, orient = 'index', 
                               columns  = ['ICD 1','ICD 2','ICD 3','Total'])
Famdf0 = pd.DataFrame.from_dict(Famfreq0, orient = 'index', 
                               columns = ['Frequency'])
Gendf0 = pd.DataFrame.from_dict(gen0, orient = 'index')

for column in list(ILdf0.columns.values):
    if 'Unnamed:' in column: del ILdf0[column]

# check whether Excel file already exists, and if so, make a new one
counter = 1
path = base + 'processed datasets/Mapped Illumina data processed '
while os.path.exists(path + str(counter) + '.xlsx'):
    counter += 1
path = path + str(counter) + '.xlsx'
        
writer = pd.ExcelWriter(path, engine = 'xlsxwriter') # create new Excel spreadsheet
write2excel(ILdf0,ICDdf0,Famdf0,Gendf0,sheet_name = unselected)

'''
ILdf0.to_excel(writer, sheet_name = unselected)

# write stats to STATS worksheet
ICDdf0 = pd.DataFrame.from_dict(feat_freq0, orient = 'index', 
                               columns  = ['ICD 1','ICD 2','ICD 3','Total'])
ICDdf0.name = 'ICD frequencies'
Famdf0 = pd.DataFrame.from_dict(Famfreq0, orient = 'index', 
                               columns = ['Frequency'])
Famdf0.name = 'ICD family frequencies'
Gendf0 = pd.DataFrame.from_dict(gen0, orient = 'index')
Gendf0.name = 'CAR generation frequencies'

workbook = writer.book
worksheet = workbook.add_worksheet('{} STATS'.format(unselected))
writer.sheets['{} STATS'.format(unselected)] = worksheet

worksheet.write_string(0,0, ICDdf0.name)    
ICDdf0.to_excel(writer,sheet_name = '{} STATS'.format(unselected),startrow = 1, 
               startcol = 0)
worksheet.write_string(0, ICDdf0.shape[1] + 2, Famdf0.name)
Famdf0.to_excel(writer,sheet_name = '{} STATS'.format(unselected),startrow = 1, 
               startcol = ICDdf0.shape[1] + 2)
worksheet.write_string(0, ICDdf0.shape[1] + Famdf0.shape[1] + 4, Gendf0.name)
Gendf0.to_excel(writer,sheet_name = '{} STATS'.format(unselected),startrow = 1, 
               startcol = ICDdf0.shape[1] + Famdf0.shape[1] + 4)
'''


for sel in selections:
    
    print('Calculating frequencies and fold changes for {} dataset...'.format(sel))
    ILdf = pd.read_excel(dataset_filename,sheet_name = sel)
    for column in list(ILdf.columns.values):
        if 'Unnamed:' in column: del ILdf[column]
        
    ILdf,gen = CARgen(ILdf) # count and annotate instances of each CAR generation
    ILdf = count2freq(ILdf) # add frequencies
    ILdf = log2fc(ILdf0,ILdf) # calculate fold changes
    ILdf.sort_values(by = 'Fold change',ascending = False) # sort by fold change
    
    ICDfreq = feature_filter(ILdf,ICDs) # count instances of each ICD at each position and across all positions

    Famfreq = famfilter(ICDfreq,ICDdict) # count instances of each ICD family across all positions (last element)
    
    # copy stats to dataframes
    ICDdf = pd.DataFrame.from_dict(ICDfreq, orient = 'index', 
                                   columns  = ['ICD 1','ICD 2','ICD 3','Total'])
    Famdf = pd.DataFrame.from_dict(Famfreq, orient = 'index', 
                                   columns = ['Frequency'])
    Gendf = pd.DataFrame.from_dict(gen, orient = 'index')
    
    # Plot heatmap of log2 frequency of each ICD
    ICDdf_log2freq = np.log2(ICDdf[['ICD 1','ICD 2','ICD 3']])
    heatmap(ICDdf_log2freq,feat_dict = ICDdict,scaletitle = 'log_2(frequency)',
            vmin = -15, vmax = 0)        
    
    # Plot heatmap of log2 fold-change of each ICD
    ICD_log2fc = np.log2(np.divide(ICDdf,ICDdf0))
    heatmap(ICD_log2fc[['ICD 1','ICD 2','ICD 3']],feat_dict = ICDdict,scaletitle = 'log_2(fold-change)')
        
    # write stats in worksheet
    write2excel(ILdf,ICDdf,Famdf,Gendf,sheet_name = sel)
        
# record settings used to generate processed data
settings = pd.DataFrame.from_dict({'dataset_filename':dataset_filename,
                                   'date and time':datetime.strftime(datetime.now(),
                                                   "%m/%d/%Y %H:%M:%S")},
                                    orient = 'index')    
                                   
settings.to_excel(writer, sheet_name = 'SETTINGS')

writer.save()

""" separate data by CAR generation """

'''

for sel in selections:
    
    for n in range(positions):
    
        subset = ILdf.loc['Gen' == n+1]
        
'''    
        




