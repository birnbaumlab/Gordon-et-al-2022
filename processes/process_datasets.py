#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:26:17 2021

@author: khloegordon
"""

# -- import standard libraries -- #
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
plt.rcParams['figure.dpi'] = 1200

import xlsxwriter
import seaborn as sns
from datetime import datetime
from scipy.stats import pearsonr
from matplotlib import patches

#%% -- SETTINGS -- %%#

'''
base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'

dataset_filename = '{}processed datasets/dataset_mapped_0000_CPcounts.xlsx'.format(base)
unselected = 'EGFP+ sorted'
selections = ['CD69+ R1','CD69+ R2','CD69+ R3','CD69+PD-1- R1','CD69+PD-1- R2']
ICDfile = 'ICD lists and sequences updated 06.24.21.xlsx'
positions = 3 # domain positions in CAR intracellular region
'''

#%% -- tools -- %%#

def count2freq(ILdf):
    '''converts counts to frequencies normalized by all counts in a round of selection'''
    
    ILdf['Frequency'] = ''
    
    #total = ILdf['Count'].sum(axis = 0)
    ILdf['Frequency'] = ILdf['Count'].divide(ILdf['Count'].sum(axis = 0))
            
    return ILdf

def log2fc(df1,df2,on = 'BC',cols = ['BC','Frequency']):
    '''converts frequencies to fold-changes relative to unselected frequencies
    for CARs found in both the selected and unselected populations'''

    og_length = len(df2.index) # number of unique barcodes in selected population
    
    # merge datasets to get intersection of unselected and selected population by barcode
    df2 = df2.merge(df1[cols], how = 'outer', on = on,
                    suffixes = ('','_unselected')).sort_values(by = on)
    
    '''
    filtered_length = len(df2.index) # number of unique barcodes in both populations
    proportion = filtered_length/og_length*100 # fraction of barcodes in selected population that were also in unselected
    print('{} of {} sequences ({}%) found in unselected population'.format(
            filtered_length,og_length,proportion))
    '''
    
    # calculate fold-change
    df2['Fold change'] = np.log2(np.divide(df2['Frequency'],df2['Frequency_unselected']))
    # fold-changes take frequency of each CAR clone relative to all things identified, 
    # not just clones in both populations
        
    return df2

#%% -- filters -- %%#
    
def feature_filter(df,features,colnames=['ICD 1','ICD 2','ICD 3'],norm = True):    
    '''count frequency of a feature (ICD) in each position of the CAR'''
    
    featuredict = {feat:[] for feat in features}
    
    for key in featuredict.keys():
        
        for j,col in enumerate(colnames):
            
            featuredict[key].append(df['Count'].loc[df[col] == key].sum(axis = 0))
            
        featuredict[key].append(sum(featuredict[key][0:len(colnames)]))
    
    if norm:
        
        # normalize ICD counts at a given position to all ICD counts at all positions
        total = sum([featuredict[key][-1] for key in featuredict.keys()])
        
        for key in featuredict.keys():
            
            for j in range(len(colnames)):
                
                featuredict[key][j] /= total
                            
    return featuredict

def famfilter(feature_freq,feature_dict):
    '''count frequency of each family of features (ICDs) using a dictionary 
    dictating {ICD:ICD family}'''
    
    fams = list(set([val for val in feature_dict.values()]))
    famdict = {fam:[] for fam in fams} # collect count of each ICD from a given family, separated by family
    totals = {fam:[] for fam in fams} # collect total count of all ICDs from a given family
                
    for fam in fams:
                
        for feat,freq in feature_freq.items():
            
            if feature_dict[feat] == fam:
                
                famdict[fam].append(freq[-1]) # add total frequency of that ICD in that round
        
        #print(famdict[fam])
                                                                        
        totals[fam] = sum(famdict[fam]) # add all frequencies of all ICDs in a given family
        
    total_feat = sum(totals.values()) # denominator for normalizing ICD family frequencies
    
    for k,v in totals.items():
        
        totals[k] = v/total_feat # normalizing ICD family frequencies
            
    #print(totals)
        
    return totals

def CARgen(df,positions = 3):
    '''assigns each CAR a generation based on number of ICDs in the df
    and counts frequency of each generation of CAR in gen'''
    
    # count number of barcodes with any domains present
    total_counts = df['Count'].sum(axis = 0)
    #mapped = df['Count'].loc[df['ICD 1'].notnull()].sum(axis = 0)
    # initialize dictionary to collect total number of each generation of CAR
    gen = {'Gen {}'.format(i+1):0 for i in range(positions)}
    df['Gen'] = ''
    
    # count frequency of each generation of CAR for non-C-terminal positions
    for i in range(positions-1):
        
        df.loc[(df['ICD {}'.format(i+1)].notnull()) & (df['ICD {}'.format(i+2)].isnull()),['Gen']] = i+1
        #df['Gen'].loc[(df['ICD {}'.format(i+1)].notnull()) & (df['ICD {}'.format(i+2)].isnull())] = i+1
        gen['Gen {}'.format(i+1)] = df['Count'].loc[(df['ICD {}'.format(i+1)].notnull()) & 
            (df['ICD {}'.format(i+2)].isnull())].sum(0)/total_counts        
    
    # count frequency of generation of CAR with the maximum number of positions
    df.loc[df['ICD {}'.format(positions)].notnull(),['Gen']] = positions
    #df['Gen'].loc[df['ICD {}'.format(positions)].notnull()] = positions
    gen['Gen {}'.format(positions)] = df['Count'].loc[df['ICD {}'.format(positions)].notnull()].sum(0)/total_counts
    
    # count frequency of unmapped CARs
    unmapped = 1 - sum([value for value in gen.values()])
    gen['Unmapped'] = unmapped
          
    return df,gen

#%% -- plotting tools -- %%#
    
def heatmap(ICDdf,feat_dict,scaletitle,rnd,saveto,cols = ['ICD 1','ICD 2','ICD 3'],
            remove = False,vmin = -8, vmax = -4,bins = 5):
    ''' plot ICD frequencies or fold-changes at each position
    ICD dataframe contains frequencies at each position and 
    total counts across all positions for each ICD'''
    
    ICDs = ICDdf.index.values.tolist() # extract ICD names
    ICDdf['Fam'] = ''
    
    # assign family to each ICD in the ICD dataframe
    for ICD in ICDs:
        ICDdf.loc[ICD,['Fam']] = feat_dict[ICD]
        
    # remove ICDs that don't appear in this round of selection    
    if remove:
        ICDdf = ICDdf.loc[~(ICDdf == 0).all(axis = 1)]
    
    # replace nan with 0 for plotting
    ICDdf.replace(np.nan,0) 
    
    # sort values by ICD family firstly, and ICD name secondly
    ICDdf = ICDdf.rename_axis('ICD').sort_values(by = ['Fam', 'ICD'], ascending = [True, True])
    
    # plot heatmap
    plt.figure(figsize = (24,3))
    plt.rcParams.update({'font.size': 10}) # adjust font size
    xticks = [i+1 for i in range(len(cols))] # create a column for each ICD position
    
    # create discrete color bar scale
    #cmap = plt.cm.get_cmap('Blues',5)
    cmap = plt.cm.Blues
    #cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = [cmap(i) for i in range(round(cmap.N/8),cmap.N)]
    cmap = colors.LinearSegmentedColormap.from_list('custom cmap',cmaplist,cmap.N)
    
    # define the bins and normalize color bar scale
    bounds = np.linspace(vmin, vmax, bins)
    print(bounds)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # plot heatmap
    ax1 = sns.heatmap(ICDdf[cols].transpose(),
                      cmap = cmap, #plt.cm.get_cmap('Blues'),
                      cbar_kws = {'label': scaletitle},
                      vmin = vmin, vmax = vmax,
                      norm = norm,
                      xticklabels = ICDdf.index.tolist(),
                      yticklabels = xticks) #plot
    
    ax1.tick_params(left = False, bottom = False) # remove tick lines
    ax1.set_title(rnd) # label with selection round
    ax1.figure.axes[-1].yaxis.label.set_size(12)
    
    for _,spine in ax1.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    cbar = ax1.collections[0].colorbar
    cbar.set_ticks(bounds)
    cbar.set_ticklabels(bounds)
                
    plt.xlabel('')
    plt.ylabel('')
    
    plt.tight_layout() # must use with newest version of matplotlib
    
    plt.savefig(saveto+'/'+rnd+' '+scaletitle+'.pdf', transparent = True)
    #plt.show()   
    
    print('{} {}.pdf exported!'.format(rnd,scaletitle))
        
def piechart(ILdf,rnd,saveto,
             colors = {'':'grey',
                       1:'gold',
                       2:'r',
                       3:'b',
                       'bottom':'silver'},
                       top = 30):
    '''plot CAR barcode frequencies in pie chart colored by generation of CAR'''
    
    # sort by CAR barcode frequency
    ILdf = ILdf.sort_values(by = ['Frequency'],ascending = False)
    ILdf = ILdf.dropna(subset = ['Frequency'])
    freqs = ILdf['Frequency'].values.tolist()

    # normalize frequency
    total = sum(freqs)
    
    for i,freq in enumerate(freqs): 
        freqs[i] /= total
            
    # get list of CAR generations corresponding to frequencies
    gen = ILdf['Gen'].values.tolist() 
    
    # create color list based on CAR generation in top most enriched CARs
    color_list = [colors[key] for key in gen[0:top]]
    
    for i in range(top,len(gen)): # color remaining CARs a neutral color
        color_list.append(colors['bottom'])
                
    # plot pie chart
    fig2,ax2 = plt.subplots()
    piechart,_ = ax2.pie(freqs,
                         colors = color_list,
                         wedgeprops = {'edgecolor':'whitesmoke',
                                       'linewidth': 0.5})    
    
    # make outline white
    center = piechart[0].center # get center of piechart
    r = piechart[0].r # get radius of piechart
    circle = patches.Circle(center, r, fill=False,
                                       edgecolor= 'white', linewidth=2) # make circle
    ax2.add_patch(circle) # add circle to piechart
        
    ax2.axis('equal')  # equal aspect ratio ensures that pie is drawn as a circle        
    ax2.set_title(rnd)

    plt.savefig(saveto+'/'+rnd+' '+'pie chart.pdf', transparent = True)
    #plt.show()    
    
    print('{} pie chart exported!'.format(rnd))
    
#%% -- key functions -- %%#
    
def write2excel(writer,ILdf,ICDdf,Famdf,sheet,Gendf = False):
    '''writes all stats for a round of selection to an Excel file'''
        
    ILdf.sort_values(by = 'Frequency', ascending = False)
    ILdf.to_excel(writer, sheet_name = sheet)
        
    workbook = writer.book
    worksheet = workbook.add_worksheet('{} STATS'.format(sheet))
    writer.sheets['{} STATS'.format(sheet)] = worksheet
    
    ICDdf.name = 'ICD frequencies'
    Famdf.name = 'ICD family frequencies'
    
    worksheet.write_string(0,0, ICDdf.name)    
    ICDdf.to_excel(writer,sheet_name = '{} STATS'.format(sheet),startrow = 1, 
                   startcol = 0)
    worksheet.write_string(0, ICDdf.shape[1] + 2, Famdf.name)
    Famdf.to_excel(writer,sheet_name = '{} STATS'.format(sheet),startrow = 1, 
                   startcol = ICDdf.shape[1] + 2)
    
    if isinstance(Gendf, pd.DataFrame):
        
        Gendf.name = 'CAR generation frequencies'
        worksheet.write_string(0, ICDdf.shape[1] + Famdf.shape[1] + 4, Gendf.name)
        Gendf.to_excel(writer,sheet_name = '{} STATS'.format(sheet),startrow = 1, 
                       startcol = ICDdf.shape[1] + Famdf.shape[1] + 4)
    
    return

def process_data(ILdf,ICDs,ICDdict,cols,sel,norm = False, gen = False):
    '''Adds barcode frequencies, log2 fold-changes relative to the unselected 
    population, ICD frequencies at all positions, ICD family frequencies, and
    CAR generation to dataset dataframes
    
    Also creates pie charts and heatmaps for a given round of selection
    '''

    ILdf = count2freq(ILdf) # add barcode frequencies
    if isinstance(norm, pd.DataFrame):
        ILdf = log2fc(norm,ILdf) # calculate fold changes
        ILdf.sort_values(by = 'Fold change',ascending = False) # sort by fold change
        ILdf.fillna({'Frequency':0,'Frequency_unselected':0})
    
    ICDfreq = feature_filter(ILdf,ICDs) # count instances of each ICD at each position and across all positions

    Famfreq = famfilter(ICDfreq,ICDdict) # count instances of each ICD family across all positions (last element)
    
    # copy stats to dataframes
    ICDdf = pd.DataFrame.from_dict(ICDfreq, orient = 'index', 
                                   columns  = cols + ['Total'])
    Famdf = pd.DataFrame.from_dict(Famfreq, orient = 'index', 
                                   columns = ['Frequency'])
    
    if gen: # split by CAR generation
        
        ILdf,gen = CARgen(ILdf)
        Gendf = pd.DataFrame.from_dict(gen, orient = 'index')
    
    else:
        
        Gendf = False
                
    return ILdf,ICDdf,Famdf,Gendf

#%% -- MAIN -- %%#
    
def run(dataset_path,**kwargs):
    
    unselected = ''
    gen_split = True
    freq_plot = True
    FC_plot = False
    pie_plot = False
    positions = 3
    vmin = -8
    vmax = -4
    
    if 'normalize_to' in kwargs:
        unselected = kwargs['normalize_to']

    if 'dataset_sheetnames' in kwargs:
        sheetnames = kwargs['dataset_sheetnames']
                    
    if 'split_by_CAR_generation' in kwargs:
        gen_split = kwargs['split_by_CAR_generation']
        
    if 'domain_reference_filename' in kwargs:
        ICDfile = kwargs['domain_reference_filename']
        
    if 'ICD_positions' in kwargs:
        positions = kwargs['ICD_positions']
        
    if 'frequency_heatmap' in kwargs:
        freq_plot = kwargs['frequency_heatmap']
        
    if 'log2FC_heatmap' in kwargs:
        FC_plot = kwargs['log2FC_heatmap']
        
    if 'pie_plot' in kwargs:
        pie_plot = kwargs['pie_plot']
        
    if 'vmin' in kwargs:
        vmin = kwargs['vmin']
        
    if 'vmax' in kwargs:
        vmax = kwargs['vmax']
        
    # make new Excel file export path
    counter = 1
    path = './processed datasets/processed_dataset_'
    while os.path.exists(path + str(counter).zfill(4) + '.xlsx'):
        counter += 1
    export_path = path + str(counter).zfill(4) + '.xlsx'
            
    # create new Excel spreadsheet
    writer = pd.ExcelWriter(path + str(counter).zfill(4) + '.xlsx', engine = 'xlsxwriter')
        
    # import ICD metadata
    ICDmeta = pd.read_excel('./' + ICDfile, 
                            usecols = 'A:B')
    ICDdict = {row['ICD']:row['Domain family'] for i,row in ICDmeta.iterrows()}
    ICDs = [key for key in ICDdict.keys()]
    #print('ICD metadata imported!')
    
    # set list of columns to use for ICDs in DataFrame
    cols = ['ICD {}'.format(i+1) for i in range(positions)]
    #cols_total = cols + ['Total']
        
    # import unselected Illumina data
    ILdf0 = pd.read_excel(dataset_path,sheet_name = unselected)
    # remove Unnamed columns from DataFrame
    for column in list(ILdf0.columns.values):
        if 'Unnamed:' in column: del ILdf0[column]

    # process unselected Illumina data
    print('Processing unselected data...')
    ILdf0, ICDdf0, Famdf0, Gendf0 = process_data(ILdf0,ICDs,ICDdict,cols,unselected,
                                            norm = False, gen = gen_split)
    print('Unselected data imported')
    
    '''
    ILdf0 = count2freq(ILdf0) # convert counts to frequencies
    ICDfreq0 = feature_filter(ILdf0,ICDs) # count ICD frequencies
    Famfreq0 = famfilter(ICDfreq0,ICDdict) # count ICD family frequencies
    if kwargs['split_by_CAR_generation']:
        ILdf0,gen0 = CARgen(ILdf0) # identify and assign CAR generation
    
    # copy unselected stats to dataframes
    ICDdf0 = pd.DataFrame.from_dict(ICDfreq0, orient = 'index', 
                                   columns  = cols_total)
    Famdf0 = pd.DataFrame.from_dict(Famfreq0, orient = 'index', 
                                   columns = ['Frequency'])
    Gendf0 = pd.DataFrame.from_dict(gen0, orient = 'index')
    '''
        
    # write unselected data to Excel spreadsheet
    write2excel(writer,ILdf0,ICDdf0,Famdf0,sheet = unselected,Gendf = Gendf0)
    
    print('Unselected Illumina data processed and exported!')
    
    # make new directory for figures
    for i in range(10000):
        figure_folder = './figures/statistics_plots_'+ str(i).zfill(4)
        if os.path.isdir(figure_folder):
            continue
        else:
            os.mkdir(figure_folder)
            break

    # generate visuals for unselected Illumina data
    print('Generating visuals for unselected Illumina data')
    if pie_plot:
        # create CAR frequency pie chart colored by generation for unselected population
        piechart(ILdf0,unselected,saveto = figure_folder)

    if freq_plot:
        ICDdf0_log2freq = np.log2(ICDdf0)
        heatmap(ICDdf0_log2freq,
                feat_dict = ICDdict,
                scaletitle = '$log_2$(frequency)',
                saveto = figure_folder,
                vmin = vmin, vmax = vmax,
                rnd = unselected)
        
    # filter out unselected sheet
    sheetnames = [sheet for sheet in sheetnames if unselected not in sheet]
    print(sheetnames)
    
    # loop through each round of selection
    for sheet in sheetnames:
        
        print('Processing {} dataset...'.format(sheet))
        ILdf = pd.read_excel(dataset_path,sheet_name = sheet) # read dataset
        for column in list(ILdf.columns.values): # remove Unnamed columnss
            if 'Unnamed:' in column: del ILdf[column]
        
        '''
        ILdf,gen = CARgen(ILdf) # count and annotate instances of each CAR generation
        ILdf,ICDdf,Famdf,Gendf = process_data(ILdf0,ILdf,gen,ICDs,ICDdict,cols,sheet)
        '''
        if unselected:
            ILdf, ICDdf, Famdf, Gendf = process_data(ILdf,ICDs,ICDdict,cols,sheet,
                                                norm = ILdf0, gen = gen_split)
        else:
            ILdf, ICDdf, Famdf, Gendf = process_data(ILdf,ICDs,ICDdict,cols,sheet,
                                                norm = False, gen = gen_split)

        ILdf = ILdf.sort_values(by = 'Frequency',ascending = False)
                
        # write stats in worksheet
        write2excel(writer,ILdf,ICDdf,Famdf,sheet = sheet,Gendf = Gendf)   
        
        # generate visuals for unselected Illumina data
        print('Generating visuals for {} Illumina data'.format(sheet))
        if pie_plot:
            # create CAR frequency pie chart colored by generation for unselected population
            piechart(ILdf,sheet,saveto = figure_folder)

        if freq_plot:
            ICDdf_log2freq = np.log2(ICDdf)
            heatmap(ICDdf_log2freq,
                    feat_dict = ICDdict,
                    cols = cols,
                    scaletitle = '$log_2$(frequency)',
                    saveto = figure_folder,
                    vmin = vmin, vmax = vmax,
                    rnd = sheet)
            
        if FC_plot:
            ICDdf_log2fc = np.log2(np.divide(ICDdf,ICDdf0))
            heatmap(ICDdf_log2fc, 
                    feat_dict = ICDdict, 
                    cols = cols,
                    scaletitle = '$log_2$(fold-change)',
                    saveto = figure_folder,
                    rnd = sheet)

        print('{} Illumina data processed and figures exported!'.format(sheet))
            
    """ record settings used to generate processed data """
    today = datetime.strftime(datetime.now(),"%m/%d/%Y %H:%M:%S")
    settings = {'dataset_filename':  [dataset_path],
                'normalized to':     [unselected],
                'date and time':     [today]}
    
    settings = pd.DataFrame(list(settings.values()),index = list(settings.keys()),columns=['Input'])
                                       
    settings.to_excel(writer, sheet_name = 'SETTINGS')
    
    writer.save()
        
    return export_path, figure_folder

'''
### Run script ###
        
base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'
os.chdir(base)

### --- analysis and visuals settings --- ###
settings = {# mapped dataset filename, if already processed -- overrides domain_mapping
            # must be placed in processed datasets folder
            'mapped_dataset_filename':   
                                         'mapped_dataset_CD4_sel1.xlsx',
            'dataset_sheetnames':                       ['1A unstim d0',
                                                         '1B unstim d1',
                                                          '1C stim d12',
                                                        '1D unstim d13',
                                                          '1E stim d13',
                                                            '1F CM d13',
                                                           '1G SCM d13',
                                                            '1H EM d13',
                                                           '1I EFF d13'
                                                                      ],
            'domain_reference_filename': 
                        'ICD lists and sequences updated 12.21.21.xlsx',
            'ICD_positions':                                          3,
            'normalize_to':                              '1B unstim d1',
            'frequency_heatmap':                                   True,
            'log2FC_heatmap':                                     False,
            'pie_plot':                                           False,
            # (generates pie charts to show CAR generation breakdown, 
               #adds CAR generation column to excel sheet)
            'split_by_CAR_generation':                             True,
            }
        
# analyze ICD frequencies and generate plots
print('Analyzing dataset...')
processed_path = run(base + 'processed datasets/' + settings['mapped_dataset_filename'],**settings)
print('Finished! Exported to ' + processed_path[0])
'''


### Archived code ###
        
'''
# filter each dataset for each round of selection for given CAR generation 
for i in range(positions):
    
    pos_subset = list(range(0,i+1)) # ICD columns to keep for i+1_th generation
    cols_subset = ['ICD {}'.format(j+1) for j in range(i+1)]
    
    # create ICD frequencies for i_th generation of CARs for unselected population
    ILdf0_subset = ILdf0[ILdf0['Gen'] == i+1] # filter for i+1_th generation
    ICDfreq0_subset = feature_filter(ILdf0_subset,ICDs) # count ICD frequencies
    ICDdf0_subset = pd.DataFrame.from_dict(ICDfreq0_subset, orient = 'index')
    # keep ICD frequencies for positions included in i+1_th generation CARs only
    ICDdf0_subset = ICDdf0_subset.iloc[:,pos_subset] 
    ICDdf0_subset.columns = cols_subset

    # create ICD frequencies for i+1_1_th generation of CARs for selected population
    ILdf_subset = ILdf[ILdf['Gen'] == i+1]
    ICDfreq_subset = feature_filter(ILdf_subset,ICDs)
    ICDdf_subset = pd.DataFrame.from_dict(ICDfreq_subset, orient = 'index')
    ICDdf_subset = ICDdf_subset.iloc[:,pos_subset]
    ICDdf_subset.columns = cols_subset
    
    # Plot heatmap of log2 frequency of each ICD for i+1_th generation of CARs
    ICDdf_subset_log2freq = np.log2(ICDdf_subset)
    heatmap(ICDdf_subset_log2freq, feat_dict = ICDdict, cols = cols_subset,
            scaletitle = '$log_2$(frequency)',saveto = figure_folder,
            vmin = -12, vmax = 0,bins = 7,rnd = sel + ' Gen {} CARs'.format(i+1))        
    
    # Plot heatmap of log2 fold-change of each ICD for i+1_th generation of CARs
    ICD_subset_log2fc = np.log2(np.divide(ICDdf_subset,ICDdf0_subset))
    heatmap(ICD_subset_log2fc, feat_dict = ICDdict, cols = cols_subset,
            scaletitle = '$log_2$(fold-change)',saveto = figure_folder,
            rnd = sel + ' Gen {} CARs'.format(i+1))

# filter for 2nd and 3rd gen CARs combineder
ILdf0_subset = ILdf0[ILdf0['Gen'] != 1]
ICDfreq0_subset = feature_filter(ILdf0_subset,ICDs)
ICDdf0_subset = pd.DataFrame.from_dict(ICDfreq0_subset, orient = 'index', 
                           columns  = cols_subset + ['Total'])

ILdf_subset = ILdf[ILdf['Gen'] != 1]
ICDfreq_subset = feature_filter(ILdf_subset,ICDs)
ICDdf_subset = pd.DataFrame.from_dict(ICDfreq_subset, orient = 'index', 
                           columns  = cols_subset + ['Total'])

# Plot heatmap of log2 frequency of each ICD
ICDdf_subset_log2freq = np.log2(ICDdf_subset[cols_subset])
heatmap(ICDdf_subset_log2freq[cols_subset],feat_dict = ICDdict,
        scaletitle = '$log_2$(frequency)',saveto = figure_folder,
        vmin = -20, vmax = 0,rnd = sel + ' Gen 2 + 3 CARs'.format(i+1))        

# Plot heatmap of log2 fold-change of each ICD
ICD_subset_log2fc = np.log2(np.divide(ICDdf_subset,ICDdf0_subset))
heatmap(ICD_subset_log2fc[cols_subset],feat_dict = ICDdict,scaletitle = '$log_2$(fold-change)',
        rnd = sel + ' Gen 2 + 3 CARs',saveto = figure_folder)

def plot(export_path):
    
    # make new directory for figures
    for i in range(10000):
        figure_folder = os.path.join(base,'figures/heatmaps_piecharts_',str(i).zfill(4))
        if os.path.isdir(figure_folder):
            continue
        else:
            os.mkdir(figure_folder)
            break
    
    # create CAR frequency pie chart colored by generation for unselected population
    piechart(ILdf0,unselected,saveto = figure_folder)
    
    # Plot pie chart of clonal CAR frequency colored by CAR generation
    piechart(ILdf,sel,saveto = figure_folder)
    
    # Plot heatmap of log2 frequency of each ICD
    ICDdf_log2freq = np.log2(ICDdf[cols])
    heatmap(ICDdf_log2freq,feat_dict = ICDdict,scaletitle = '$log_2$(frequency)',
            vmin = -20, vmax = 0,rnd = sel,saveto = figure_folder)     
    
    # Plot heatmap of log2 fold-change of each ICD
    ICD_log2fc = np.log2(np.divide(ICDdf,ICDdf0))
    heatmap(ICD_log2fc[cols],feat_dict = ICDdict,scaletitle = '$log_2$(fold-change)',
            rnd = sel,saveto = figure_folder)

'''
