#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 15:32:08 2021

@author: khloegordon
"""

# -- import standard libraries -- #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -- import nonstandard libraries -- #
from scipy.stats import rankdata, pearsonr, binom_test
from random import shuffle
from statistics import mean
from process_datasets import feature_filter

#%% -- settings -- %%#

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'

dataset_filename = '{}/processed datasets/dataset_processed_00011.xlsx'.format(base)
unselected = 'Unselected'
selections = ['CD69+ R1','CD69+ R2','CD69+ R3','CD69+PD-1- R2']
ICDfile = 'ICD lists and sequences updated 06.24.21.xlsx'
positions = 3

#%% -- tools -- #

def FDR(x,q):
    '''
    Applies Benjamini Hochberg multiple hypothesis testing correction
    
    Parameters
    ----------
    x 
        A dataframe containing metadata, counts, and relevant p-values or a
        list or numpy array containing p-values
    q 
        The maximum false discovery rate cut off
        
    Outputs
    ----------
    TF 
        A matrix of True/False values indicating whether null hypothesis was 
        rejected (True) or not (False)
    x 
        The dataframe or list with an appended column indicating whether the null
        hypothesis was rejected after multiple hypothesis testing correction
        was applied
    '''
    
    if isinstance(x, pd.DataFrame):
        x['Significant?'] = ''
        ranks = rankdata(x['p-values']) # assign ranks to p-values
        pstar = ranks/len(ranks)*q # calculate p-value cut off
        TF = [x['p-values'] <= pstar] # create Boolean matrix of indices that meet criteria
        x[TF]['Significant?'] = 'Y' # fill dataframe with Y/N based on whether p-values meet criteria for significance
        x[~TF]['Significant?']= 'N'
    elif isinstance(x,(list,np.ndarray)):
        ranks = rankdata(x) # assign ranks to p-values
        pstar = ranks/len(ranks)*q # calculate p-value cut off
        if isinstance(x,list): 
            x = np.array(x) # convert list type to numpy array
        TF = x <= pstar # create Boolean matrix of indices that meet criteria
        TF_str = ['Y' if val else 'N' for val in TF] # fill list with Y/N based on whether p-values meet criteria for significance
        x = [[p,sig] for p,sig in zip(x,TF_str)] # pack p-values and their corresponding Y/N values into a new list
    else:
        raise ValueError('Invalid input type')
        
    return x,TF

def log2FC(x1,x2):
    '''
    Calculates elementwise log2 fold change between two vectors
    
    Parameters
    ----------
    x1
        A list or array containing counts of each element in sequenced population 1
        If array, dimensions must be (number of constructs or ICDs) x n technical replicates
    x2 are
        A list or array containing counts of each element in sequenced population 2
        (order of must match that in x1)
    '''
    
    if len(x1) != len(x2):
        raise ValueError('x1 and x2 must be of the same length')

    fc = np.log2(np.divide(x2,x1))
    
    return fc

def pearsoncorr(df, order_by = 'Barcode', col1 = 'Frequency', col2 = 'Frequency_unselected'):
    '''
    Determine the Pearson correlation coefficient for count values between two
    populations when pre-filtered for constructs that appear in both
        
    Parameters
    ----------    
    df
        A dataframe containing construct frequencies for population 1 
        (unselected) and population 2 (selected)
    col1, col2
        The column name from which to extract count or frequency values to compute
        the correlation coefficient
        
    Outputs
    ----------
    PCC
        The Pearson correlation coefficient for frequencies of constructs
        in populations 1 and 2
    '''
        
    df1 = df1.sort_values(by = order_by, ascending = 'False')
    df2 = df2.sort_values(by = order_by, ascending = 'False')
    
    PCC = pearsonr(df[col1].values.tolist(),df[col2].values.tolist())
    
    return PCC

#%% -- construct / ICD enrichment studies -- #

def binomial(x1,x2,q=0.05):
    '''
    Calculates the probability of observing each construct x times upon 
    performing selections and Illlumina sequencing, given an initial frequency
    provided by sequencing of the unselected population 

    Parameters
    ----------
    x1 
        A list that contains counts in unselected population
    x2 
        A list that contains counts in selected population 
        Note: order of constructs must match x1
    q
        The false discovery rate cut off for multiple hypothesis testing
        using a Benjamini Hochberg procedure
    '''
    
    if len(x1) != len(x2):
        raise ValueError('x1 and x2 must be of the same length')

    n1,n2 = sum(x1),sum(x2)
    freq = [x/n1 for x in x1]
    #n1,n2 = sum(x1).item(),sum(x2).item()
    #freq = [x.item()/n1 for x in x1]
    #x2 = [x.item() for x in x2]
    pvals = [binom_test(k,n2,p) for k,p in zip(x2,freq)]
    
    # adjust for multiple hypotheses
    if q:
        pvals = FDR(pvals,q)
        
    return pvals

def position(x1,x2,test = 'Wald',positions=3,q=0.05):
    '''
    Calculates the two-tailed p-value for observed occurrences of features in 
    particular structural configurations or positions in a selected population
    
    Parameters
    ----------    
    x1 
        An array containing counts for each feature at each position
        in population 1 (unselected)
        Has dimensions (number of features) x (number of positions) x (number of technical replicates)
        Can be prepared using _,positions = feature_filter(df,features)
    x2 
        An array containing counts for each feature at each position
        in population 1 (unselected)
        Has dimensions (number of features) x (number of positions) x (number of technical replicates)
        Can be prepared using _,positions = feature_filter(df,features)
    test 
        The statistical test to be used to test the null hypotheses that 
        the distributions for each feature at each position the
        distributions between positions are the same
        Can be 'Wald' (default) or 'Fisher' 
    positions
        The number of positions or configuration that a given feature can
        have
    q
        The false discovery rate cut off for multiple hypothesis testing
        using a Benjamini Hochberg procedure
        
    Outputs
    ----------      
    pval1 
        An array giving the two-tailed p-values for the null hypothesis that 
        the underlying distribution of frequencies of a given ICD
        at a given position is the same between the two populations
        Has dimensions (number of ICDs) x (number of positions)
    pval2 
        An array giving the two-tailed p-values for the null hypothesis that
        the underlying distriubtion of frequencies of a given ICD
        at a given position is the same as those at the other positions
        weighted by frequency in unselected population
        has dimensions (number of ICDs) x (number of positions)
        column 1 compares position 1 vs. 2, column 2 compares position 1 vs. 3
        column 3 compares position 2 vs. 3'''
        
    if len(x1) != len(x2):
        raise ValueError('x1 and x2 must be of the same length')
    
    X1,X2 = sum(x1,axis = 2),sum(x2,axis = 2)
    n1,n2 = len(x1.shape[2]),len(x2.shape[2])
    #nm = np.divide(x2,x1)
        
    if test == 'Wald':
        # compare each ICD at each position between populations
        FC = mean(log2FC(x1,x2),axis = 2)
        denom = [[sqrt((2 + n1/n2 + n2/n1)/(col)) for col in row] for row in np.add(X1,X2)]
        T_Wald = np.divide(FC,denom)
        pvals1 = [[norm.pdf(abs(col))*2 for col in row] for row in T_Wald]
        # compare each ICD between positions, normalized for unselected frequency
        idx = 0
        FC_pos = np.zeros()
        for i in range(positions-1):
            for j in range(i+1,positions):
                FC_pos[:,idx] = mean(log2FC(norm[:,i],norm[:,j]))
                sums[:,idx] = np.add(norm_mean[:,i],norm_mean[:,j])
                idx += 1
        '''    
        FC_pos = np.hstack(mean(log2FC(norm[:,0],norm[:,1]),axis = 2),
                           mean(log2FC(norm[:,0],norm[:,2]),axis = 2),
                           mean(log2FC(norm[:,1],norm[:,2]),axis = 2))
        sums = np.hstack(np.add(norm_mean[:,0],norm_mean[:,1]),
                         np.add(norm_mean[:,0],norm_mean[:,2]),
                         np.add(norm_mean[:,1],norm_mean[:,2]))
        '''
        denom_pos = sqrt(np.divide((2 + n1/n2 + n2/n1),sums))
        T_Wald_pos = np.divide(FC_pos,denom_pos)
        pvals2 = [[norm.pdf(abs(col))*2 for col in row] for row in T_Wald_pos]
    
    elif test == 'Fisher':
        # compare each ICD at each position between populations
        m1,m2 = np.sum(x1),np.sum(x2)
        pvals1 = [[fisher_exact([[X1[col],m1-X1[col]],[X2[col],m2-X2[col]]])[1] for col in row] for row in X1]
        # compare each ICD between positions, normalized for unselected frequency
        m_pos = [np.sum(norm[:,i,:]) for i in range(positions)]
        idx = 0
        for i in range(positions-1):
            for j in range(i+1,positions):
                for row in norm:
                    pvals2[row,idx] = fisher_exact([norm[row,i],m[i]-norm[row,i]],[norm[row,j],m[j]-norm[row,j]])
                    idx += 1
        '''
        pvals2[row,0] = fisher_exact([norm[row,0],m11-norm[row,0]],[norm[row,1],m33-norm[row,1]])
        pvals2[row,1] = fisher_exact([norm[row,0],m11-norm[row,0]],[norm[row,2],m33-norm[row,2]])
        pvals2[row,2] = fisher_exact([norm[row,1],m22-norm[row,1]],[norm[row,2],m33-norm[row,2]])
        '''
        
    # adjust for multiple hypotheses    
    if q:
        pvals1 = [FDR(pvals1[:,col],q) for col in range(pvals1.shape[1])]
        pvals2 = [FDR(pvals2[:,col],q) for col in range(pvals2.shape[1])] 
   
    return pvals1, pvals2

#%% -- covariation studies -- #

def SCA(df1,df2,features):
    '''
    Determines the thermodynamic coupling of signaling elements in a 
    combinatorially designed library following functional selections
    
    Parameters
    ----------    
    df1
        A dataframe containing construct and feature metadata as well as
        sequencing counts for population 1 (unselected)
    df2 
        A dataframe containing construct and feature metadata as well as
        sequencing counts for population 2 (selected)
        Note: order must match df1
    features
        A list of strings containing feature names

    Outputs
    ----------    
    ddGij
        An array containing free energy changes reflecting perturbations in 
        frequency for feature j if df2 is filtered by those
        to include only those elements containing feature i
        Has dimensions (number of features) x (number of features) with columns
        corresponding to the perturbation in which feature i is fixed to create
        a subdataframe and rows corresponding to the affected feature
    ax    
        A heatmap plot figure for the ddGij values
        '''
    f1,f2 = feature_filter(df1,features)[0],feature_filter(df2,features)[0]
    freq = np.divide(f1,sum(f1))
    N = len(df2.index)
    Pi = [binom_test(x,N,y) for x,y in zip(f2,freq)]
    Pij = np.zeros(len(features),len(features))
    for i,feature in enumerate(features):
        subdf = df2[df2.isin(feature).any(axis = 1)]
        subcounts = feature_filter(subdf,feature)[0]
        n = len(subdf.index)
        # columns correspond to feature that is fixed, rows reflect binomial probability 
        # of each other feature in response to that perturbation
        Pij[:,i] = [binom_test(x,n,y) for x,y in zip(subcounts,freq)]
        ddGij[:,i] = sqrt(np.subtract(np.log2(Pij[:,i]),np.log2(Pi))**2)
    
    # plot heat map
    fig, ax = plt.subplots()
    im, cbar = heatmap(ddGij,features,features, ax=ax,
                       cmap="viridis", cbarlabel="\delta\deltaG_ij [kT]")
    #texts = annotate_heatmap(im)
    fig.tight_layout()
    plt.show()
    plt.savefig('SCA heatmap',transparent = True)
    
    return ddGij, ax
    
#%% -- functional family studies -- #  

def GSEA(df1,df2,features,sets,p,q):
    '''
    Adapted from Subramanian et al., PNAS 102(43):1545-11550
    Assesses enrichment of a set of features in an enriched population
    
    Parameters
    ----------    
    df1
        A dataframe containing construct and feature metadata as well as
        sequencing counts for population 1 (unselected)
    df2 
        A dataframe containing construct and feature metadata as well as
        sequencing counts for population 2 (selected)
        Note: order must match df1
    features
        A list of strings containing feature names
    sets
        A dictionary containing feature sets as keywords and its corresponding
        features as values
    p
        Weighting factor for correlation with gene set
        When p = 0, enrichment score reduces to a Kolmogorov-Smirnov statistic
        When p = 1, features are weighted by their correlation to the phenotype
        of interest and normalized by the sum of the correlations over all of 
        the features        
    q
        The false discovery rate cut off for multiple hypothesis testing
        using a Benjamini Hochberg procedure

    Outputs
    ----------    
    scores
        An array of scores for each feature set with dimensions
        (number of features + 1) x (number of feature sets)
    ES
        A list of enrichment scores for each feature set
    pvals
        A list of p-values for the observed enrichment scores corresponding to
        each feature set. Computed by randomly permuting the features and their
        corresponding enrichment ratios 1000 times, finding ES_null,
        and calculating the frequency that ES > ES_null
    ax
        A plot figure tracking the running enrichment scores
    '''
    
    # count cumulative normalized frequencies of each feature in enriched population
    feature_counts1 = feature_filter(df1,features)
    feature_counts2 = feature_filter(df2,features)
    # normalize the count frequencies in the enriched population by the
    # count frequencies in the unselected population
    enrichment_ratio = np.divide(np.divide(feature_counts2,sum(feature_counts2)),
                                 np.divide(feature_counts1,sum(feature_counts1)))
    # get indices of features after sorting in descending order
    ranked_idx = np.argsort(enrichment_ratio)[::-1]
    ranked_ratio = enrichment_ratio[ranked_idx]
    ranked_features = features[ranked_idx]
    
    scores = np.zeros(len(ranked_features+1,len(sets.keys)))   
    ES = np.zeros(len(sets.keys))
    for i,key in enumerate(sets.keys): # iterate over unique values in sets
        N_r = np.sum(ranked_ratio[np.in1d(ranked_features,sets[key])])**p
        P_miss = 1/(len(ranked_features)-len(sets[key])) 
        runsum = 0
        for j,f in enumerate(ranked_feature):
            if f in sets[key]:
                P_hit = ranked_ratio[j]^p/N_r
                runsum += P_hit
                scores[j+1,i] = runsum
            else:
                runsum = runsum - P_miss
                scores[j+1,i] = runsum
        if abs(max(scores)) > abs(min(scores)):
            ES[i] = max(scores)
        else:
            ES[i] = min(scores)
        
        # plot running enrichment score for gene set i
        fig, ax = plt.subplots()
        plt.plot(list(range(len(ranked_features+1))),scores[:,i])
        ax.set_xlabel('Feature')
        ax.set_ylabel('Score')
        ax.set_title('Running enrichment score for feature set {} (ES = {})'.format(key,ES[i]))
        fig.tight_layout()
        plt.show()
        plt.savefig('GSEA feature set {}'.format(key),transparent = True)

    # generate null distribution
    ES_null = np.zeros(1000,len(sets.keys))
    for rep in range(1000):
        null_features = shuffle(features)
        scores_null = np.zeros(len(null_features+1,len(sets.keys)))
        for i,key in enumerate(sets.keys): # iterate over unique values in sets
            N_r = np.sum(ranked_ratio[np.in1d(null_features,sets[key])])**p
            P_miss = 1/(len(null_features)-len(sets[key])) 
            runsum = 0
            for j,f in enumerate(null_features):
                if f in sets[key]:
                    P_hit = ranked_ratio[j]^p/N_r
                    runsum += P_hit
                    scores_null[j+1,i] = runsum
                else:
                    runsum = runsum - P_miss
                    scores_null[j+1,i] = runsum
            if abs(max(scores_null)) > abs(min(scores_null)):
                ES_null[rep,i] = max(scores_null)
            else:
                ES_null[rep,i] = min(scores_null)
    
    # normalize enrichment scores to adjust for set size
    u_null = mean(ES_null,axis = 0)
    NES = [ES[i]/u for i,u in enumerate(u_null)]
    for i,u in enumerate(u_null):
        NES_null[:,i] = np.divide(ES_null[:,i],u)
    
    # calculate p-values
    pvals = np.ones(len(sets.keys))
    for i,sc in enumerate(NES):
        if sc > 0:
            pvals[i] = sum(nsc for nsc in NES_null[:,i] if nsc > sc)/1000
        else:
            pvals[i] = sum(nsc for nsc in NES_null[:,i] if nsc > sc)/1000
            
    return scores,ES,pvals,ax
        
#%% -- plotting tools -- #

def volcano(x1,x2,pvals):
    '''
    Creates a volcano plot
    
    Parameters
    ----------
    x1 
        A list containing counts for elements (constructs or features)
        in population 1 (unselected)
    x2 
        A list containing counts for elements (constructs or features)
        in population 2 (selected)
        Note: order must match x1
    p-values
        A list of p-values reflecting the distributions in x1 and x2
    '''
        
    if not ax:
        ax = plt.gca()

    # plot volcano plot
    fig, ax = plt.subplots()
    plt.scatter(log2FC(x1,x2),-np.log10(pvals))
    ax.set_xlabel('log_2FC')
    ax.set_ylabel('-log_10(p-values)')
    fig.tight_layout()
    plt.show()
    plt.savefig('Volcano',transparent = True)

    return ax

#%% -- MAIN -- %%#
    
# make new Excel file export path
counter = 1
path = base + 'processed datasets/dataset_statistical_analysis_'
while os.path.exists(path + str(counter).zfill(4) + '.xlsx'):
    counter += 1

# import ICD lists
ICDdict = pd.read_excel(ICDfile, usecols = 'A:B')

# create new Excel spreadsheet
writer = pd.ExcelWriter(path + str(counter).zfill(4) + '.xlsx', engine = 'xlsxwriter')
    
for sel in selections:
    
    df = pd.read_excel(dataset, sheet_name = sel)
    
    stats = {}
    
    # calculate and record Pearson correlation coefficient
    pcc = pearsoncorr(df)
    stats['Pearson correlation coefficient'] = pcc
    
    # calculate and record binomial pvals for barcode enrichment
    pvals = binomial(df['Frequency_unselected','Frequency'],q = 0.05)
    
    # calculate and record 
    
