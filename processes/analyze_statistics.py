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
import seaborn as sns
import os

# -- import nonstandard libraries -- #
from scipy.stats import rankdata, fisher_exact, norm, pearsonr, binom_test
from scipy.special import comb
from random import shuffle
from statistics import mean
from matplotlib.collections import PatchCollection
from matplotlib import colors
from math import sqrt

#%% -- settings -- %%#

base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'
os.chdir(base)
dataset = '{}/processed datasets/dataset_processed_0005.xlsx'.format(base)
unselected = 'EGFP+ sorted'
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

def binom_prob(f1,N,x2):
        
    prob = comb(N,x2)*f1**x2*(1-f1)**(N-x2)
    
    return prob


#%% -- construct / ICD enrichment studies -- #

def binomial(f1,x2,q=0.05):
    '''
    Calculates the probability of observing each construct x times upon 
    performing selections and Illlumina sequencing, given an initial frequency
    provided by sequencing of the unselected population 

    Parameters
    ----------
    f1 
        A list or array that contains frequency in unselected population
    x2 
        A list or array that contains counts in selected population 
        Note: order of constructs must match x1
    q
        The false discovery rate cut off for multiple hypothesis testing
        using a Benjamini Hochberg procedure
    '''
    
    n2 = sum(x2)
    pvals = [binom_test(k,n2,p,alternative = 'greater') for k,p in zip(x2,f1)]
    
    # adjust for multiple hypotheses
    if q:
        pvals = FDR(pvals,q)
        
    return pvals

def Waldlog(x1,x2,q=0.05):
    '''
    Calculates the two-tailed p-value for observed read counts in two 
    different populations. Assumes a negative binomial distribution describing
    the technical variation in the sequencing count data. The null hypothesis
    is that the means for the counts of the two populations are the same. 
    Adapted from the DESeq2 manuscript and Chen et al. BMC Systems Biology, 5 (2011)
    
    Parameters
    ----------  
    x1 
        A list or array with the counts for each construct or feature in 
        the unselected population
        Has dimensions (number of elements) rows x (n1 replicates) columns
        Elements can be constructs or features
    x2 
        A list or array with the counts for each construct or feature in 
        the selected population
        Has dimensions (number of elements) rows x (n2 replicates) columns
        Elements can be constructs or features
        Note: order must match x1
    q
        The false discovery rate cut off for multiple hypothesis testing
        using a Benjamini Hochberg procedure

    Outputs
    ----------
    pvals
        An array of p-values for each element in x2 given x1
    '''
    
    if len(x1) != len(x2):
        raise ValueError('x1 and x2 must be of the same length')

    if isinstance(x1[0],list):
        X1 = np.sum(x1,axis = 1)
        n1 = len(x1[0])
    else:
        X1 = x1
        n1 = 1 # 1 technical replicate

        
    if isinstance(x2[0],list):
        X2 = np.sum(x2,axis = 1)
        n2 = len(x2[0])
    else:
        X2 = x2
        n2 = 1 # 1 technical replicate
    
    FC = log2FC(X1,X2)
    
    '''            
    FC = np.mean(log2FC(x1,x2),axis = 1) 
    # sum counts across technical replicates for each element
    X1,X2 = sum(x1,axis = 1),sum(x2,axis = 1) 
    n1,n2 = len(x1[1]),len(x2[1])
    '''

    # calculate denominator for Wald statistic
    denom = [sqrt((2 + n1/n2 + n2/n1)/(a + b)) if (a != 0 and b != 0) else None for a,b in zip(X1,X2)] 
    T_Wald = [fc/d if d else None for fc,d in zip(FC,denom)] # calculate Wald statistic
    # compare Wald statistic to normal distribution, as is done for DESeq2
    pvals = [norm.pdf(abs(ts))*2 if ts else 1 for ts in T_Wald] 
    #pvals = [nbinom.ppf(abs(ts))*2 for ts in T_Wald]

    # adjust for multiple hypotheses    
    if q:
        pvals = FDR(pvals,q)
    
    return pvals
        
    '''
    if mode == 'construct':
    if mode == 'ICD':
        if ICDs:
            X1 = list()
            for i in range(n1):
                filtered = ICD_filter(x1[:,i],ICDs)
                np.hstack(X1,filtered)
            X2 = ICD_filter(x2,ICDs)
            FC = mean(log2FC(X1,X2),axis = 1)
        else:
            raise ValueError('Must provide ICD list')
    else:
        raise ValueError('Must choose construct or ICD enrichment mode')'''

'''
def position(x1,x2,test = 'Wald',positions=3,q=0.05):
    
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
        column 3 compares position 2 vs. 3
        
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
        sums = np.zeros()
        for i in range(positions-1):
            for j in range(i+1,positions):
                FC_pos[:,idx] = mean(log2FC(norm[:,i],norm[:,j]))
                sums[:,idx] = np.add(norm_mean[:,i],norm_mean[:,j])
                idx += 1
          
        FC_pos = np.hstack(mean(log2FC(norm[:,0],norm[:,1]),axis = 2),
                           mean(log2FC(norm[:,0],norm[:,2]),axis = 2),
                           mean(log2FC(norm[:,1],norm[:,2]),axis = 2))
        sums = np.hstack(np.add(norm_mean[:,0],norm_mean[:,1]),
                         np.add(norm_mean[:,0],norm_mean[:,2]),
                         np.add(norm_mean[:,1],norm_mean[:,2]))
        
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
                    pvals2[row,idx] = fisher_exact([norm[row,i],m_pos[i]-norm[row,i]],[norm[row,j],m_pos[j]-norm[row,j]])
                    idx += 1
        
        pvals2[row,0] = fisher_exact([norm[row,0],m11-norm[row,0]],[norm[row,1],m33-norm[row,1]])
        pvals2[row,1] = fisher_exact([norm[row,0],m11-norm[row,0]],[norm[row,2],m33-norm[row,2]])
        pvals2[row,2] = fisher_exact([norm[row,1],m22-norm[row,1]],[norm[row,2],m33-norm[row,2]])
        
        
    elif test == 'binomial':
        
        for i in range(positions-1):
            for j in range(i+1,positions):
                for row in norm:
                    pvals2[row,idx] = fisher_exact([norm[row,i],m_pos[i]-norm[row,i]],[norm[row,j],m_pos[j]-norm[row,j]])
                    idx += 1

        
        
    # adjust for multiple hypotheses    
    if q:
        pvals1 = [FDR(pvals1[:,col],q) for col in range(pvals1.shape[1])]
        pvals2 = [FDR(pvals2[:,col],q) for col in range(pvals2.shape[1])] 
   
    return pvals1, pvals2

'''

#%% -- covariation studies -- #

def SCA(df1,df2,features,rnd,saveto):
    '''
    Determines the thermodynamic coupling of signaling elements in a 
    combinatorially designed library following functional selections
    Adapted from Lockless, S.W. and Ranganathan, R. Science 286. 1999.
    
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
    # extract total counts for each ICD in each population
    x1,x2 = feature_filter(df1,features),feature_filter(df2,features)
    x1 = [x[-1] for x in list(x1.values())]
    x2 = [x[-1] for x in list(x2.values())]
    
    # normalize frequency in unselected population
    f1 = [x/sum(x1) for x in x1]
    
    N = sum(x2)
    Pi = [binom_prob(x,N,y) for x,y in zip(x2,f1)]
    Pij = np.zeros((len(features),len(features)))
    ddGij = Pij
    
    for i,feature in enumerate(features):
        subdf = df2[df2.isin([feature]).any(axis = 1)]
        subcounts = feature_filter(subdf,features)
        subcounts = [x[-1] for x in list(subcounts.values())]
        n = sum(subcounts)
        # columns correspond to feature that is fixed, rows reflect binomial probability 
        # of each other feature in response to that perturbation
        pij = [binom_prob(x,n,y) for x,y in zip(subcounts,f1)]
        Pij[:,i] = pij[:]
        M = (np.subtract(np.log2(Pi),np.log2(pij)))
        ddgij = [-1/N*m for m in M]
        ddGij[:,i] = ddgij[:]
    
    ddGij[ddGij < -1e308] = 0 # remove infinite values
    
    # plot heat map
    '''
    fig, ax = plt.subplots()
    im, cbar = heatmap(data = ddGij,col_labels = features,row_labels = features, 
                       xlabel = 'ICDs',ylabel = 'ICDs',ax=ax,
                       title = 'Statistical coupling analysis of ICDs',
                       cmap="viridis", cbarlabel=r'$\Delta\Delta G_{ij}$ [kT]')
    plt.xticks(fontsize=10, rotation=45)
    plt.yticks(fontsize=10, rotation=0)
    
    #texts = annotate_heatmap(im)
    #fig.tight_layout()
    '''
    
    cmap = plt.cm.Blues
    cmaplist = [cmap(i) for i in range(cmap.N)]
    #cmaplist = [cmap(i) for i in range(round(cmap.N/8),cmap.N)]
    cmap = colors.LinearSegmentedColormap.from_list('custom cmap',cmaplist,cmap.N)
    
    # define the bins and normalize color bar scale
    vmin = 0
    vmax = 1
    bins = 5
    
    bounds = np.linspace(vmin, vmax, bins)
    print(bounds)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    # plot heatmap
    plt.figure(figsize = (24,28))
    plt.rcParams.update({'font.size': 12}) # adjust font size

    ax = sns.heatmap(ddGij,
                      cmap = cmap, #plt.cm.get_cmap('Blues'),
                      cbar_kws = {'label': r'$\Delta\Delta G_{ij}$ [kT]'},
                      vmin = vmin, vmax = vmax,
                      norm = norm,
                      xticklabels = features,
                      yticklabels = features) #plot
    
    ax.tick_params(left = False, bottom = False) # remove tick lines
    ax.set_title('Statistical coupling analysis of ICDs ' + rnd) # label with selection round
    
    for _,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks(bounds)
    cbar.set_ticklabels(bounds)

                    
    plt.tight_layout() # must use with newest version of matplotlib
    plt.savefig(saveto + '/SCA heatmap ' + rnd + '.pdf',transparent = True)
    plt.show()
    
    return ddGij, Pij, Pi, ax
    
#%% -- functional family studies -- #  

def GSEA(rnd,f1,f2,features,sets,p,saveto,q = 0.05):
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
    
    # normalize the count frequencies in the enriched population by the
    # count frequencies in the unselected population
    enrichment_ratio = np.divide(f2,f1)
    # get indices of features after sorting in descending order
    ranked_idx = np.argsort(enrichment_ratio)[::-1]
    ranked_ratio = enrichment_ratio[ranked_idx]
    ranked_features = np.array(features)[ranked_idx.astype(int)]
    
    n = len(ranked_features)
    m = len(sets.keys())
    
    scores = np.zeros((n+1,m))
    ES = np.zeros(m)
    for i,(key,value) in enumerate(sets.items()): # iterate over unique values in sets
        for j,bc in enumerate(ranked_ratio):
            N_r = np.sum(ranked_ratio[np.in1d(ranked_features,value)])**p
            P_miss = 1/(n-len(value)) 
            runsum = 0
            for j,f in enumerate(ranked_features):
                if f in value:
                    P_hit = ranked_ratio[j]**p/N_r
                    runsum += P_hit
                else:
                    runsum = runsum - P_miss
                scores[j+1,i] = runsum
            if abs(np.amax(scores)) > abs(np.amin(scores)):
                ES[i] = np.amax(scores)
            else:
                ES[i] = np.amin(scores)
        
    # generate null distribution
    ES_null = np.zeros((1000,m))
    for rep in range(1000):
        null_features = shuffle(features)
        scores_null = np.zeros((n+1,m))
        for i,(key,value) in enumerate(sets.items()): # iterate over unique values in sets
            N_r = np.sum(ranked_ratio[np.in1d(null_features,value)])**p
            P_miss = 1/(n-len(value)) 
            runsum = 0
            for j,f in enumerate(null_features):
                if f in value:
                    P_hit = ranked_ratio[j]**p/N_r
                    runsum += P_hit
                else:
                    runsum = runsum - P_miss
                scores_null[j+1,i] = runsum
            if abs(np.amax(scores_null)) > abs(np.amin(scores_null)):
                ES_null[rep,i] = np.amax(scores_null)
            else:
                ES_null[rep,i] = np.amin(scores_null)
    
    # normalize enrichment scores to adjust for set size
    u_null = mean(ES_null,axis = 0)
    NES = [ES[i]/u for i,u in enumerate(u_null)]
    NES_null = np.zeros((1000,m))
    for i,u in enumerate(u_null):
        NES_null[:,i] = np.divide(ES_null[:,i],u)
    
    # calculate p-values
    pvals = np.ones(m)
    for i,sc in enumerate(NES):
        if sc > 0:
            pvals[i] = sum(nsc for nsc in NES_null[:,i] if nsc > sc)/1000
        else:
            pvals[i] = sum(nsc for nsc in NES_null[:,i] if nsc < sc)/1000
            
    # adjust for multiple hypotheses    
    if q:
        pvals = FDR(pvals,q)        
            
    # plot running enrichment score for gene set i
    plt.rcParams["figure.figsize"] = (15,8)
    fig, axs = plt.subplots(2,2)
    fig.suptitle(rnd) # label with round of selection
    
    fams = list(sets.keys())
    
    for i in range(2):
        for j in range(2):
            axs[i,j].plot(scores[:,i+j])
            axs[i,j].set_xlabel('ICDs')
            axs[i,j].set_ylabel('Running enrichment score')
            axs[i,j].set_title(fams[i+j])
                    
    plt.savefig('figures/'+rnd+' GSEA')
    plt.show()

    plt.plot(list(range(n+1)),scores[:,i])
    axs.set_xlabel('Feature')
    axs.set_ylabel('Score')
    axs.set_title('Running enrichment score for feature set {} (ES = {})'.format(key,ES[i]))
    fig.tight_layout()
    plt.savefig(saveto+'/GSEA {}'.format(rnd),transparent = True)
    plt.show()
    
    return pvals,scores,NES,axs

        
#%% -- plotting tools -- #

def volcano(data,title,rnd,saveto,ax = []):
    '''
    Creates a volcano plot
    
    Parameters
    ----------
    data
        A list of lists containing [x1,x2,p-values]
    x1 
        A list containing counts for elements (constructs or features)
        in population 1 (unselected)
    x2 
        A list containing counts for elements (constructs or features)
        in population 2 (selected)
        Note: order must match x1
    p-values
        A list of p-values reflecting the distributions in x1 and x2
    title
        A list of strings containing subplot titles
    '''
        
    if not ax:
        ax = plt.gca()

    # plot volcano plot
    if len(data) == 4:
        fig, ax = plt.subplots(2,2)
        plt.rcParams["figure.figsize"] = (15,8)
        for i,d in enumerate(data):
            plt.scatter(log2FC(d[0],d[1]),-np.log10(d[3]))
            ax.set_xlabel('log_2FC')
            ax.set_ylabel('-log_10(p-values)')
            ax.set_title(titles[i])
    if len(data) == 2:
        fig, ax = plt.subplots(1,2)
        plt.rcParams["figure.figsize"] = (15,6)
        for i,d in enumerate(data):
            plt.scatter(log2FC(d[0],d[1]),-np.log10(d[3]))
            ax.set_xlabel('log_2FC')
            ax.set_ylabel('-log_10(p-values)')
            ax.set_title(titles[i])
    else:
        fig,ax = plt.subplot(1,len(data))
        plt.rcParams["figure.figsize"] = (8*len(data),6)
        for i,d in enumerate(data):
            plt.scatter(log2FC(d[0],d[1]),-np.log10(d[3]))
            ax.set_xlabel('log_2FC')
            ax.set_ylabel('-log_10(p-values)')
            ax.set_title(titles[i])
        
            
    fig.suptitle(rnd)      
    fig.tight_layout()
    plt.show()
    plt.savefig(rnd + title,transparent = True)

    return ax

def circlemap(ICDdf,feat_dict,pvals,scaletitle,rnd,saveto,pvalmethod,cols = ['ICD 1','ICD 2','ICD 3'],
            remove = False,vmin = 0, vmax = 4,bins = 5):
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
        
    x, y = np.meshgrid(np.arange(len(ICDs)), np.arange(len(cols)))
    fig, ax = plt.subplots()

    R = pvals/pvals.max()/2
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=ICDdf[cols].transpose().flatten(), cmap=cmap,
                          cbar_kws = {'label': scaletitle},
                          vmin = vmin, vmax = vmax,
                          norm = norm)
    ax.add_collection(col)
    
    ax.set(xticks=np.arange(len(ICDs)), yticks=np.arange(len(cols)),
       xticklabels=ICDs, yticklabels=cols)
    ax.set_xticks(np.arange(len(ICDs)+1)-0.5, minor=True)
    ax.set_yticks(np.arange(len(cols)+1)-0.5, minor=True)
    ax.grid(which='minor')
    ax.tick_params(left = False, bottom = False) # remove tick lines
    ax.set_title(rnd) # label with selection round
    ax.figure.axes[-1].yaxis.label.set_size(12)
    
    for _,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    cbar = ax.collections[0].colorbar
    cbar.set_ticks(bounds)
    cbar.set_ticklabels(bounds)
                
    fig.colorbar(col)
    
    plt.tight_layout() # must use with newest version of matplotlib
    
    plt.savefig(saveto+'/'+rnd+' '+scaletitle+pvalmethod+'circle map.pdf', transparent = True)
    plt.show()   
    
    print('{} {}.pdf exported!'.format(rnd,scaletitle))
    
def heatmap(data,xlabel,ylabel,title,col_labels,row_labels,
            cmap,cbar_kw={},cbarlabel="",ax = [],**kwargs):
    '''
    Copied from https://www.google.com/search?q=import+heatmap+matplotlib&rlz=1C5CHFA_enUS854US854&oq=import+heatmap+matplotlib&aqs=chrome..69i57j0.3439j0j7&sourceid=chrome&ie=UTF-8
    
    Creates a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    '''

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar
    
#%% -- MAIN -- %%#
    
# make new Excel file export path
counter = 1
path = base + 'processed datasets/dataset_statistical_analysis_'
while os.path.exists(path + str(counter).zfill(4) + '.xlsx'):
    counter += 1

# import ICD lists
ICDmeta = pd.read_excel(ICDfile, usecols = 'A:B')
ICDdict = {row['ICD']:row['Domain family'] for i,row in ICDmeta.iterrows()}
ICDs = [key for key in ICDdict.keys()]

# create ICD dictionary with families as keywords and ICDs as values
fams = list(set([val for val in ICDdict.values()]))
famdict = {fam:[] for fam in fams} # collect count of each ICD from a given family, separated by family
for k,v in ICDdict.items():
    famdict[v].append(k)

# create new Excel spreadsheet
writer = pd.ExcelWriter(path + str(counter).zfill(4) + '.xlsx', engine = 'xlsxwriter')

# import unselected data
df0 = pd.read_excel(dataset, sheet_name = unselected)
df0 = df0.sort_values(by = 'BC')
ICD_counts0 = feature_filter(df0,ICDs)
ICDdf0 = pd.DataFrame.from_dict(ICD_counts0, orient = 'index', 
                               columns  = ['ICD 1','ICD 2','ICD 3','Total'])
ICDdf_totals0 = [x/sum(ICDdf0['Total']) for x in list(ICDdf0['Total'].values)]

# make new directory for figures
for i in range(10000):
    figure_folder = base + 'figures/statistics_plots_' + str(i).zfill(4)
    if os.path.isdir(figure_folder):
        continue
    else:
        os.mkdir(figure_folder)
        break

pcc = {}

for sel in selections:
    
    df = pd.read_excel(dataset, sheet_name = sel)
    df_filtered = df[df['Fold change'].notnull()]
    df_filtered = df_filtered.sort_values(by = 'BC')
    df0_merged = df0.merge(df_filtered, how = 'inner', on = 'BC', suffixes = ('_unselected','_selected'))
    df0_merged = df0_merged.sort_values(by = 'BC')
    x0 = df0_merged['Count_unselected']
    f0 = df_filtered['Frequency_unselected']
    x = df_filtered['Count']
    f = df_filtered['Frequency']
    
    # count frequency of each ICD at each position and total ICD count across all positions
    ICD_counts = feature_filter(df,ICDs)
    ICDdf = pd.DataFrame.from_dict(ICD_counts, orient = 'index', 
                                   columns  = ['ICD 1','ICD 2','ICD 3','Total'])
    ICDdf_totals = [x/sum(ICDdf['Total']) for x in list(ICDdf['Total'].values)]
        
    # calculate and record Pearson correlation coefficient
    pcc[sel] = pearsonr(f0,f)
    
    # calculate and record pvals for barcode enrichment
    binom_pvals = binomial(f0,x,q = 0.05)

    df_filtered['Binom p'] = [pval[0] for pval in binom_pvals[0]]
    df_filtered['Binom significant?'] = [pval[1] for pval in binom_pvals[0]]
    '''
    Wald_pvals = Waldlog(x0,x,q = 0.05)
    df_filtered['Wald p'] = [pval[0] for pval in Wald_pvals[0]]
    df_filtered['Wald significant?'] = [pval[1] for pval in Wald_pvals[0]]
    '''
    
    # calculate and record pvals for ICD enrichment
    binom_ICD_pvals = binomial(ICDdf_totals0,ICDdf['Total'])
    ICDdf['Binom p'] = [pval[0] for pval in binom_ICD_pvals[0]]
    ICDdf['Binom significant?'] = [pval[1] for pval in binom_ICD_pvals[0]]
    '''
    Wald_ICD_pvals = Waldlog(ICDdf0['Total'],ICDdf['Total'])
    ICDdf['Wald p'] = [pval[0] for pval in Wald_ICD_pvals[0]]
    ICDdf['Wald significant?'] = [pval[1] for pval in Wald_ICD_pvals[0]]
    '''
    
    # export pvals to Excel
    df_filtered.to_excel(writer, sheet_name = sel)
    ICDdf.to_excel(writer, sheet_name = sel + ' ICDs')
    writer.save()
    
    ''' FIXME
    # calculate and record pvals for ICD positional preference
    Wald_pos_pvals = positions(ICD_counts0,ICD_counts,test = 'Wald')
    Fisher_pos_pvals = positions(ICD_counts0,ICD_counts,test = 'Fisher')
    '''
    
    # assess ICD covariation
    ddGij,Pij,Pi,ax = SCA(df0,df,ICDs,rnd = sel,saveto = figure_folder)
    
    # assess ICD family enrichment
    GSEA(rnd = sel, f1 = ICDdf_totals0, f2 = ICDdf_totals,
         features = ICDs,sets = famdict,p = 1,saveto = figure_folder,q = 0.05)
    
    # plot volcano plots
    '''
    data = [[f0,f,binom_pvals],[f0,f,Wald_pvals],
            [ICDdf_totals0,ICDdf_totals,binom_ICD_pvals],
            [ICDdf_totals0,ICDdf_totals,Wald_ICD_pvals]]
    titles = ['Barcode bimomial pvals','Barcode Wald log pvals',
          'ICD binomial pvals','ICD Wald log pvals']
    '''
    data = [[f0,f,binom_pvals],
            [ICDdf_totals0,ICDdf_totals,binom_ICD_pvals]]
    titles = ['Barcode bimomial pvals','ICD binomial pvals']
    volcano(data,title = titles,rnd = sel,saveto = figure_folder)
    
    '''
    # plot circle heatmap of ICD enrichment for each generation of CAR
    idx = 0
    binom_ICDpos_pvals = np.ones(len(ICDs),positions)
    comparisons = comb(positions,2)
    for i in range(positions-1):
        for j in range(i+1,comparisons)):
            binom_ICDpos_pvals[:,i] = binomial(f1 = ICDdf0['ICD {}'.format(i+1)])
            idx += 1
    
    ICDdf,feat_dict,pvals,scaletitle,rnd,saveto,pvalmethod = 'binomial ',cols = ['ICD 1','ICD 2','ICD 3'],
            remove = False,vmin = 0, vmax = 4,bins = 5
    circle_heatmap(ICDdf,ICDdict,binom_ICD_pvals,scaletitle = 'log_2(fold change)',
                   rnd = sel,saveto = figure_folder, pvalmethod = 'binomial',
                   cols = ['ICD 1','ICD 2','ICD 3'],remove = False,vmin = 0, 
                   vmax = 4,bins = 5, pvals = binom_ICD_pvals)
    circle_heatmap(ICDdf,ICDdict,Wald_ICD_pvals,scaletitle = 'log_2(fold change)',
                   rnd = sel,saveto = figure_folder, pvalmethod = 'Wald log',
                   cols = ['ICD 1','ICD 2','ICD 3'],remove = False,vmin = 0, 
                   vmax = 4,bins = 5, pvals = binom_ICD_pvals)
    '''