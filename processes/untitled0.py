#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 23:59:32 2021

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

path = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/processed datasets/'

def count2freq(ILdf):
    '''converts counts to frequencies normalized by all counts in a round of selection'''
    
    ILdf['Frequency'] = ''
    
    #total = ILdf['Count'].sum(axis = 0)
    ILdf['Frequency'] = ILdf['Count'].divide(ILdf['Count'].sum(axis = 0))
            
    return ILdf


# import unselected Illumina data
ILdf0 = pd.read_excel(path + 'dataset_mapped_0000_CPcounts.xlsx',sheet_name = 'Unselected')
ILdf0 = count2freq(ILdf0) # convert counts to frequencies

# create new Excel spreadsheet
writer = pd.ExcelWriter(path + 'Unselected CP counts.xlsx', engine = 'xlsxwriter')
ILdf0.to_excel(writer,
               sheet_name = 'Unselected')
