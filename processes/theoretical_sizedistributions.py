#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:19:37 2020

@author: khloegordon
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

plt.rcParams["figure.figsize"] = (6,4)

path = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'
os.chdir(path)

df = pd.read_excel('ICD lists and sequences updated 06.24.21.xlsx')
print('{:d} unique ICDs'.format(len(df.index)))

#%% Plot theoretical ICD length distribution

plt.hist(df['Length (bp)'], color = 'forestgreen', density = 'True', bins = 10)
plt.ylabel('Frequency')
plt.xlabel('Length (bps)')
plt.title('ICD length distribution')
plt.show()

#%% Plot theoretical CAR length distribution %%#

from itertools import product

a = 307 # total length of constant regions in amplicons, including linkers and etc.
# 310 for pGGA cut with BamHI + PacI
# 387 for amplification from pHIV with IL13zk_seq2 and N_IRES_R
# 325 for PacBio samples
# 2210 for pHIV cut with EcoRI + XbaI
# 307 for CD19 CAR library

construct_lengths = product(df['Length (bp)'].values,
                               df['Length (bp)'].values,
                               df['Length (bp)'].values)
total_construct_lengths = [sum(combo)+a for combo in construct_lengths]
print('{:d} unique permutations'.format(len(total_construct_lengths)))

plt.hist(total_construct_lengths, color = 'forestgreen', density = True, bins = 10) # plot histogram

plt.ylabel('Frequency')
plt.xlabel('Length (bps)')
plt.title('Theoretical construct length distribution')
plt.show()