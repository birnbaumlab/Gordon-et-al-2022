#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import distance
import matplotlib.pyplot as plt
from itertools import combinations,combinations_with_replacement
import time
import os


# In[3]:


# Utility functions for parsing FASTQs
def get_ID_from_line(line):
    return line.strip().split()[0]

def get_seq_from_line(line):
    return line.strip()

def get_sequence_generator(fname):
    with open(fname,'r') as f:
        while True:
            for i in range(4):
                try:
                    line = next(f)
                except StopIteration:
                    return
                if i == 0: ID = get_ID_from_line(line)
                if i == 1: seq = get_seq_from_line(line)
                if i == 3: yield (ID, seq)

# Creates dictionary of sequences a max Hamming distance away from input sequence, 
# allowing clustering of similar barcodes
def get_mutation_dict(name,seq,hamming_distance = 1):
    """ makes a dictionary that maps sequence variants to the feature,hamming distance """
    mutation_dict = {seq: (name,0)} # set up empty dictionary

    def mutate_seq(seq,dist):
        """ adds mutated sequences of set hamming distance to list"""
        position_sets = list(combinations(range(len(seq)),dist))
        base_sets = list(combinations_with_replacement('AGCTN',dist))
        for positions in position_sets:
            for bases in base_sets:
                mutated_seq = seq # create copy of sequence
                #
                for pos,b in zip(positions,bases):
                    mutated_seq = mutated_seq[:pos] + b + mutated_seq[pos + 1:]
                # if there isn't a better sequence already stored
                if not mutated_seq in mutation_dict:
                    mutation_dict[mutated_seq] = (name,dist)

    for dist in range(1,hamming_distance + 1):
        mutate_seq(seq,dist)

    return mutation_dict


# In[4]:


# Returns BC; if barcode is not found, returns None
def find_barcode(seq, pre_hash, post_hash, bc_length, max_dist):
    
    # Constraining positions of flanking regions and BCs:
    pre = list(pre_hash.keys())[0]
    post = list(post_hash.keys())[0]
    BC_index = 37 # updated
    pre_index = 17
    post_index = 75
    
    pre_check = seq[pre_index:BC_index]
    post_check = seq[BC_index+bc_length:post_index]
    
    if pre_check in pre_hash and post_check in post_hash:
        return seq[BC_index:BC_index+bc_length]
    else:
        return None

'''
    max_index = len(seq)-len(pre)-bc_length-len(post)
    for i in range(max_index+1):
        current = seq[i:i+len(pre)]
        if current in pre_hash.keys():
            bc_start = i+len(pre)
            post_check = seq[bc_start+bc_length:bc_start+bc_length+len(post)]
            if post_check in post_hash.keys():
                bc = seq[bc_start:bc_start+bc_length]
                return [bc, current, post_check, seq]


    return None
'''


# In[14]:


# Returns a dictionary of {barcodes : counts}
# Inputs = input_file name, 5' flanking region, 3' flanking region, bc_length, hamming distance for flanking regions
def count_barcodes(input_file, bc_5p, bc_3p, bc_length, max_dist):

    bc_5p_hash = get_mutation_dict(bc_5p, bc_5p, max_dist)
    bc_3p_hash = get_mutation_dict(bc_3p, bc_3p, max_dist)

    bc_freq = {}
    seq_count = 0
    flank_count = 0
    t0 = time.process_time()

    for ID,seq in get_sequence_generator(input_file):
        seq_count += 1

        if seq_count % 1000000 == 0:
            print(str(seq_count) + ' sequences processed. (Elapsed time: ' + str(time.process_time()-t0) + ' secs.)')
            
        bc = find_barcode(seq, bc_5p_hash, bc_3p_hash, bc_length, max_dist)
        if bc != None:
            flank_count += 1
            if bc in bc_freq:
                bc_freq[bc] += 1
            else:
                bc_freq[bc] = 1
            '''
                collapsed = False
                for k in bc_freq:
                    if distance.hamming(k,bc) < dist_threshold:
                        bc = k
                        collapsed = True
                        break
                if collapsed:
                    bc_freq[bc] += 1
                else:
                    bc_freq[bc] = 1
            '''
    
    bc_purity = (flank_count / seq_count) * 100.0
    
    print(str(seq_count) + ' total sequences processed.')
    print('Flanking regions found in ' + str(round(bc_purity,2)) + '% of sequences.')
    
    df = pd.Series(bc_freq)
    return df.sort_values(ascending = False)
    #df = pd.DataFrame.from_dict(bc_freq, orient='index', columns = ['Counts'])
    #return df.sort_values(by = 'Counts', ascending = False)


# In[15]:


# Collapse similar barcodes post-hoc
# Input = barcode dataframe/series containing sorted raw counts of uncollapsed barcodes, and Hamming dist threshold
# Returns sorted dataframe containing collapsed barcode counts and dictionary of {collapsed BC : [original BCs]}
def collapse_barcodes(barcode_df, max_dist):
    barcodes = list(barcode_df.index) # retains sorted order
    bc_dict = barcode_df.to_dict()
    i = 0
    collapsed_total = []
    collapsed_dict = {}
    while i + len(collapsed_total) < len(barcodes):
        bc = barcodes[i]
        if bc not in collapsed_total:
            similar_hits = []
            similar_bcs = list(get_mutation_dict(bc, bc, max_dist).keys())[1:]
            for similar_bc in similar_bcs:
                if similar_bc in bc_dict:
                    similar_hits.append(similar_bc)
                    similar_bc_count = bc_dict.pop(similar_bc)
                    bc_dict[bc] += similar_bc_count
                    collapsed_total.append(similar_bc)
            if len(similar_hits) > 0:
                collapsed_dict[bc] = similar_hits
        i += 1
    print(str(len(collapsed_total)) + ' barcodes collapsed.')
    print(str(len(barcodes)) + ' starting barcodes folded into ' + str(len(bc_dict)) + '.')
    df = pd.DataFrame.from_dict(bc_dict, orient='index', columns = ['Counts'])
    collapsed_bc_df = pd.DataFrame.from_dict(collapsed_dict, orient='index')
    return df.sort_values(by = 'Counts', ascending = False), collapsed_bc_df


# In[16]:


# Save output file
def save_output(sample, df, suffix = None):
    if suffix == None:
        suffix = ''
    else:
        suffix = '_' + suffix
    output_file = './datasets/' + sample + suffix + '.xlsx'
    df.to_excel(output_file)
    return output_file


# In[17]:


### MAIN ###

#20 bp flanking regions, max Hamming = 1 for flanking regions and barcode collapsing
bc_5p = 'GGAGAACCACCTTGTTGG'
bc_3p = 'GTTTAAGAGCTAA' # updated with Taeyoon's BC flanking sequences
dist_threshold = 1
bc_length = 18

cwd = os.getcwd()
for filename in os.listdir(cwd):
    if filename.endswith('.fastq'):
        print('Processing file: ' + filename)
        uncollapsed_bcs = count_barcodes(filename, bc_5p, bc_3p, bc_length, dist_threshold)
        collapsed_bc_counts, collapsed_bcs = collapse_barcodes(uncollapsed_bcs, dist_threshold)
        output_file = save_output(filename, collapsed_bc_counts)
        print('Collapsed barcode counts saved in ' + output_file + '.')


# In[ ]:


# Checking collapsing results
check = True
for index, row in collapsed_bcs.iterrows():
    for i in range(len(row)):
        if row[i] != None:
            check = distance.hamming(index,row[i]) <= 1
print(check)


# In[ ]:


# Uncollapsed analysis
### MAIN ###

#20 bp flanking regions, max Hamming = 1 for flanking regions and barcode collapsing
bc_5p = 'TTGGAGAACCACCTTGTTGG'
bc_3p = 'TTTGTACCCCGTATTCCGTT'
dist_threshold = 1
bc_length = 18

cwd = os.getcwd()
for filename in os.listdir(cwd):
    if filename.endswith('.fastq'):
        print('Processing file: ' + filename)
        uncollapsed_bcs = count_barcodes(filename, bc_5p, bc_3p, bc_length, dist_threshold)
        collapsed_bcs = collapse_barcodes(uncollapsed_bcs, dist_threshold)
        output_file = save_output(filename, uncollapsed_bcs)
        print('Collapsed barcode counts saved in ' + output_file + '.')


# In[15]:


# Returns table of barcode counts, frequencies, and BC ranks for input sample
# By default processes collapsed reads, without a BC rank cutoff
def read_barcode_table(sample, collapsed = True, rank_cutoff = None):
    if collapsed:
        sample += '_collapsed'
    input_file = './output/' + sample + '.xlsx'
    df = pd.read_excel(input_file, index_col = 0, header = 0, usecols = 'A:D', nrows = rank_cutoff)
    return df


# In[8]:


# Comparison of selected to control barcodes
# Input: Sample names (as strings) of selected and control samples, Hamming distance cutoff to define overlapping
# barcodes, and number of top barcodes in selected sample to use as a query (default = None, processes all BCs)
# Returns dataframe of: shared BCs, read counts, frequencies, and BC ranks in both selected and control datasets
def barcode_overlap(selected_sample, control_sample, max_dist, read_cutoff, rank_cutoff = None, collapsed = True):
    selected_df = read_barcode_table(selected_sample, collapsed = collapsed, rank_cutoff = rank_cutoff)
    control_df = read_barcode_table(control_sample, collapsed = collapsed)
    selected_bcs = selected_df.to_dict(orient='index')
    control_bcs = control_df.to_dict(orient='index')
    shared_bcs = []
    for bc in selected_bcs:
        similar_bcs = list(get_mutation_dict(bc, bc, max_dist).keys())
        for bc_check in similar_bcs:
            if bc_check in control_bcs:
                selected_count = selected_df['Counts'][bc]
                selected_freq = selected_df['Freq'][bc]
                selected_rank = selected_df['BC Num'][bc]
                control_count = control_df['Counts'][bc_check]
                control_freq = control_df['Freq'][bc_check]
                control_rank = control_df['BC Num'][bc_check]
                if selected_count >= read_cutoff and control_count >= read_cutoff:
                    shared_bcs.append([bc, selected_count, selected_freq, selected_rank, 
                                       bc_check, control_count, control_freq, control_rank])
    output = pd.DataFrame(data = shared_bcs, columns = ['Selected BC', 'Selected Counts', 'Selected Freq', 'Selected Rank',
                                                       'Control BC', 'Control Counts', 'Control Freq', 'Control Rank'])
    return output.sort_values(by='Selected Rank')


# In[9]:


# Comparison of selected to control barcodes (Hamming distance of 1)
# Returns dataframe of: shared BCs, read counts, frequencies, and BC ranks in both selected and control datasets
rank_cutoff = 1000
hamming_cutoff = 1
unselected_1 = 'D21-9816'
selected_1 = 'D21-9817'
unstim_1 = 'D21-9818'
unselected_2 = 'D21-9819'
selected_2 = 'D21-9820'
unstim_2 = 'D21-9821'

selected_v_unstim_2 = barcode_overlap(selected_2, unstim_2, hamming_cutoff, read_cutoff = 5)
selected_v_unstim_2
save_output('Rep2_selected_unstim_overlap', selected_v_unstim_2)


# In[18]:


selected_2_v_1 = barcode_overlap(selected_2, selected_1, hamming_cutoff, read_cutoff = 5)
selected_2_v_1
save_output('Rep1_selected_Rep2_selected_overlap', selected_2_v_1)


# In[19]:


unstim2_v_stim2 = barcode_overlap(unstim_2, selected_2, hamming_cutoff, read_cutoff = 5)
unstim2_v_stim2
save_output('Rep2_unstim_v_stim', unstim2_v_stim2)


# In[20]:


unstim2_v_stim1 = barcode_overlap(unstim_2, selected_1, hamming_cutoff, read_cutoff = 5)
unstim2_v_stim1
save_output('Rep2_unstim_v_Rep1_stim', unstim2_v_stim1)


# In[ ]:


triple_overlap = selected_2_v_1.merge(selected_v_unstim_2, on = 'Selected BC')
triple_overlap


# In[ ]:


triple_common_bcs = list(triple_overlap['Selected BC'])
shared_selection_bcs = selected_2_v_1
for common_bc in triple_common_bcs:
    index = shared_selection_bcs[shared_selection_bcs['Selected BC'] == common_bc].index.values
    shared_selection_bcs = shared_selection_bcs.drop(index,'index')
shared_selection_bcs
save_output('selection_1_and_2_overlap_without_unstim', shared_selection_bcs)


# In[ ]:




