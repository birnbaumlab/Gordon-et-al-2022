#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 09:17:59 2022

@author: khloegordon
"""
# standard libraries
import pickle

# nonstandard libraries
import openpyxl

# homegrown libraries
from methods.loader import get_generator, create_data_generator,filter_sequence
from methods.utilities import hamming_dist,levenshtein_dist
from methods.output2 import publish_reference_to_excel

def merge_references(domain2seq,bc2seq,**kwargs):
    #domain2seq is a dictionary with ICDs as keys and DNA sequence as values
    #bc2seq is a dictionary with bcs as keys and each PacBio sequence containing that barcode

    settings = {
            'dist_function':      'hamming',
            'dist_threshold':           0.8,
            }

    # use user settings
    if 'dist_function' in kwargs: 
        settings['dist_function'] = kwargs['dist_function']
    if 'dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['dist_threshold']
    
    # select distance function
    if settings['dist_function'] == 'levenshtein': dist_func = levenshtein_dist
    elif settings['dist_function'] == 'hamming':   dist_func = hamming_dist 
    else: raise KeyError('Distance function not recognized! Try levenshtein or hamming.')
    
    threshold = settings['dist_threshold'] # get threshold

    reference = {} # initialize dictionary with bc as key, ICD sequences as values
    min_seq = min([len(seq) for seq in domain2seq.values()]) # find length of shortest ICD

    for bc_count,(bc,seqs) in enumerate(bc2seq.items()): # iterate through bcs identified in PacBio
            
        print('{}/{} sequences processed...'.format(bc_count,len(bc2seq)))

        '''# premature termination bridge
        if bc_count > 10: 
            break
        #'''#

        reference[bc] = [] # initialize list of domain combinations

        for seq_count,seq in enumerate(seqs): # iterate over each PacBio sequence that was identified with that bc
            
            matched_domains_with_location = [] # create list with items [matched ICD, length of ICD, location of matched ICD, back of sliding window]
            domain_list = [] # create list of final matched ICD combinations for each sequence
            
            if seq_count >= 25: break # don't need to iterate over more than 25 sequence occurrences of a barcode
            front,domain_list = 0,[] # start scanning at front
            for back in range(min_seq,len(seq)): # back of sliding window  will start from the minimum ICD length to the end of the sequence
                if back - front < min_seq: continue # if segment smaller than smallest domain

                # get the domains that matched within given threshold, scanning within window [front:back]
                matched_domains = _endswith(seq[front:back],domain2seq,dist_func,threshold)

                # iterate through each matched domain and create an entry with location
                for (domain,length) in matched_domains:
                    matched_domains_with_location += [(domain,length,back - length,back)]

            # sort entries by start position, then by length (ascending, descending, respectively)
            matched_domains_with_location.sort(key = lambda x: (x[2],-x[1]))

            # take the largest domains that do not overlap and save in domain_list
            current_pos = -1
            for (domain,length,start,end) in matched_domains_with_location:
                if start > current_pos:
                    domain_list.append(domain)
                    current_pos = end
                    
            #print('matched domains:',matched_domains_with_location)
            #print('domain list:',domain_list)

            # if at least one sequence was aligned, add to reference list with barcode as key
            if len(domain_list) > 0:
                reference[bc].append(tuple(domain_list))

        #print(bc,': ',reference[bc])

    return reference

#---------------------------------#

def _endswith(seq,domain_dict,dist_func,threshold):
    """ Checks for domain matching end of input sequence at certain hamming distance """
    matched_domains = []
    
    for domain,domain_seq in domain_dict.items(): # loop through all ICDs and their sequences
        if len(seq) < len(domain_seq): # if current window is shorter than ICD length, skip
            continue
        if 1. - float(dist_func(seq[-len(domain_seq):],domain_seq))/len(domain_seq) >= threshold: 
            # find distance between back end of sequence in window that is the length of current ICD from that of the ICD sequence 
            # normalize number of edits to length of ICD sequence to get error rate
            # subtract from 1 to get % match and compare to threshold
            # if above threshold, add that ICD and it's length to matched_domains
            matched_domains += [(domain,len(domain_seq))]
    
    return matched_domains