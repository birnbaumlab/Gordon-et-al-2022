#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:48:24 2019

@author: khloegordon
"""

import os
import pandas as pd
import numpy as np
import xlsxwriter
import editdistance
import distance


"""
# Workflow for code

my_list_of_rows = load_data_from_excel(filename)
my_list_of_illumina_bcs = get_illumina_bcs()

# list comprehension
my_bc_data = dict([(bc,[]) for bc in my_list_of_illumina_bcs])

# iterate through each row
for index,row in enumerate(my_list_of_rows):
    # iterate through each bc
    for bc in my_list_of_illumina_bcs:
        # if the distance is below threshold
        if levenshtein(row[my_bc],bc) < 5: # some arbitrary number
            my_bc_data[bc].append(('seq{}'.format(index),row[my_bc]))

# something like the code from above, that writes a dictionary to separate FASTA files
save_bc_data(my_bc_data)

#"""

def barcodes2fasta(path,file1,sheet1,file2,sheet2,criteria):
    
    ''' searches for every instance of barcodes in file 1 (Illumina data) in 
    file 2 (PacBio data) given a criteria dictionary containing keywords and 
    values that describe sequences of interest
    
    keywords = selection ( scheme), round (of selection), 
    top enriched (cut off)''' 

    os.chdir(path)
    
    df1 = pd.read_excel(file1,
                        sheet_name = sheet1)
        
    df2 = pd.read_excel(file2,
                        sheet_name = sheet2)
    
    top = int(criteria['top enriched'])
    
    data = df1[df1['Dataset'] == criteria['selection']]
    data = data[data['Round'].values == criteria['round']]
    data.sort_values(by = ['Count'], ascending = False)
    
    bc_list = data['BC'].astype(str).values.tolist()
    bc_list = bc_list[0:top]  
    bc_dict = {key: [] for key in bc_list}
    #bc_dict = dict.fromkeys(bc_list,[])
    
    for bc_key in bc_dict:
        for i,row in df2.iterrows():
            if editdistance.eval(bc_key,row['Barcode']) <= criteria['n_edits']:
                seq_str = row['Sequence']
                seq_list = seq_str.split(';')
                for seq in seq_list:
                    bc_dict[bc_key].append(seq)
                    
    matched = 0
    for bc_key in bc_dict:
        if bc_dict[bc_key]:
            matched += 1
    print('\n{}/{} top {} enriched barcodes found in PacBio data\n'.format(matched,len(bc_list),criteria['top enriched']))

                   
    '''
    print('\nbarcodes2fasta dictionary:')                
    for k,v in bc_dict.items():
        print('{}: {}'.format(k,v))
    '''
    
    path = path + '/Nested BCs'
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        os.chdir(path)
            
    count = 0
    total = len(bc_list)
    for bc,seqs in bc_dict.items():
        if seqs:
            count += 1
            with open('{}.fasta'.format(bc),'w') as f:
                for ii,seq in enumerate(seqs):
                    f.write('> sequence{}\n{}\n\n'.format(ii,seq))
                f.close
            print('{} FASTA files written!'.format(count))
                        
def findnestedBCs(path1,file1,sheet1,path2,file2,sheet2,criteria,exportpath):
    
    os.chdir(path1)
    
    df1 = pd.read_excel(file1,
                        sheet_name = sheet1)
    
    os.chdir(path2)
    
    df2 = pd.read_excel(file2,
                        sheet_name = sheet2)
        
    data = df1[df1['Dataset'] == criteria['selection']]
    data = data[data['Round'].values == criteria['round']]
    data.sort_values(by = ['Count'], ascending = False)

    bc_list = data['BC'].astype(str).values.tolist()
    bc_dict = {key: [] for key in bc_list}
    #bc_dict = dict.fromkeys(bc_list,[])
            
    for bc_key in bc_dict:
        for i,row in df2.iterrows():
            if editdistance.eval(bc_key,row['Barcode']) <= criteria['n_edits']:
                bc_dict[bc_key].append(row['Barcode'])
    
    print('\nfindnestedBCs dictionary:')                
    for k,v in bc_dict.items():
        print('{}: {}\n'.format(k,v))
        
    matched = 0
    for bc_key in bc_dict:
        if bc_dict[bc_key]:
            matched += 1
    print('\n{}/{} Illumina barcodes found in PacBio data\n'.format(matched,len(bc_list)))

    
    """ Publish barcodes to excel .xlsx """
    
    exportpath = exportpath + '/Nested BCs'
    if not os.path.exists(exportpath):
        os.makedirs(exportpath)
    else:
        os.chdir(exportpath)

    workbook = xlsxwriter.Workbook('Nested barcodes.xlsx')
    worksheet = workbook.add_worksheet('DATA')
    cell_format = workbook.add_format()
    cell_format.set_bold()

    worksheet.write(0,0,'Settings',cell_format)

    row = 0
    col = 0
    
    for k,v in criteria.items():
        row += 1
        worksheet.write(row,0,k)
        worksheet.write(row,1,v)
    
    row += 2    
    worksheet.write(row,0,'Enriched Barcode',cell_format)
    worksheet.write(row,1,'Nested Barcodes',cell_format)


    for key in bc_dict.keys():
        row += 1
        col = 0
        worksheet.write(row, col, key)
        for item in bc_dict[key]:
            col += 1
            worksheet.write(row, col, item)
    
    worksheet.save()
    workbook.close()
    
    print('Nested BCs exported to Excel!')
    
def levenshtein(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in range(size_x):
        matrix [x, 0] = x
    for y in range(size_y):
        matrix [0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    return (matrix[size_x - 1, size_y - 1])

'''
# test levenshtein algorithm
def unit_test():
    
    a = 'AACTCGATCG'
    b = 'TAACACGATCG'
    dist1 = levenshtein(a,b)
    dist2 = editdistance.eval(a,b)
    dist3 = distance.levenshtein(a,b)

    if dist1 < 5:
        print('\nLevenshtein function works')
    else:
        print('\nLevenshtein function does not work')
        
    if dist2 < 5:
        print('\nEditdistance works')
    else:
        print('\nEditdistance does not work')
        
    if dist3 < 5:w
        print('\nDistance.levenshtein works')
    else:
        print('\nDistance.levenshtein does not work')

unit_test()'''
    
# execute search and create FASTA files
settings = {'selection': 'CD69+PD', 
            'round': 2, 
            'top enriched': 1000,
            'n_edits': 2}

'''
barcodes2fasta(path = '/Users/KhloeGordon/Dropbox (MIT)/Research/Code/DomainSeq/',
               file1 = 'dataset_0001.xlsx',
               sheet1 = 'DATA',
               file2 = 'reference_0001.xlsx',
               sheet2 = 'DATA',
               criteria = settings)

print('\nbarcodes2fasta finished!')
'''

findnestedBCs(path1 = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/datasets',
               file1 = 'dataset_0001.xlsx',
               sheet1 = 'DATA',
               path2 = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/references',
               file2 = 'reference_0001-PB1 run.xlsx',
               sheet2 = 'DATA',
               criteria = settings,
               exportpath = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq')

print('\nfindnestedBCs finished!')