#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:48:12 2022

@author: khloegordon
"""

# standard libraries
import os

# nonstandard libraries

# custom libraries
from tests import barcode_count
from tests import create_reference
from processes import levenshtein_clustering
from processes import domain_mapping
from processes import process_datasets
from processes import track_enrichment

def main():
    
    print('Initializing DomainSeq...')

    """ Main overhead function to run pipeline """

    """ --------------------- """ 
    """ MODIFY ANYTHING BELOW """
    """ --------------------- """ 

    base = '/Users/khloegordon/Dropbox (MIT)/Code/DomainSeq/'
    os.chdir(base)

    # must be placed in fastq files folder
    fastq_filenames = [base + 'fastq files/' + fname for fname in (
                        ''
                        )]
    # note: unmapped Illumina sequencing data can be de-multiplexed by 
    # round-specific indices using demultiplex.py in /processes

    # dataset and reference generation settings
    # mapping settings
    settings = {
            
            # general settings
            'overwrite':                                          False,
            'silent':                                             False,
            
            ### --- dataset settings --- ###
            #'fastq_filenames':                         fastq_filenames,
            'read':                                                   1,
            'barcode_5p':                                            '',
            'barcode_3p':                                            '',
            'pre_index':                                             38, 
            # distance from beginning of read to beginning of 5p flanking region; 
            # use reverse complement if read 2
            'barcode_length':                                        18,
            'barcode_dist_function':                          'hamming',
            'barcode_dist_threshold':                                 1, # allowable edit distance
            'levenshtein_threshold':                                  0, # set to 0 to skip
            'count_threshold':                                       10,
            'dataset_sheetnames':                                  [''],
            # dataset filename, if already processed -- overrides dataset settings
            # must be placed in datasets folder
            'dataset_filename':                   
                                                                     '',
            
            ### --- reference settings --- ###
            'sequence_reference_filename':
                                                                     '', 
            'domain_reference_filename': 
                                                                     '', 
            'ref_dist_function':                              'hamming',
            'ref_dist_threshold':                                   0.8, # percent match, NOT edit distance
            'pre_amplicon':                      'GCATCTTCTTCGGCGGAGTG',
            'post_amplicon':                     'CAGATTGGACTGGTCACCTT',
            # reference filename, if already processed -- overrides reference settings
            # must be placed in references folder
            'reference_filename':                                    '',
            
            ### --- mapping settings --- ###
            # domain to barcode mapping edit distance / levenshtein threshold
            'mapping_distance':                                       2,
            'ICD_positions':                                          3,
    
            ### --- analysis and visuals settings --- ###
            # mapped dataset filename, if already processed -- overrides domain_mapping
            # must be placed in processed datasets folder
            'mapped_dataset_filename':                              
                                                                     '',
            'normalize_to':                                          '',
            # adds CAR generation column to excel sheet)
            'split_by_CAR_generation':                             True, 
            'frequency_heatmap':                                   True,
            'vmin':                                                   0,
            'vmax':                                                   4,
            'log2FC_heatmap':                                      True,
            'pie_plot':                                            True,
            
            ### -- tracking settings --- ###
            'tracking':                                            True,
            # processed dataset filename, if already processed -- 
            # will override output of processed dataset in a given run
            'processed_dataset_filename':
                                                                     '',
            'consensus domains':                          ['4-1BB,CD3z',
                                                            'CD28,CD3z',
                                                     'CD28,4-1BB,CD3z'],
            'domains':                                         ['4-1BB',
                                                                 'CD28',
                                                                'CD3z'],
            'barcodes':                                              [],
            # track enrichment of top X barcodes [0] in a specific round [1] of 
            # selection individually based on 'Frequency' or 'Fold change' [2]
            'topBCs':                      (10,'2nd CD69+','Frequency'),
            # track aggregate proportion of enrichment of top Y barcodes
            'topN':                                   [10,25,50,100,500]
            }

    """ --------------------- """ 
    """ MODIFY ANYTHING ABOVE """
    """ --------------------- """ 
    
    # collect barcode frequencies from Illumina data
    if settings['dataset_filename']:
        dataset_path = './datasets/' + settings['dataset_filename']
        print('Using {} as dataset'.format(dataset_path))
    else:
        print('Counting barcode frequencies...')
        _, dataset_path = barcode_count.run(fastq_filenames,**settings)
        if settings['levenshtein_threshold'] > 0:
            _, dataset_path = levenshtein.run(dataset_path,settings['dataset_filename'],**settings)
        print('Finished! Dataset exported to {}'.format(dataset_path))

    # identify domains in PacBio data
    if settings['reference_filename']:
        reference_path = './references/' + settings['reference_filename']
        print('Using {} as reference'.format(reference_path))
    else:
        print('Generating reference...')
        _, reference_path = create_reference.run(**settings)
        print('Finished! Reference exported to {}'.format(reference_path))

    # map domains to barcodes
    if settings['mapped_dataset_filename']:
        mapped_path = './processed datasets/' + settings['mapped_dataset_filename']
        print('Using {} as mapped dataset'.format(mapped_path))
    else:
        print('Mapping barcode frequencies to PacBio reference...')
        _,mapped_path = domain_mapping.run(dataset_path,reference_path,**settings)
        print('Finished! Mapped dataset exported to {}'.format(mapped_path))
    
    # analyze ICD frequencies and generate plots
    if settings['processed_dataset_filename']:
        processed_path = './processed datasets/' + settings['processed_dataset_filename']
        print('Using {} as processed dataset'.format(processed_path))
    else:
        print('Processing dataset...')
        processed_path,figure_folder = process_datasets.run(mapped_path,**settings)
        print('Finished! Exported to {}. Figures exported to {}.'.format(processed_path[0],figure_folder))
        
    # track enrichment
    if settings['tracking']:
        print('Tracking enrichment...')
        tracked_path = track_enrichment.run(processed_path,**settings)

def _check_settings(settings):
    type_errors = {
            'overwrite':bool,
            'silent':bool,
            'clustering_threshold':int,
            'count_threshold':int,
            'barcode_5p':str,
            'barcode_3p':str,
            'sequence_reference_filename':str,
            'domain_reference_filename':str,
            'dist_function':str,
            'dist_threshold':float,
            }
          
    # general case
    for var,var_type in type_errors.items():
        if not isinstance(settings[var],type_errors[var]):
            raise TypeError('settings[{}] is not {}'.format(var,var_type))

    # special case
    if not isinstance(settings['dataset_filenames'],dict):
        raise TypeError('settings[dataset_filenames] is not {}'.format(dict))

    for key,val in settings['dataset_filenames'].items():
        if not isinstance(settings['dataset_filenames'][key],(list,tuple)):
            raise TypeError('dataset_filenames[{}] is not {}'.format(key,list))

'''

ARCHIVED CODE

    """
    There exists four major methods here:

        load_dataset.run(**settings)

            Inputs: keyword arguments
                clustering_threshold (default: 3)
                silent (default: False)
                dataset_filenames (default: dataset_filenames)
                count_threshold (default: 10)
                barcode_5p (default: None)
                barcode_3p (default: None)
                overwrite (default: False)

            Outputs: 

        create_reference.run(**settings)
                
            Inputs: keyword arguments
                sequence_reference_filename (default: 'CCS_Run2.fastq')
                domain_reference_filename (default: 'ICD lists and sequences.xlsx')
                dist_function (default: 'hamming')
                dist_threshold (default: 0.8)

            Outputs: 

    """
    
    _check_settings(settings)
    settings['results_directory'] =_get_results_directory()

    print('Collecting datasets...')
    datasets = load_dataset.run(**settings)
    print('Finished!')

    #comparison = compare_datasets.run(datasets,reference=reference,**settings)
    print('Finished!')

    print('Generating visuals...')
    #generate_visuals.run(datasets=datasets,reference=reference,comparison=comparison,**settings)
    print('Finished!'

    def _get_results_directory():
    
        """ Generates a directory with a unique index """
    
        directory = 'results'
        overwrite = False
    
        if not os.path.isdir(directory):
            print('Making {} directory...'.format(directory))
            os.mkdir(directory)
    
        # change in directory
        home = os.getcwd()
    
        for i in range(1,10000):
    
            results_folder = os.path.join(home,directory,str(i).zfill(4))
    
            if os.path.isdir(results_folder): continue
            else: 
                os.mkdir(results_folder) # KG added
                break
    
        return results_folder
'''
# if you call "python main.py", this will run
if __name__ == "__main__":
    
    main()
