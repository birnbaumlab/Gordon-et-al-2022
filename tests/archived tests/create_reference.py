
# standard libraries
import os,glob
import pickle

# nonstandard libraries
import openpyxl

# homegrown libraries
from methods.loader import get_generator, create_data_generator,filter_sequence,get_xlsx_settings
from methods.utilities import hamming_dist,levenshtein_dist
from methods.output2 import publish_reference_to_excel # KG changed from output
#from methods.reference import *
from methods.reference3 import * # KG changed from reference



def run(**kwargs):

    directory = './references' 
    overwrite = False

    if not os.path.isdir(directory):
        print('Making {} directory...'.format(directory))
        os.mkdir(directory)

    settings = {
            'sequence_reference_filename':             'CCS_Run2.fastq', 
            'domain_reference_filename': 'ICD lists and sequences.xlsx', 
            'dist_function':                                  'hamming',
            'dist_threshold':                                       0.8,
            }

    # adjust filtering settings
    if 'sequence_reference_filename' in kwargs: 
        settings['sequence_reference_filename'] = kwargs['sequence_reference_filename']
    if 'domain_reference_filename' in kwargs: 
        settings['domain_reference_filename'] = kwargs['domain_reference_filename']
    if 'ref_dist_function' in kwargs: 
        settings['dist_function'] = kwargs['ref_dist_function'] 
    if 'ref_dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['ref_dist_threshold'] 
    if 'overwrite' in kwargs: 
        overwrite = kwargs['overwrite']
    if 'barcode_5p' in kwargs:
        pre = kwargs['barcode_5p']
    if 'barcode_3p' in kwargs:
        post = kwargs['barcode_3p']
    if 'pre_amplicon' in kwargs:
        start = kwargs['pre_amplicon']
    if 'post_amplicon' in kwargs:
        stop = kwargs['post_amplicon']


    # check for existing reference objects
    if overwrite == False:
        for fname in os.listdir(directory):

            if not fname.endswith('.xlsx'): continue # escape sequence for non-excel

            wb = openpyxl.load_workbook(os.path.join(directory,fname))

            loaded_settings = _get_settings(wb)
            
            if loaded_settings == settings: # if match between current settings
                decision = input('Discovered existing reference file {}, use? (Y/n)'.format(fname))
                if decision.upper() == 'Y' or decision == '':
                    reference = _get_data(wb)
                    print('\nReference loaded with {} barcodes.'.format(len(reference)))
                    return reference
                
    
    domain2seq = domain_reference_dictionary(settings['domain_reference_filename'])
    domain2seq = {k:v for k,v in domain2seq.items() if k is not None}
    bc2seq = sequence_reference_dictionary('./fastq files/' + settings['sequence_reference_filename'],pre,post)
    #bc2seq = sequence_reference_dictionary2('./fastq files/' + settings['sequence_reference_filename'],pre,post,start,stop)
    # FIX : sequence_reference_dictionary2 uses quality scores from get_generator2 in loader.py 
    # to sort sequences and trims amplicon region of interest for ICD scanning in sequence_reference_dictionary2
    bc2domain = merge_references(domain2seq,bc2seq,**settings)
    
    reference_dict = {
            'DATA':bc2domain,
            'SETTINGS':settings,
            }
        
    # publish to excel
    fname = publish_reference_to_excel(reference_dict,directory=directory)

    # return result
    return reference_dict, fname


def _get_settings(wb):

    loaded_settings = {}
    # pull single sheet from xlsx 
    ws = wb['SETTINGS']
    # extract settings from file
    for row in ws: loaded_settings[row[0].value] = row[1].value
    return loaded_settings


def _get_data(wb):

    loaded_data = {}
    # pull single sheet from xlsx 
    ws = wb['DATA']

    # extract settings from file
    for i,row in enumerate(ws): 
        if i == 0:
            index = [cell.value for cell in row].index('All Domain Sets')
            continue
        loaded_data[row[0].value] = [cell.value for cell in row[index:] if cell.value != None]

    return loaded_data


