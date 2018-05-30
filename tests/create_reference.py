
# standard libraries
import os,glob
import pickle

# nonstandard libraries
import openpyxl

# homegrown libraries
from methods.loader import create_data_generator,filter_sequence,get_xlsx_settings
from methods.utilities import hamming_dist,levenshtein_dist
from methods.output import publish_reference_to_excel
from methods.reference import *



def run(**kwargs):

    directory = './references' 
    overwrite = False

    if not os.path.isdir(directory):
        print('Making {} directory...'.format(directory))
        os.makedir(directory)

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
    if 'dist_function' in kwargs: 
        settings['dist_function'] = kwargs['dist_function']
    if 'dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['dist_threshold']
    if 'overwrite' in kwargs: 
        overwrite = kwargs['overwrite']


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
    bc2seq =     sequence_reference_dictionary(settings['sequence_reference_filename'])

    bc2domain = merge_references(domain2seq,bc2seq,**settings)
    
    reference_dict = {
            'DATA':bc2domain,
            'SETTINGS':settings,
            }
    
    # publish to excel
    publish_reference_to_excel(reference_dict,directory=directory)

    # return result
    return reference_dict


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


