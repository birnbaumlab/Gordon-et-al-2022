
# standard libraries
import os

# nonstandard libraries

# homegrown libraries
from methods.enrichment import *
from methods.loader import *
from methods.output import *
from methods.visuals import *

def run(**kwargs):

    directory = './datasets/'
    overwrite = False

    if not os.path.isdir(directory):
        print('Making {} directory...'.format(directory))
        os.mkdir(directory)

    # datasets submited as dictioanry of lists
    dataset_filenames = {
            'Selection Condition 1': ['non-selected.rtf','selected.rtf']
            }

    settings = {
            'clustering_threshold':                                   3,
            'silent':                                             False,
            'dataset_filenames':                      dataset_filenames,
            'count_threshold':                                       10,
            'barcode_5p':                                          None,
            'barcode_3p':                                          None,
            }

    # update dictionary
    if 'clustering_threshold' in kwargs: 
        settings['clustering_threshold'] = kwargs['clustering_threshold']
    if 'count_threshold' in kwargs: 
        settings['count_threshold'] = kwargs['count_threshold']
    if 'dataset_filenames' in kwargs: 
        settings['dataset_filenames'] = kwargs['dataset_filenames']
    if 'silent' in kwargs: 
        settings['silent'] = kwargs['silent']
    if 'barcode_5p' in kwargs: 
        settings['barcode_5p'] = kwargs['barcode_5p']
    if 'barcode_3p' in kwargs: 
        settings['barcode_3p'] = kwargs['barcode_3p']
    if 'overwrite' in kwargs: 
        overwrite = kwargs['overwrite']

    # check for existing reference objects
    if overwrite == False:
        for fname in os.listdir(directory):

            if not fname.endswith('.xlsx'): continue # escape sequence for non-excel

            wb = openpyxl.load_workbook(os.path.join(directory,fname))

            loaded_settings = get_xlsx_settings(wb)

            # weird corner case >.<
            if 'TRUE' in loaded_settings['silent']:
                loaded_settings['silent'] = True 
            if 'FALSE' in loaded_settings['silent']:
                loaded_settings['silent'] = False

            if loaded_settings == settings: # if match between current settings

                decision = input('Discovered existing reference file {}, use? (Y/n)'.format(fname))

                if decision.upper() == 'Y' or decision == '':
                    datasets_dict = _get_data(wb)
                    print('\nDataset loaded containing {} selections:'.format(len(datasets_dict)))
                    for selection_name,datasets in datasets_dict.items():
                        print('-- {} across {} rounds:'.format(selection_name,len(datasets)))
                        for i,dataset in enumerate(datasets):
                            print('---- Round {} - {} barcodes'.format(i,len(dataset)))
                            
                    return datasets_dict 

    print ('Starting analysis...')

    dataset_generators = {}
    dataset_dicts = {}

    for dataset_name,dataset_filename_list in settings['dataset_filenames'].items():

        dataset_generators[dataset_name] = []
        dataset_dicts[dataset_name] = []

        for dataset_filename in dataset_filename_list: 
            # begin processing sequencing data
            seq_data = create_data_generator(dataset_filename)
            seq_dict = frequency_dictionary(seq_data,**settings)
            adjusted_seq_dict = adjust_frequency_dictionary(seq_dict,**settings)
            # store data
            dataset_generators[dataset_name].append(seq_data)
            dataset_dicts[dataset_name].append(adjusted_seq_dict)

    data_dict = {
            'DATA':dataset_dicts,
            'SETTINGS':settings,
            }

    # Publish data
    print ('Publishing data...')
    publish_datasets_to_excel(data_dict,directory)
    print ('Finished!')

    print ('Job Finished!')

    return dataset_dicts

def _get_data(wb):

    loaded_data = {}
    # pull single sheet from xlsx 
    ws = wb['DATA']

    # extract settings from file
    for i,row in enumerate(ws): 

        values = [cell.value for cell in row]

        if i == 0:
            # acquire necessary indices
            name_ind = values.index('Dataset')
            round_ind = values.index('Round')
            bc_ind = values.index('BC')
            count_ind = values.index('Count')
            continue

        name,rnd,bc,count = [values[i] for i in (name_ind,round_ind,bc_ind,count_ind)]

        try:
            loaded_data[name][rnd][bc] = count

        except KeyError:
            loaded_data[name] = [{} for _ in range(rnd+1)]
            loaded_data[name][rnd][bc] = count

        except IndexError:
            loaded_data[name] += [{} for _ in range(rnd+1-len(loaded_data[name]))]
            loaded_data[name][rnd][bc] = count


    return loaded_data









