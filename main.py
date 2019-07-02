
# standard libraries
import os

# custom libraries
from tests import load_dataset
from tests import create_reference 
from tests import compare_datasets 
from tests import generate_visuals 


def main():

    """ Main overhead function to run pipeline """

    """ --------------------- """ 
    """ MODIFY ANYTHING BELOW """
    """ --------------------- """ 

    base = '/run/media/pholec/Seagate Backup Plus Drive/Sequence_Files_PVH/'

    dataset_filenames = {
            'CD69+PD-1-': [base + fname for fname in (
                'Non-Selected.fastq',
                '1st CD69+ sorted.fastq',
                '2nd CD69+PD-1- sorted.fastq')],
            'CD69+PD': [base + fname for fname in (
                'Non-Selected.fastq',
                '1st CD69+ sorted.fastq',
                '2nd CD69+ sorted.fastq',
                '3rd CD69+ sorted.fastq')]
            }

    settings = {
            # general settings
            'overwrite':                                          False,
            'silent':                                             False,
            # dataset settings
            'clustering_threshold':                                   3,
            'count_threshold':                                       10,
            'dataset_filenames':                      dataset_filenames,
            'barcode_5p':                     'CCCTTGGAGAACCACCTTGTTGG',
            'barcode_3p':              'GTTTAAGAGCTAAGCTGGAAAAAGTGGCAC',
            # reference settings
            'sequence_reference_filename':             'CCS_Run2.fastq', 
            'domain_reference_filename': 'ICD lists and sequences.xlsx', 
            'dist_function':                                  'hamming',
            'dist_threshold':                                       0.8,
            # comparison settings
            # visual settings
            }

    """ --------------------- """ 
    """ MODIFY ANYTHING ABOVE """
    """ --------------------- """ 

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

    '''
    print('Generating reference...')
    reference = create_reference.run(**settings)
    print('Finished!')

    print('Analyzing datasets using reference...')
    comparison = compare_datasets.run(datasets,reference=reference,**settings)
    print('Finished!')

    print('Generating visuals...')
    generate_visuals.run(datasets=datasets,reference=reference,comparison=comparison,**settings)
    print('Finished!')
    '''

    print('Finished!')

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
        else: break

    return results_folder

# script catch
if __name__ == "__main__":
    main()


