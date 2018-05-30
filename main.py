
from tests import load_dataset
from tests import create_reference 
from tests import compare_datasets 
from tests import generate_visuals 

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
        'dataset_filenames':                      dataset_filenames,
        'count_threshold':                                       10,
        # reference settings
        'sequence_reference_filename':             'CCS_Run2.fastq', 
        'domain_reference_filename': 'ICD lists and sequences.xlsx', 
        'dist_function':                                  'hamming',
        'dist_threshold':                                       0.8,
        # comparison settings
        # visual settings
        }

datasets = load_dataset.run(**settings)
reference = create_reference.run(**settings)
comparison = compare_datasets.run(datasets,reference=reference,**settings)

generate_visuals.run(datasets=datasets,reference=reference,comparison=comparison,**settings)

print('Finished!')


