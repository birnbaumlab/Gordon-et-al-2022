
from tests import load_dataset 
from tests import create_reference 
from tests import compare_datasets 

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
        'dataset_filenames':dataset_filenames,
        'count_threshold':10,
        }

datasets = load_dataset.run(**settings)
reference = create_reference.run()

comparison = compare_datasets.run(datasets,reference=reference)

