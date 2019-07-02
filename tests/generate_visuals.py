
# standard libraries
import os

# nonstandard libraries

# homegrown libraries
from methods.visuals import *

def run(**kwargs):

    # datasets submited as dictioanry of lists
    settings = {
            'datasets':None,
            'reference':None,
            'comparison':None,
            }

    # update dictionary
    if 'datasets' in kwargs: 
        settings['datasets'] = kwargs['datasets']
    if 'reference' in kwargs: 
        settings['reference'] = kwargs['reference']
    if 'comparison' in kwargs: 
        settings['comparison'] = kwargs['comparison']
    if 'results' in kwargs: 
        settings['results'] = kwargs['results']
    if 'overwrite' in kwargs: 
        overwrite = kwargs['overwrite']

    # local namespace
    datasets = settings['datasets']
    reference = settings['reference']
    comparison = settings['comparison']

    # change directory
    os.mkdir(settings['results_directory'])
    os.chdir(settings['results_directory'])

    print('Saving figures in folder: {}...'.format(results_folder))

    # generate figure batches for each folder
    if datasets:
        datasets_visuals(datasets) 

    if reference:
        reference_visuals(reference) 

    if comparison:
        comparison_visuals(comparison) 











