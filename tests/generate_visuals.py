
# standard libraries
import os

# nonstandard libraries

# homegrown libraries
from methods.visuals import *

def run(**kwargs):

    directory = './datasets/'
    overwrite = False

    # datasets submited as dictioanry of lists
    settings = {
            'datasets':None,
            'reference':None,
            'comparison':None,
            'directory':'results',
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
    directory = settings['directory']

    # change in directory
    home = os.getcwd() 
    os.chdir(os.path.join(home,directory))

    for i in range(1,10000):

        results_folder = os.path.join(home,directory,str(i).zfill(4))

        if os.path.isdir(results_folder): 
            continue
        else:
            os.makedirs(results_folder)
            os.chdir(results_folder)
            break

    print('Saving figure in folder: {}...'.format(results_folder))

    # generate figure batches for each folder
    if datasets:
        datasets_visuals(datasets) 

    if reference:
        reference_visuals(reference) 

    if comparison:
        comparison_visuals(comparison) 











