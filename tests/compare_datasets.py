
# standard libraries

# nonstandard libraries

# homegrown libraries
#from methods.visuals import *
from methods.enrichment import *

def run(datasets_dict,**kwargs):

    directory = './datasets/'
    overwrite = False

    # datasets submited as dictioanry of lists
    dataset_filenames = {
            'reference':None,
            'enrichment':True,
            'counts':True,
            'targets':('BC','Domains')
            }

    settings = {
            }

    # update dictionary
    if 'reference' in kwargs: 
        settings['reference'] = kwargs['reference']
    if 'overwrite' in kwargs: 
        overwrite = kwargs['overwrite']

    wb = openpyxl.Workbook() # open new workbook

    write_counts(datasets_dict,workbook=wb)
    write_counts(datasets_dict,workbook=wb)

    try: del wb['Sheet']
    except KeyError: pass

    wb.save('results.xlsx')


             
