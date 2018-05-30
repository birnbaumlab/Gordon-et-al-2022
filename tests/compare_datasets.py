
# standard libraries

# nonstandard libraries

# homegrown libraries
#from methods.visuals import *
from methods.enrichment import *

def run(datasets_dict,**kwargs):

    overwrite = False

    # datasets submited as dictioanry of lists
    settings = {
            'reference':None,
            #'enrichment':True,
            #'counts':True,
            #'targets':('BC','Domains')
            }

    # update dictionary
    if 'reference' in kwargs: 
        settings['reference'] = kwargs['reference']

    # local namespace
    reference = settings['reference']

    wb = openpyxl.Workbook() # open new workbook

    comparison = write_counts(datasets_dict,workbook=wb)

    if reference:
        comparison = write_counts(datasets_dict,workbook=wb,reference=reference)

    try: del wb['Sheet']
    except KeyError: pass

    wb.save('results.xlsx')

    return comparison


             
