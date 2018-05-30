
# nonstandard libraries
import numpy as np

# homegrown libraries
from methods.utilities import hamming_dist_np

#---------------------------------#

def create_data_generator(fname):
    """ Opens data as a generator (i.e. don't load data to RAM) """
    with open(fname,'r') as infile:
        for line in infile:
            yield line

#---------------------------------#

def filter_sequence(mystr,pre=None,post=None):
    """ Returns filtered string, otherwise False """
    if pre == None:
        pre  = 'CCCTTGGAGAACCACCTTGTTGG'
    if post == None:
        post = 'GTTTAAGAGCTAAGCTGGAAAAAGTGGCAC'
    if pre in mystr and post in mystr:
        return mystr[mystr.index(pre)+len(pre):mystr.index(post)]
        
#---------------------------------#

def frequency_dictionary(data_iter,**kwargs):
    """ Iterates across data, outputs dictionary with frequencies """

    settings = {
            'pre':None,
            'post':None,
            }

    if 'barcode_5p' in kwargs: 
        settings['pre'] = kwargs['barcode_5p']
    if 'barcode_3p' in kwargs: 
        settings['post'] = kwargs['barcode_3p']

    pre,post = settings['pre'],settings['post']
        
    freq_dict = {}

    for i,line in enumerate(data_iter):
        my_seq = filter_sequence(line,pre=pre,post=post)
        if my_seq:
            try:
                freq_dict[my_seq] += 1
            except KeyError:
                freq_dict[my_seq] =  1

        if i != 0 and i % 500000 == 0:
            print ('{} lines processed...'.format(i))

    return freq_dict
    
#---------------------------------#

def adjust_frequency_dictionary(freq_dict,**kwargs):

    # default settings
    characters = 'ATCG'
    settings = {
            'clustering_threshold': 1,
            'count_threshold': 1,
            'silent':False,
            }

    # adjust filtering settings
    if 'clustering_threshold' in kwargs: 
        settings['clustering_threshold'] = kwargs['clustering_threshold']
    if 'count_threshold' in kwargs: 
        settings['count_threshold'] = kwargs['count_threshold']
    if 'silent' in kwargs: 
        settings['silent'] = kwargs['silent']
        
    # unload in local namespace
    cluster_threshold = settings['clustering_threshold']
    count_threshold = settings['count_threshold']
    silent = settings['silent']

    freq_items = [(k,v) for k,v in freq_dict.items() if v >= count_threshold]

    print('Count threshold ({}) reduced barcodes from {} -> {}'.format(
        count_threshold,len(freq_dict),len(freq_items)))
    
    freq_items = sorted(freq_items,key=lambda x: -x[1])

    # convert frequency items to minimum processible form
    seq2array = dict([(seq[0],np.array([characters.index(s) 
        for s in seq[0]])) for seq in freq_items])

    filtered_freq_dict = {} 
    
    while len(freq_items) > 0:

        if len(freq_items)%1000 == 0: print('{} sequences remaining, {} accepted...'.format(
            len(freq_items),len(filtered_freq_dict)))

        item = freq_items[0]
        exit = False

        for j,item2 in enumerate(filtered_freq_dict):
            
            if hamming_dist_np(seq2array[item[0]],seq2array[item2]) <= cluster_threshold:
                filtered_freq_dict[item2] += item[1]
                exit = True 
                break

        if exit == False:
            filtered_freq_dict[item[0]] = item[1]

        freq_items.pop(0)

    if not silent:
        print('Filtered barcode count: {} -> {}...'.format(
              len(freq_dict),len(filtered_freq_dict)))

    return filtered_freq_dict         
    
#---------------------------------#

def get_xlsx_settings(wb):
    """ Extract attributes from settings sheet """
    loaded_settings = {}
    # pull single sheet from xlsx 
    ws = wb['SETTINGS']
    # extract settings from file
    for row in ws: 

        values = [v.value for v in row]
       
        if values[1] == 'dict->':

            subdict_key = values[0] 
            loaded_settings[subdict_key] = {} 

        elif values[0] == '->':

            loaded_settings[subdict_key][values[1]] = [v for v in values[2:] if v != None]
        
        else:
        
            loaded_settings[values[0]] = values[1] 

    return loaded_settings
















