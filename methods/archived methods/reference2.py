# standard libraries
import pickle

# nonstandard libraries
import openpyxl

# homegrown libraries
from methods.loader import create_data_generator,filter_sequence
from methods.utilities import hamming_dist,levenshtein_dist
from methods.output2 import publish_reference_to_excel

    
#---------------------------------#

def domain_reference_dictionary(fname):
    """ Creates domain reference dictionary """
    # 
    domain2seq = {}
    wb = openpyxl.load_workbook(fname)
    ws = wb.active
    
    # iterate through rows of file
    for i,row in enumerate(ws):
        if i == 0: continue # skip header
        domain2seq[row[0].value] = row[1].value

    return domain2seq


#---------------------------------#

def sequence_reference_dictionary(fname):
    """ Iterates across data, outputs dictionary with frequencies """
    data_iter = create_data_generator(fname)

    seq_dict = {}

    for i,line in enumerate(data_iter):
        my_seq = filter_sequence(line)
        if my_seq:
            try:
                seq_dict[my_seq].append(line)
            except KeyError:
                seq_dict[my_seq] = [line]

    return seq_dict
            
#---------------------------------#

def merge_references(domain2seq,bc2seq,**kwargs):

    settings = {
            'dist_function':  'levenshtein',
            'dist_threshold':             1,
            }

    # use user settings
    if 'dist_function' in kwargs: 
        settings['dist_function'] = kwargs['dist_function']
    if 'dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['dist_threshold']
    
    # select distance function
    if settings['dist_function'] == 'levenshtein': dist_func = levenshtein_dist
    elif settings['dist_function'] == 'hamming':   dist_func = hamming_dist 
    else: raise KeyError('Distance function not recogonized! Try levenshtein or hamming.')
    
    threshold = settings['dist_threshold'] # get threshold

    reference = {} # initialize final object
    min_seq = max([len(seq) for seq in domain2seq.values()])

    for bc_count,(bc,seqs) in enumerate(bc2seq.items()):
            
        print('{}/{} sequences processed...'.format(bc_count,len(bc2seq)))

        '''# premature termination bridge
        if bc_count > 10: 
            break
        #'''#

        reference[bc] = []

        for seq_count,seq in enumerate(seqs): # iterate sequences
            
            matched_domains_with_location = []
            domain_list = []
            
            if seq_count >= 25: break
            front,domain_list = 0,[]
            for back in range(min_seq,len(seq)):
                if back - front < min_seq: continue # if segment smaller than smallest domain

                # get the domains that matched
                matched_domains = _endswith(seq[front:back],domain2seq,dist_func,threshold)

                # iterate through each matched domain and create an entry with location
                for (domain,length) in matched_domains:
                    matched_domains_with_location += [(domain,length,back - length,back)]

            # sort entries by start position, then by length (ascending, descending, respectively)
            matched_domains_with_location.sort(key = lambda x: (x[2],-x[1]))

            # take the largest domains that do not overlap and save in domain_list
            current_pos = -1
            for (domain,length,start,end) in matched_domains_with_location:
                if start > current_pos:
                    domain_list.append(domain)
                    current_pos = end
                    
            #print('matched domains:',matched_domains_with_location)
            #print('domain list:',domain_list)

            # if atleast one sequence was aligned
            if len(domain_list) > 0:
                reference[bc].append(tuple(domain_list))

        #print(bc,': ',reference[bc])

    return reference

#---------------------------------#

def _endswith(seq,domain_dict,dist_func,threshold):
    """ Checks for domain matching end of input sequence at certain hamming distance """
    matched_domains = []
    
    for domain,domain_seq in domain_dict.items():
        if len(seq) < len(domain_seq):
            continue
        if 1. - float(dist_func(seq[-len(domain_seq):],domain_seq))/len(domain_seq) >= threshold:
            matched_domains += [(domain,len(domain_seq))]
    
    return matched_domains
