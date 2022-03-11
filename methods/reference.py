# standard libraries
import pickle

# nonstandard libraries
import openpyxl

# homegrown libraries
from methods.loader import get_generator, get_generator2, create_data_generator,filter_sequence, filter_sequence2
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
        domain2seq[row[0].value] = row[2].value

    return domain2seq


#---------------------------------#

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in dna[::-1])
    return reverse_complement

def sequence_reference_dictionary(fname,pre,post):
    """ Iterates across data, outputs dictionary with frequencies """
    #data_iter = create_data_generator(fname)
    data_iter = get_generator(fname)

    seq_dict = {}

    for i,line in enumerate(data_iter):
        
        #print('Line {}: {}'.format(i,line))
        #line = line.replace("\n", "")
        
        #print('DNA?: {}'.format(all([i in 'ATGC' for i in line]))
        
        #if all([i in 'ATGC' for i in line.replace("\n", "")]): # only process DNA encoding lines
                     
        my_seq = filter_sequence(line,pre,post)

        if my_seq:
            #print('Barcode found: {}'.format(my_seq))
            try:
                seq_dict[my_seq].append(line)
            except KeyError:
                seq_dict[my_seq] = [line]
        else: # try reverse complement
            #line = str(Seq(line).reverse_complement())
            line = reverse_complement(line)
            my_seq = filter_sequence(line,pre,post)
            if my_seq:
                #print('Barcode found in reverse complement: {}'.format(my_seq))
                try:
                    seq_dict[my_seq].append(line)
                except KeyError:
                    seq_dict[my_seq] = [line]

    return seq_dict

def sequence_reference_dictionary2(fname,pre,post,start,stop):
    """ Iterates across data, outputs dictionary with frequencies and quality score """
    #data_iter = create_data_generator(fname)
    data_iter = get_generator2(fname)

    seq_dict = {}
    phred = {chr(x+33).encode('ascii'):x for x in range(0,94)}


    for i,(line,qseq) in enumerate(data_iter):
        
        #print('Line {}: {}'.format(i,line))
        #line = line.replace("\n", "")
        
        #print('DNA?: {}'.format(all([i in 'ATGC' for i in line]))
        
        #if all([i in 'ATGC' for i in line.replace("\n", "")]): # only process DNA encoding lines
                     
        my_seq,line = filter_sequence2(line,pre,post,start,stop)
        score = sum([phred[k] for k in qseq])/len(qseq)

        if my_seq:
            #print('Barcode found: {}'.format(my_seq))
            try:
                seq_dict[my_seq].append([line,score,len(line)])
            except KeyError:
                seq_dict[my_seq] = [line,score,len(line)]
        else: # try reverse complement
            #line = str(Seq(line).reverse_complement())
            line = reverse_complement(line)
            my_seq = filter_sequence(line,pre,post)
            if my_seq:
                #print('Barcode found: {}'.format(my_seq))
                try:
                    seq_dict[my_seq].append([line,score,len(line)])
                except KeyError:
                    seq_dict[my_seq] = [line,score,len(line)]
                    
    # sort by descending quality score, then descending length
    for k in seq_dict.keys():
        seq_dict[k] = seq_dict[k].sort(key = lambda x: (-x[1],-x[2]))
        
    # take first value (sequences) only, in sorted order
    sorted_seqs = {}
    for k,vals in seq_dict.items():
        for v in vals:
            try:
                seqs[k].append(v[0])
            except KeyError:
                seqs[k] = v[0]

    return sorted_seqs
            
#---------------------------------#

def merge_references(domain2seq,bc2seq,**kwargs):

    settings = {
            'dist_function':      'hamming',
            'dist_threshold':           0.8,
            }

    # use user settings
    if 'dist_function' in kwargs: 
        settings['dist_function'] = kwargs['dist_function']
    if 'dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['dist_threshold']
    
    # select distance function
    if settings['dist_function'] == 'levenshtein': dist_func = levenshtein_dist
    elif settings['dist_function'] == 'hamming':   dist_func = hamming_dist 
    else: raise KeyError('Distance function not recognized! Try levenshtein or hamming.')
    
    threshold = settings['dist_threshold'] # get threshold

    reference = {} # initialize final object
    min_seq = min([len(seq) for seq in domain2seq.values()])

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

def merge_references2(domain2seq,bc2seq,**kwargs):

    settings = {
            'dist_function':      'hamming',
            'dist_threshold':           0.8,
            }

    # use user settings
    if 'dist_function' in kwargs: 
        settings['dist_function'] = kwargs['dist_function']
    if 'dist_threshold' in kwargs: 
        settings['dist_threshold'] = kwargs['dist_threshold']
    
    # select distance function
    if settings['dist_function'] == 'levenshtein': dist_func = levenshtein_dist
    elif settings['dist_function'] == 'hamming':   dist_func = hamming_dist 
    else: raise KeyError('Distance function not recognized! Try levenshtein or hamming.')
    
    threshold = settings['dist_threshold'] # get threshold

    reference = {} # initialize final object
    min_seq = min([len(seq) for seq in domain2seq.values()])

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
            
            if seq_count >= 20: break
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
