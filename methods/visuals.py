
# standard libraries
from random import shuffle
import collections

# nonstandard libraries
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def datasets_visuals(datasets):
    
    for name,dataset in datasets.items():
        _show_saturation(dataset,name=name)


def reference_visuals(reference):

    # create consensus
    bc2domain = [(k,collections.Counter(v).most_common(1)) 
                       for k,v in reference.items()]
    bc2domain = dict([(k,v[0][0]) for k,v in bc2domain if len(v) > 0])

    ### BARCODE DOMAIN COUNTS ###
    # number of functional reads
    domain_freq = collections.Counter([len(v.split(',')) for v in bc2domain.values()])
    x,y = domain_freq.keys(),domain_freq.values()
    plt.figure()
    plt.bar(x,y,align='center',alpha=0.5)

    # label
    plt.title('Barcode Reference Reads')
    plt.xlabel('Domain Count (#)')
    plt.ylabel('Occurances (#)')
    plt.savefig('domains_per_barcode.png')

    ### BARCODE READ COUNTS ###
    # number of functional reads
    domain_freq = collections.Counter([len(v) for v in reference.values()])
    x,y = domain_freq.keys(),domain_freq.values()
    fig = plt.figure()
    ax = plt.gca()
    plt.bar(x,y,align='center',alpha=0.5)

    # label
    plt.title('Reads per Barcode')
    plt.xlabel('Recognized Reads/Barcode (#)')
    plt.ylabel('Occurances (#)')
    plt.savefig('reads_per_barcode.png')


def comparison_visuals(comparison):


    ### FREQUENCY OF DOMAINS ###    
    for name,domain_sets in comparison.items():
        
        filtered_domain_sets = dict([(tuple(k.split(',')),v) 
            for k,v in domain_sets.items()])
        fig,ax = plt.subplots(1,1,figsize=(4,20))

        index,scores = _show_domain_frequency(filtered_domain_sets,name,ax)
    
        plt.savefig('{}_domain_frequency.png'.format(name))

        fig,ax = plt.subplots(1,3,figsize=(12,20))

        index,scores = [index],[scores]

        _show_domain_enrichment(scores,index,name,positional = False)

        plt.savefig('{}_domain_enrichment.png'.format(name))

        sitewise_index,sitewise_scores = [],[]

        for i in range(3):
            filtered_domain_sets = [((k.split(',')[i],),v) 
                for k,v in domain_sets.items() if len(k.split(',')) == 3]
            sum_freq = np.array([d[1] for d in filtered_domain_sets]).sum(axis=0)
            filtered_domain_sets = dict([(k,np.array(v)/sum_freq) for k,v in filtered_domain_sets])
            index,scores = _show_domain_frequency(
                    filtered_domain_sets,name+'- Site {}'.format(i+1),ax[i])
            sitewise_index.append(index) 
            sitewise_scores.append(scores) 

        fig.savefig('{}_positional_domain_frequency.png'.format(name))

        _show_domain_enrichment(sitewise_scores,sitewise_index,name,positional = True)

        plt.savefig('{}_positional_domain_enrichment.png'.format(name))

def _show_domain_enrichment(scores,indices,name,positional = False):

    round_total = scores[0].shape[1]
    set_total = len(scores)

    if round_total <= 1:
        print('{} round(s) of data passed for {}, skipping enrichment...'.format(round_total))
        return None

    fig,ax = plt.subplots(1,set_total,figsize=(4*set_total,20),squeeze=False)
    rounds = ['Round {}/{}'.format(i+1,i) for i in range(round_total-1)]

    for i,(score,index) in enumerate(zip(scores,indices)):

        enriched_score = np.log2(score[:,1:]/score[:,:-1])

        cax = ax[0][i].imshow(enriched_score)

        ax[0][i].set_xticks(np.arange(len(rounds)))
        ax[0][i].set_yticks(np.arange(len(index)))

        ax[0][i].set_xticklabels(rounds)
        ax[0][i].set_yticklabels(index.keys())

        plt.setp(ax[0][i].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        if positional == True: ax[0][i].set_title(name+'- Position {}'.format(i))
        elif positional == False: ax[0][i].set_title(name)

def _show_domain_frequency(domain_sets,name,ax):

    round_total = len(next(iter(domain_sets.values())))

    domain_index = dict([(d,i) for i,d in enumerate(sorted(list(set(
        [a for b in domain_sets.keys() for a in b]))))])
    domain_scores = np.zeros((len(domain_index),round_total))

    for domain_set,freqs in domain_sets.items():
        for domain in domain_set:
            domain_scores[domain_index[domain]][:] += freqs

    cax = ax.imshow(domain_scores)

    rounds = ['Round {}'.format(i) for i in range(round_total)]

    ax.set_xticks(np.arange(len(rounds)))
    ax.set_yticks(np.arange(len(domain_index)))

    ax.set_xticklabels(rounds)
    ax.set_yticklabels(domain_index.keys())

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title(name)

    ticks = [0,0.5,1]

    #cbar = fig.colorbar(cax,ticks=ticks)
    #cbar.ax.set_yticklabels(['0%','50%','100%'])
    return domain_index,domain_scores


def _show_saturation(seq_dicts,name='Sequence Data'):
    """ Show scatters for uniques """

    plt.figure()

    for i,seq_dict in enumerate(seq_dicts):


        uniques = [index for index,count in enumerate(seq_dict.values()) 
                   for _ in range(count)]

        shuffle(uniques)
        unique_counter,memory = [0],set()

        for unique in uniques: 

            if unique in memory: 
                unique_counter.append(unique_counter[-1]) 

            else:
                unique_counter.append(unique_counter[-1] + 1) 
                memory.add(unique)
        
        plt.plot(range(len(uniques)+1),unique_counter,label = 'Round {}'.format(i))

    plt.title(name)
    plt.xlabel('Read Count (#)')
    plt.ylabel('Unique Read Count (#)')

    plt.legend(loc=4)

    plt.savefig('{}_saturation.png'.format(name))
