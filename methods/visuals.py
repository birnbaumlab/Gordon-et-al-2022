
# standard libraries
from random import shuffle

# nonstandard libraries
import matplotlib.pyplot as plt


def show_saturation(seq_dicts):
    """ Show scatters for uniques """

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

    plt.xlabel('Read count')
    plt.ylabel('Unique read count')

    plt.legend()
    plt.show(block=False)

    input('Press enter to close...')
