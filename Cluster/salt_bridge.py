#! /usr/bin/env python3
"""Compare the salt bridges present in each cluster made by TTClust"""

import os
import numpy as np



def read_out(file):

    with open(file) as f:
        data = f.readlines()
    data = data[18:-2]
    data = [sb.replace('\n','') for sb in data]
    return data

def compare_interaction(directory):
    list_dir = os.listdir(directory)
    list_dir.sort()
    table = []
    for i, name in enumerate(list_dir):
        table.append(
            read_out(
                os.path.join(directory,name,f'sb_c{i+1}')
            )
        )
    
    compare = []
    for cluster_1 in table:
        for cluster_2 in table:
            compare.append(
                list(
                    set(cluster_1) & set(cluster_2)
                )
            )
    
    compare = np.array(compare,dtype=object)
    compare = compare.reshape((len(table),len(table)))
    
    return compare

if __name__ == '__main__':
    directory='/home/ccattin/Documents/Cluster/traj15/all_dt2000/clustering/analysis/'
    file='C1/sb_c1'
    data = read_out(directory+file)
    compare = compare_interaction(directory)