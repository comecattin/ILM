#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

def load_log(log_file):
    clusters = []
    with open(log_file) as f:
        for line in f.readlines():
            if 'Frames' in line:
                list_frame = line.replace('\n','').replace(' ','').replace('Frames','').replace(':','')
                list_frame = [int(s) for s in list_frame.split(',')[1:-1]]
                clusters.append(list_frame)
    return clusters

def get_dict(clusters, temperatures, trajectory_lim):
    dict_cluster = {}
    for temperature in temperatures:
        for id_cluster, cluster in enumerate(clusters):
            if temperature == temperatures[0]:
                dict_cluster[(id_cluster,temperature)] = [frame for frame in cluster if frame <= trajectory_lim]
            if temperature == temperatures[1]:
                dict_cluster[(id_cluster,temperature)] = [frame for frame in cluster if frame >= trajectory_lim]
    return dict_cluster

def get_number_frames(dict_cluster):
    cluster_number = []
    for key, list_of_frame in dict_cluster.items():
        cluster_number.append(len(list_of_frame))
    cluster_number = np.array(cluster_number).reshape((2,5))
    return cluster_number


def plot_barplot(cluster_number):
    fig, ax = plt.subplots()
    cluster_id = np.arange(5)

    for i, number in enumerate(cluster_number):
        ax.bar(cluster_id + i*0.25,number,width=0.25)




    plt.show()

if __name__ == '__main__':
    log_file = '/home/ccattin/Documents/Cluster/total/clustering/clustering.log'
    trajectory_lim = 4001
    temperatures = (278,300)
    clusters = load_log(log_file)

    dict_cluster = get_dict(clusters=clusters,
                            temperatures=temperatures, 
                            trajectory_lim=trajectory_lim)
    
    cluster_number = get_number_frames(dict_cluster=dict_cluster)

    plot_barplot(cluster_number=cluster_number)