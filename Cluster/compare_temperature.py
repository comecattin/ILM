#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

def normalize(cluster_number):
    cluster_number_normalized = cluster_number/np.sum(cluster_number,axis=1).reshape((2,1))
    return cluster_number_normalized

def plot_barplot(cluster_number):
    fig, ax = plt.subplots()
    color_palette = sns.color_palette('cool',12)
    color = [color_palette[6],color_palette[2]]
    cluster_id = np.arange(5)
    temperatures=(278,300)

    for i, number in enumerate(cluster_number):
        ax.bar(cluster_id + i*0.25,
               number,
               width=0.25,
               color=color[i],
               label=f"{temperatures[i]}K")

    ax.set_xlabel("Cluster number")
    ax.set_ylabel("Normalized population")
    ax.set_xticks(cluster_id+0.125)
    ax.set_xticklabels(cluster_id+1)
    ax.legend()


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
    cluster_number_normalized = normalize(cluster_number=cluster_number)
    plot_barplot(cluster_number=cluster_number_normalized)