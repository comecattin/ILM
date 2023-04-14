#! /usr/bin/env python3
"""Extract and plot the relative population for each cluster at different temperature"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_log(log_file):
    """Extract the frame number in each cluster

    Parameters
    ----------
    log_file : str
        Path to the log file of TTClust (clustering.log usually)

    Returns
    -------
    cluster : list
        List containing all the frame in every cluster. A sublist is associated with one cluster
    """
    # Init
    clusters = []

    # Open and read
    with open(log_file) as f:
        for line in f.readlines():
            if 'Frames' in line:
                # Remove unecessary stuff
                list_frame = line.replace('\n','').replace(' ','').replace('Frames','').replace(':','')
                # Convert into list and append
                list_frame = [int(s) for s in list_frame.split(',')[1:-1]]
                clusters.append(list_frame)
    
    return clusters

def get_dict(clusters, temperatures, trajectory_lim):
    """Convert the list of frame into a dictionary to be clearer.

    Parameters
    ----------
    cluster : list
        List containing all the frame in every cluster. A sublist is associated with one cluster
    temperatures : tuple
        Tuple containing the different temperature
    trajectory_lim : int
        Number of frame of trajectory

    Returns
    -------
    dict_cluster : dict
        Dictionary containing all the frame at a certain temperature in a certain cluster. The dictionary keys are in format (id_cluster, temperature)
    """

    # Init
    dict_cluster = {}

    # Loop over all the temperature
    for temperature in temperatures:
        # Loop over all the cluster
        for id_cluster, cluster in enumerate(clusters):
            # Get the upper part of the cluster or the lower part
            if temperature == temperatures[0]:
                dict_cluster[(id_cluster,temperature)] = [frame for frame in cluster if frame <= trajectory_lim]
            if temperature == temperatures[1]:
                dict_cluster[(id_cluster,temperature)] = [frame for frame in cluster if frame >= trajectory_lim]
    
    return dict_cluster

def get_number_frames(dict_cluster):
    """Get the number of frame in each cluster

    Parameters
    ----------
    dict_cluster : dict
        Dictionary containing all the frame at a certain temperature in a certain cluster. The dictionary keys are in format (id_cluster, temperature)

    Returns
    -------
    cluster_number : np.array
        Array containing for each cluster and temperature the number of frame inside. A row represent a temperature. Columns are clusters
    """
    # Init
    cluster_number = []
    # Loop over all the cluster and temperature
    for key, list_of_frame in dict_cluster.items():
        cluster_number.append(len(list_of_frame))
    # Convert and reshape in the correct shape
    cluster_number = np.array(cluster_number).reshape((2,5))
    
    return cluster_number

def normalize(cluster_number):
    """Get the relative population in each cluster

    Parameters
    ----------
    cluster_number : np.array
        Array containing for each cluster and temperature the number of frame inside. A row represent a temperature. Columns are clusters

    Returns
    -------
    cluster_number_normalized : np.array
        Array containing the relative population of each cluster at different temperature.
    """
    cluster_number_normalized = cluster_number/np.sum(cluster_number,axis=1).reshape((2,1))
    
    return cluster_number_normalized

def plot_barplot(cluster_number,temperatures, output):
    """Plot the result as barplot

    Parameters
    ----------
    cluster_number : np.array
        Array containing the relative population of each cluster at different temperature.
    temperatures : tuple
        Tuple containing the different temperature
    output : str
        Path to the .pdf file to output
    """
    # Init
    fig, ax = plt.subplots()
    # Color, cluster numbers and temperatures
    color_palette = sns.color_palette('cool',12)
    color = [color_palette[6],color_palette[2]]
    cluster_id = np.arange(5)

    # Loop over the different temperatures
    for i, number in enumerate(cluster_number):
        ax.bar(cluster_id + i*0.25,
               number,
               width=0.25,
               color=color[i],
               label=f"{temperatures[i]}K")

    # Aesthetic
    ax.set_xlabel("Cluster number")
    ax.set_ylabel("Normalized population")
    ax.set_xticks(cluster_id+0.125)
    ax.set_xticklabels(cluster_id+1)
    ax.legend()

    # Save
    plt.savefig(output, dpi=300, bbox_inches="tight")
    # Display
    plt.show()

if __name__ == '__main__':
    log_file = '/home/ccattin/Documents/Cluster/total/clustering/clustering.log'
    output = '/home/ccattin/Documents/Code/outputs/clustering_temperature.pdf'
    trajectory_lim = 4001
    temperatures = (278,300)

    #Get the cluster from thye log
    clusters = load_log(log_file)

    #Convert to dict
    dict_cluster = get_dict(clusters=clusters,
                            temperatures=temperatures, 
                            trajectory_lim=trajectory_lim)
    
    # Get the population of each cluster
    cluster_number = get_number_frames(dict_cluster=dict_cluster)
    
    # Get the normalized population
    cluster_number_normalized = normalize(cluster_number=cluster_number)
    
    # Plot and save the result
    plot_barplot(cluster_number=cluster_number_normalized,
                 temperatures=temperatures,
                 output=output)