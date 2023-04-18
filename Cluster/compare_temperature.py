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
            if "Frames" in line:
                # Remove unecessary stuff
                list_frame = (
                    line.replace("\n", "")
                    .replace(" ", "")
                    .replace("Frames", "")
                    .replace(":", "")
                )
                # Convert into list and append
                list_frame = [int(s) for s in list_frame.split(",")[1:-1]]
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
    for i_temp, temperature in enumerate(temperatures):
        # Loop over all the cluster
        for id_cluster, cluster in enumerate(clusters):
            # Get the upper part of the cluster or the lower part
            lower, higher = trajectory_lim[i_temp*2:(i_temp+1)*2]
            # Append to the dict
            dict_cluster[(id_cluster, temperature)] = [
                frame for frame in cluster if (frame > lower and frame <= higher)
            ]

    return dict_cluster


def get_number_frames(dict_cluster, number_temperature, number_cluster):
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
    cluster_number = np.array(cluster_number).reshape(
        (number_temperature, number_cluster)
        )

    return cluster_number


def normalize(cluster_number, number_temperature):
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
    cluster_number_normalized = cluster_number / np.sum(cluster_number, axis=1).reshape(
        (number_temperature, 1)
    )

    return cluster_number_normalized


def plot_barplot(cluster_number, temperatures, output, color):
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
    cluster_id = np.arange(len(cluster_number[0]))

    # Loop over the different temperatures
    for i, number in enumerate(cluster_number):
        ax.bar(
            cluster_id + i * 0.25,
            number,
            width=0.25,
            color=color[i],
            label=temperatures[i],
        )

    # Aesthetic
    ax.set_xlabel("Cluster number")
    ax.set_ylabel("Normalized population")
    ax.set_xticks(cluster_id + 0.125)
    ax.set_xticklabels(cluster_id + 1)
    ax.legend()

    # Save
    plt.savefig(output, dpi=300, bbox_inches="tight")
    # Display
    plt.show()


if __name__ == "__main__":

    log_file = "/home/ccattin/Documents/Cluster/total_and_data/clustering/clustering.log"
    output = "/home/ccattin/Documents/Code/outputs/clustering_temperature.pdf"
    
    trajectory_lim = (0,4001,4001,8021,8021,12022)
    
    temperatures = ("278K", "300 Amber 14", "300K Amber 19")
    
    number_temperature = len(temperatures)
    number_cluster = 5
    
    color_palette = sns.color_palette("cool", 12)
    color = [color_palette[6], color_palette[2],color_palette[10]]
    
    # Get the cluster from the log
    clusters = load_log(log_file)

    # Convert to dict
    dict_cluster = get_dict(
        clusters=clusters, temperatures=temperatures, trajectory_lim=trajectory_lim
    )

    # Get the population of each cluster
    cluster_number = get_number_frames(dict_cluster=dict_cluster,
                                       number_temperature=number_temperature,
                                       number_cluster=number_cluster)

    # Get the normalized population
    cluster_number_normalized = normalize(cluster_number=cluster_number,
                                          number_temperature=number_temperature)

    # Plot and save the result
    plot_barplot(
        cluster_number=cluster_number_normalized,
        temperatures=temperatures,
        output=output,
        color=color
    )
