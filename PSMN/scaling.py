#! /usr/bin/env python3
"""
Extract and plot the speed up
"""
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def extract_simulation_time(path):
    """Extract the wall time that one simulation takes

    Parameters
    ----------
    path : str
        Path to the production log file

    Returns
    -------
    time : float
        Wall time of the simulation
    """
    command = "grep Time: " + path
    command = command.split()
    out = str(subprocess.check_output(command))
    out = out.replace("b'", "").replace("Time:", "").replace("\\n'", "")
    out = out.split()
    time = float(out[1])
    return time


def get_time(path):
    """Get the wall time of multiple simulation

    Parameters
    ----------
    path : str
        Path were the simulations have been run

    Returns
    -------
    cores_list : list
        List containing the number of cores used
    time_list : list
        List containing the corresponding wall time
    """
    # Init
    cores_list = []
    time_list = []
    # Search inside directory
    for name in os.listdir(path):
        # Take only directories with cores in their name
        if "cores_" in name:
            cores = int(name.replace("cores_", ""))
            # Get the number of cores
            cores_list.append(cores)
            # Production log file
            log_path = path + name + "/production/prod.log"
            # Get the wall time
            time = extract_simulation_time(log_path)
            time_list.append(time)
    return cores_list, time_list


def speedup(cores_list, time_list):
    """Compute the speed-up

    Parameters
    ----------
    cores_list : list
        List containing the number of cores used
    time_list : list
        List containing the corresponding wall time

    Returns
    -------
    cores_list : np.array
        Array containing the number of cores used
    speedup : np.array
        Array containing the corresponding speedup
    """
    # Convert to array
    cores_list = np.array(cores_list)
    time_list = np.array(time_list)
    minimal_core_index = np.argmin(cores_list)
    # Compute speed-up
    speedup = time_list[minimal_core_index] / time_list
    return cores_list, speedup


def plot_speedup(cores_list, speedup, color="b", outputname="speedup.pdf"):
    """Plot the speed up as a fonction of the number of nodes

    Parameters
    ----------
    cores_list : np.array
        Array containing the number of cores used
    speedup : np.array
        Array containing the corresponding speedup
    color : str, optional
        Color of the points, by default 'b'
    outputname : str, optional
        Path where to save the file, by default 'speedup.pdf'
    """
    fig, ax = plt.subplots()
    # Plot a theoritical linear curve
    linear = cores_list / np.amin(cores_list)
    ax.plot(cores_list, linear, "--", label="Linear law")
    ax.scatter(cores_list, speedup, color=color)
    ax.set_xticks(cores_list)
    ax.set_xlabel("Number of CPU")
    ax.set_ylabel("Speed-up")
    ax.legend()
    ax.grid()
    plt.savefig(outputname, format="pdf", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    color = sns.color_palette("cool", 12)[6]
    path = "/home/ccattin/Documents/HSP90_278K/test_parallel/"
    outputname = "/home/ccattin/Documents/HSP90_278K/test_parallel/analysis/speedup.pdf"
    cores_list, time_list = get_time(path)
    cores_list, speedup = speedup(cores_list, time_list)
    plot_speedup(cores_list, speedup, color, outputname)
