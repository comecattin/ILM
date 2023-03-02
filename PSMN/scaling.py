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


def get_CPU_node_time(path):
    """Get the number of cpu, of node and their associated job time

    Parameters
    ----------
    path : str
        Path were the analysis simulations are

    Returns
    -------
    cpu : np.array
        List of the number of CPU
    node : np.array
        List of the number of nodes
    time : np.array
        List of the associated time
    """
    # Initialize
    cpu = []
    node = []
    time = []
    # Search inside directory
    for name in os.listdir(path):
        # Get only the right format
        if "CPU" in name:
            # Extract the cpu and the node
            list_CPU_NODE = (
                name.replace("CPU_", "")
                .replace("NODE_", "")
                .replace(".log", "")
                .split("_")
            )
            cpu.append(int(list_CPU_NODE[0]))
            node.append(int(list_CPU_NODE[1]))
            # Extract the time
            time.append(extract_simulation_time(path + "/" + name))
    return np.array(cpu), np.array(node), np.array(time)


def plot_CPU_node_time(
    cpu, node, time, color=["b", "r", "g"], outputname="speedup_node.pdf"
):
    """Plot the speedup as a function of node number and cpu number

    Parameters
    ----------
    cpu : np.array
        List of the number of CPU
    node : np.array
        List of the number of nodes
    time : np.array
        List of the associated time
    color : list, optional
        Color of the plot, by default ["b","r","g"]
    outputname : str, optional
        Output filename, by default "speedup_node.pdf"
    """
    fig, ax = plt.subplots()
    # Loop over the nodes
    for i, node_i in enumerate((1, 2, 4)):
        # Get the associated cpu and time
        cpu_node_i = list(cpu[np.where(node == node_i)])
        time_node_i = list(time[np.where(node == node_i)])
        # Get the reference for the speedup
        #   The reference is the time for
        #       the minimal number of node and cpu
        if node_i == 1:
            minimal_cpu_index = np.argmin(cpu_node_i)
            ref = time_node_i[minimal_cpu_index]
        # Compute the speedup
        speedup = ref / time_node_i
        ax.scatter(cpu_node_i, speedup, label=str(node_i) + " nodes", color=color[i])
    ax.legend()
    ax.grid()
    ax.set_xticks([2, 4, 8, 16, 32, 64, 96, 128])
    ax.set_xlabel("Number of CPUs")
    ax.set_ylabel("Speed-up")
    plt.savefig(outputname, format="pdf", dpi=300, bbox_inches="tight")
    plt.show()

def get_CPU_MPI_OMP_time(path,
                         outputname="MPI_OMP.pdf",
                         color=["b","r"]
                         ):
    """Get and plot the time as a function of MPI threads for various cores

    Parameters
    ----------
    path : str
        Path were the analysis simulations are
    outputname : str, optional
        Output file name, by default "MPI_OMP.pdf"
    color : list, optional
        Color of the plot, by default ["b","r"]
    """
    # Initialize
    cpu = []
    i = 0
    fig,ax = plt.subplots()
    # Search inside directory
    for name in os.listdir(path):
        # Get only the right format
        joined = os.path.join(path,name)
        #Get the CPU
        if os.path.isdir(joined):
            cpu_i = int(name.replace("_NODE_1",""))
            cpu.append(cpu_i)
            i = len(cpu) - 1
            #Get the MPI and OMP
            MPI = []
            OMP = []
            time = []
            for file in os.listdir(joined):
                #MPI = 0 represent the default selected by gromacs
                if file == "ref.log":
                    MPI.append(0)
                else:
                    MPI_OMP = file.replace("MPI_","").replace("OMP_","").replace(".log","").split("_")
                    MPI.append(int(MPI_OMP[0]))
                    OMP.append(int(MPI_OMP[1]))
                #Extract the time
                time.append(
                    extract_simulation_time(
                    os.path.join(joined,file)
                    )
                )
            ax.scatter(
                MPI,
                time,
                label="CPU= "+str(cpu_i),
                color=color[i])
    ax.legend()
    ax.set_xticks((0,1,2,3,4))
    ax.grid()
    ax.set_xlabel("Number of MPI threads")
    ax.set_ylabel("Production time (s)")
    plt.savefig(outputname, format="pdf", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    ###     One node    ###
    color_palette = sns.color_palette("cool", 12)
    path = "/home/ccattin/Documents/HSP90_278K/test_parallel/"
    outputname = "/home/ccattin/Documents/HSP90_278K/test_parallel/analysis/speedup.pdf"
    cores_list, time_list = get_time(path)
    cores_list, speedup_list = speedup(cores_list, time_list)
    #plot_speedup(cores_list, speedup_list, color_palette[6], outputname)

    ###     Various node    ###
    path_node = "/home/ccattin/Documents/HSP90_278K/test_parallel/analysis"
    outputname = (
        "/home/ccattin/Documents/HSP90_278K/test_parallel/analysis/speedup_node.pdf"
    )
    cpu, node, time = get_CPU_node_time(path_node)
    color = [color_palette[i] for i in (0, 6, 11)]
    #plot_CPU_node_time(cpu, node, time, color=color, outputname=outputname)
    
    ###     Same node MPI OMP different     ###
    outputname = (
        "/home/ccattin/Documents/HSP90_278K/test_parallel/analysis/MPI_OMP.pdf"
    )
    color = [color_palette[i] for i in (0, 6)]
    get_CPU_MPI_OMP_time(
        path_node,
        outputname=outputname,
        color=color
        )
