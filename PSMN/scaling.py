#! /usr/bin/env python3
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def extract_simulation_time(path):

    command = "grep Time: "+path
    command = command.split()
    out = str(subprocess.check_output(command))
    out = out.replace("b'","").replace("Time:",'').replace("\\n'",'')
    out = out.split()
    time = float(out[1])

    return time

def get_time(path):
    cores_list = []
    time_list = []
    for name in os.listdir(path):
        if 'cores_' in name:
            cores = int(name.replace("cores_", ""))
            cores_list.append(cores)
            log_path = path+name+'/production/prod.log'
            time = extract_simulation_time(log_path)
            time_list.append(time)
    return cores_list, time_list

def speedup(cores_list,time_list):
    cores_list = np.array(cores_list)
    time_list = np.array(time_list)
    minimal_core_index = np.argmin(cores_list)
    speedup = time_list[minimal_core_index]/time_list

    return cores_list, speedup

def plot_speedup(cores_list, speeup, color='b'):
    fig, ax = plt.subplots()
    ax.plot(cores_list,speedup,'.',color=color)
    ax.set_xlabel("Numer of CPU")
    ax.set_ylabel("Speed-up")
    ax.grid()
    plt.show()


if __name__ == "__main__":
    color = sns.color_palette("cool", 12)[6]
    path = "/home/ccattin/Documents/HSP90_278K/test_parallel/"
    cores_list, time_list = get_time(path)
    cores_list, speedup = speedup(cores_list,time_list)
    plot_speedup(cores_list,speedup,color)
