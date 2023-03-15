#! /usr/bin/env python3
"""Extract and write distance"""
import matplotlib.pyplot as plt
import MDAnalysis as mda
import scipy as sp
import numpy as np
import pyemma
import os
import seaborn as sns

import distances__64CA_130CA__119CA_24CA as distances

import sys
sys.path.append('/home/ccattin/Documents/Code/water/')
import protein_volume

def load_volume_distance(data_volume,data_distance,state,number):
    
    config = protein_volume.configuration(number,state)
    
    # Volume
    file_volume = 'no_water_md_concatenated.txt'
    if state == 'GS':
        dir_volume = f"R6OA_GS{number:02d}_2021_11_19_Amber19SB_OPC_NaCl170mM_GMX_JeanZay"
    if state == 'ES':
        dir_volume = f"R46A_ES{number:02d}_2021_11_08_Amber19SB_OPC_NaCl170mM_GMX_JeanZay"
    path_volume = os.path.join(data_volume,dir_volume,file_volume)
    (time, 
     no_water, 
     volume_mean, 
     volume_error) = config.load_txt(output_name=path_volume)
    
    
    # Distance
    file_distance = f"{state}_{number}_distance.txt"
    path_distance = os.path.join(data_distance,file_distance)
    data = np.loadtxt(path_distance)
    d_64_CA_130_CA = data[:,0]
    d_119_CA_24_CA = data[:,1]
    
    #Seems that all the .xtc files have been saved less regularly 
    # only take 1/10 of the data on the volume
    return no_water[::10], d_64_CA_130_CA, d_119_CA_24_CA, time

def plot_volume_distance_2d(volume,d_64_CA_130_CA,d_119_CA_24_CA):
    fig, ax = plt.subplots()
    heatmap = ax.scatter(d_64_CA_130_CA,
                         d_119_CA_24_CA,
                         c=volume,
                         cmap='cool',
                         s=5)

    # Set axis labels and title
    ax.set_xlabel('64CA-130CA')
    ax.set_ylabel('119CA-24CA')

    # Add colorbar
    cbar = fig.colorbar(heatmap)
    cbar.set_label(r"Volume (nm$^3$)")

    # Show the plot
    plt.show()

def smoothing(volume,d_1,d_2,time,window_size):

    config = protein_volume.configuration(None,None)
    smoothed = config.smoothing(time=time,
                                 no_water=volume,
                                 window_size = window_size)
    return (smoothed,
            d_1[window_size//2:-window_size//2+1],
            d_2[window_size//2:-window_size//2+1])




if __name__ == '__main__':
    DATA = "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC"
    data_volume = '/home/ccattin/Documents/EAU/HSP90_simulation'
    data_distance = '/home/ccattin/Documents/Markov/volume_pressure/Output_distances_64CA-130CA_119CA-24CA'
    state = 'GS'
    number = 1
    (volume,
     d_64_CA_130_CA,
     d_119_CA_24_CA,
     time) = load_volume_distance(data_volume=data_volume,
                         data_distance=data_distance,
                         state=state,
                         number=number)
    plot_volume_distance_2d(volume=volume,
                            d_64_CA_130_CA=d_64_CA_130_CA,
                            d_119_CA_24_CA=d_119_CA_24_CA)
    
    window_size = [100,1000,3000]

    for window in window_size:
        smoothed, d1, d2 = smoothing(volume,d_64_CA_130_CA,d_119_CA_24_CA,time,window)
        plot_volume_distance_2d(volume=smoothed,
                                d_64_CA_130_CA=d1,
                                d_119_CA_24_CA=d2)