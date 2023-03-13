#! /usr/bin/env python3
"""
Extract the molecular volume of protein of multiple simualtion.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import extract_gmx_energy as xtract
import eau
import os
import subprocess
import seaborn as sns


class protein:
    """Class of a protein"""

    def __init__(self, path, water_file, water_volume):
        """Init function of the protein class

        Parameters
        ----------
        path : str
            path where the simulation output are
        water_file : str
            name of the water file
        water_volume : float
            volume of a water molecule in Angstrom**3
        """
        # Path
        self.path = path
        # Water file
        self.water_file = []
        self.water_file.append(os.path.join(path, water_file[0]))
        self.water_file.append(os.path.join(path, water_file[1]))
        # Water volume
        self.water_volume = water_volume

    def extract_volume_command(self, code):
        """Runs the gmx analysis script

        Parameters
        ----------
        code : str
            Path to the script
        """
        subprocess.call([code, self.path])

    def number_water(self):
        """Get the number of water molecule

        Returns
        -------
        water : list
            list containing number of water in two lists
                0 : exited state
                1 : ground state
        """
        # Initialization
        water = [[], []]
        # ES and GS
        for i in (0, 1):
            with open(self.water_file[i], "r") as f:
                for line in f:
                    water[i].append(int(line.split()[-1]))
        return water

    def remove_water(self):
        """Remove the volume of water"""
        # Get the number of water molecules
        water_number = self.number_water()

        # For every configuration
        for name in os.listdir(self.path):
            # Don't consider .txt files
            if ".txt" in name:
                continue
            # Get the configuration number
            conf = int(name[7:9])
            # Create an empty list to hold the trajectory data
            traj_data = []
            # For every trajectory
            #   (10 trajectories each time)
            for traj in range(1, 11):

                # File names
                volume_file = os.path.join(
                    self.path, name, "volume_md_" + str(traj) + ".xvg"
                )
                error_file = os.path.join(
                    self.path, name, "errest_md_" + str(traj) + ".xvg"
                )
                output = os.path.join(
                    self.path, name, "no_water_md_" + str(traj) + ".txt"
                )

                # Extract the volume
                volume = xtract.gromacs_output(
                    file_name=volume_file, error_file=error_file
                )
                time, volume_array = volume.extract()

                # Remove the water
                if "ES" in name:
                    state='ES'
                    no_water = (
                        volume_array
                        - self.water_volume * 1e-3 * water_number[0][conf - 1]
                    )
                if "GS" in name:
                    state='GS'
                    no_water = (
                        volume_array
                        - self.water_volume * 1e-3 * water_number[1][conf - 1]
                    )
                # Append the trajectory data to the list
                traj_data.append(np.array([time, no_water]))
                # Save the result
                np.savetxt(output, np.array([time, no_water]))

            # Concatenate the trajectory data
            output_concatenated = os.path.join(
                self.path, name, "no_water_md_concatenated.txt"
            )
            configuration_i = configuration(number=conf,state=state)
            configuration_i.save_txt(traj_data=traj_data,output_name=output_concatenated)
            

class configuration():

    def __init__(self,number,state):
        self.num=number
        self.state=state
    
    def smoothing(self,time,no_water,window_size):
        smooth = []
        for window in window_size:
            weights = np.repeat(1.0, window) / window
            smoothing = np.convolve(no_water, weights, 'valid')
            smooth.append(smoothing)
        return smooth

    def save_txt(self,traj_data,output_name):
        concatenated = np.concatenate(traj_data, axis=1)
        # Compute the mean and the error
        volume_mean = np.mean(concatenated[1])
        volume_mean_error = stats.sem(concatenated[1], axis=None)
        # Save to a single file
        np.savetxt(
            output_name,
            concatenated,
            header="{} {}".format(volume_mean, volume_mean_error),
        )
    
    def load_txt(self,output_name):
        time,no_water = np.loadtxt(output_name)
        with open(output_name) as f:
            line = f.readline()
        line = line.replace("#", "").replace("\n", "").split(" ")[1:]
        volume_mean = float(line[0])
        volume_error = float(line[1])
        
        return time,no_water,volume_mean,volume_error
    
    def plot(self,time,no_water,volume_mean,volume_error,color,
             smoothing,
             window_size):
        fig,ax = plt.subplots()
        ax.plot(time,no_water,color=color[-1])
        volume_mean_array = volume_mean*np.ones(len(time))
        ax.plot(time,volume_mean_array,label="Mean",color=color[-2])
        ax.errorbar(time,volume_mean_array,yerr=volume_error,color=color[1])
        for i, smooth in enumerate(smoothing):
            ax.plot(smooth,
                    color=color[i],
                    label=f'Window size : {window_size[i]}')
        ax.legend()
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(r"Volume (nm$^{3}$)")
        plt.show()

if __name__ == "__main__":
    # Definition
    #   Pathes
    path = "/home/ccattin/Documents/EAU/HSP90_simulation"
    water_file = ["ES.txt", "GS.txt"]
    file_name = "/home/ccattin/Documents/EAU/Elisa_parameters/computed_by_Come/TIP3P/300K/analysis/density.xvg"
    error_file = "/home/ccattin/Documents/EAU/Elisa_parameters/computed_by_Come/TIP3P/300K/analysis/errest.xvg"
    code = "/home/ccattin/Documents/Code/GMX/analysis_water_protein"

    # Get the volume of water
    TIP3P_300K = eau.eau(file_name, error_file)
    water_volume = TIP3P_300K.density_to_volume()[2]
    # water_volume = 30.72669350142569

    # Initialize a protein object
    HSP90 = protein(path=path, water_file=water_file, water_volume=water_volume)

    # Run script to extract all the volume
    # HSP90.extract_volume_command(code=code)

    # Get the number of water molecule
    ES, GS = HSP90.number_water()

    # Remove the water volume
    HSP90.remove_water()

    # Plot the result
    output_name = path+'/R6OA_GS01_2021_11_19_Amber19SB_OPC_NaCl170mM_GMX_JeanZay/no_water_md_concatenated.txt'
    window_size=[3,6,12]
    #   Color
    color = sns.color_palette("cool", len(window_size)+2)
    GS_1 = configuration(number=1,state='GS')
    (time,
     no_water,
     volume_mean,
     volume_error
     ) = GS_1.load_txt(output_name=output_name)
    smooth = GS_1.smoothing(time,no_water,window_size)
    GS_1.plot(time=time,
                        no_water=no_water,
                        volume_mean=volume_mean,
                        volume_error=volume_error,
                        color=color,
                        smoothing=smooth,
                        window_size=window_size)

