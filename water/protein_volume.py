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

    def __init__(self, path = '',
                 water_file = ['',''],
                 water_volume = None,
                 water_volume_error = None):
        """Init function of the protein class

        Parameters
        ----------
        path : str
            Path where the simulation output are
        water_file : list of str
            Name of the water files
        water_volume : float
            Volume of a water molecule in Angstrom**3
        water_volume : float
            Error on the estimation of the volume
            of a water molecule in Angstrom**3
        """
        # Path
        self.path = path
        # Water file
        self.water_file = []
        self.water_file.append(os.path.join(path, water_file[0]))
        self.water_file.append(os.path.join(path, water_file[1]))
        # Water volume and its error
        self.water_volume = water_volume
        self.water_volume_error = water_volume_error

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
                    state = "ES"
                    no_water = (
                        volume_array
                        - self.water_volume * 1e-3 * water_number[0][conf - 1]
                    )
                    # Get the error
                    error = self.water_volume_error * 1e-3 * water_number[0][conf - 1]
                if "GS" in name:
                    state = "GS"
                    no_water = (
                        volume_array
                        - self.water_volume * 1e-3 * water_number[1][conf - 1]
                    )
                    # Get the error
                    error = self.water_volume_error * 1e-3 * water_number[1][conf - 1]
                # Append the trajectory data to the list
                traj_data.append(np.array([time, no_water]))
                # Save the result
                np.savetxt(output, np.array([time, no_water]))

            # Concatenate the trajectory data
            output_concatenated = os.path.join(
                self.path, name, "no_water_md_concatenated.txt"
            )
            configuration_i = configuration(number=conf, state=state)
            configuration_i.save_txt(
                traj_data=traj_data, error=error, output_name=output_concatenated
            )

    def mean_GS_ES(self, state):
        """Compute the mean over all the ground or exctited states

        Parameters
        ----------
        state : str
            State to compute the mean : 'GS' or 'ES'
        Returns
        -------
        mean_list : list
            Mean values of the volume of the protein of each configuration
        error_list : list
            Error associated
        """
        # Initialize
        mean_list = []
        error_list = []
        # For every configuration
        for name in os.listdir(self.path):
            # Ignore the .txt files
            if ".txt" in name:
                continue
            # Get only the GS configurations
            if state in name:
                # Get the configuration number
                conf = int(name[7:9])

                # Extract the mean volume and the associated error
                output_concatenated = os.path.join(
                    self.path, name, "no_water_md_concatenated.txt"
                )
                state_i = configuration(conf, state)
                (time, no_water, volume_mean, volume_error) = state_i.load_txt(
                    output_concatenated
                )

                # Append the list
                mean_list.append(volume_mean)
                error_list.append(volume_error)

        return mean_list, error_list

    def plot_state_mean(self, mean_list, error_list, output_path, color):
        """Plot of the ground or excited state mean with their associated error.

        Parameters
        ----------
        mean_list : list
            Mean values of the volume of the protein of each configuration
        error_list : list
            Error associated
        output_path : str
            Path where to save the output .pdf
        color : str
            Color of the plot
        """
        fig, ax = plt.subplots()
        # Plot the mean of the different configurations
        ax.errorbar(range(20), mean_list, yerr=error_list, fmt=".", color=color)
        ax.set_xlabel("Configuration")
        ax.set_ylabel(r"Mean volume (nm$^3$)")
        ax.set_xticks(range(20))
        ax.grid()
        # Add a text to have the mean over all the ground state
        ax.text(
            0.95,
            0.95,
            f"Mean: {np.mean(mean_list):.4f} nm$^3$",
            transform=plt.gca().transAxes,
            ha="right",
            va="top",
        )
        # Save and show
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.show()


class configuration:
    """Configuration class"""

    def __init__(self, number, state):
        """Init function for the configuration class

        Parameters
        ----------
        number : float
            Number of the configuration (between 1 and 20 usually)
        state : str
            Ground state ('GS') or excited state ('ES')
        """
        self.num = number
        self.state = state

    def smoothing(self, time, no_water, window_size):
        """Do a rolling average

        Parameters
        ----------
        time : np.array
            time of the simulation
        no_water : np.array
            volume of the protein without water
        window_size : list or int
            The different window size

        Returns
        -------
        smooth : np.array
            Smoothed values
        """

        # If only one value of window size
        if type(window_size) == int:
            weights = np.ones(window_size) / window_size
            smoothing = np.convolve(no_water, weights, "valid")

            return smoothing

        else:
            # Initialize
            smooth = []
            # For each window values
            for window in window_size:
                weights = np.repeat(1.0, window) / window
                smoothing = np.convolve(no_water, weights, "valid")
                smooth.append(smoothing)
    
                return smooth

    def save_txt(self, traj_data, error, output_name):
        """Save the value of the volume of the protein without water

        Parameters
        ----------
        traj_data : np.array
            Data containing the volume in each trajectories
        error : float
            Error related with the estimation of the water volume
        output_name : str
            Path where to save the file
        """
        concatenated = np.concatenate(traj_data, axis=1)
        # Compute the mean and the error
        # Error on the mean
        volume_mean = np.mean(concatenated[1])
        volume_mean_error = stats.sem(concatenated[1], axis=None)
        # Systematic error on the water volume
        volume_mean_error += error
        # Save to a single file
        np.savetxt(
            output_name,
            concatenated,
            header="{} {}".format(volume_mean, volume_mean_error),
        )

    def load_txt(self, output_name):
        """Load the result previously saved with the function save_txt

        Parameters
        ----------
        output_name : str
            Path where the file has been saved

        Returns
        -------
        time : np.array
            time of the simulation
        no_water : np.array
            volume of the protein without water
        volume_mean : float
            Mean of the volume of the protein without water
        volume_error : float
            Associated error to the mean
        """
        # Load the file
        time, no_water = np.loadtxt(output_name)
        # Get the mean and the error
        with open(output_name) as f:
            line = f.readline()
        line = line.replace("#", "").replace("\n", "").split(" ")[1:]
        volume_mean = float(line[0])
        volume_error = float(line[1])

        return time, no_water, volume_mean, volume_error

    def plot(
        self,
        time,
        no_water,
        volume_mean,
        volume_error,
        color,
        smoothing,
        window_size,
        output_path,
    ):
        """Plot of the volume as a function of the time

        Parameters
        ----------
        time : np.array
            time of the simulation
        no_water : np.array
            volume of the protein without water
        volume_mean : float
            Mean of the volume of the protein without water
        volume_error : float
            Associated error to the mean
        color : list
            Color of the different plot
        smoothing : list
            Smoothed value of the volume at different window size
        window_size : list
            The different windw sizes for the smoothing
        output_path : str
            Path where to save the output .pdf
        """

        fig, ax = plt.subplots()

        # Plot of the raw volume
        ax.plot(time, no_water, color=color[-1])

        # Plot of the mean and its associated error
        volume_mean_array = volume_mean * np.ones(len(time))
        ax.plot(time, volume_mean_array, label="Mean", color=color[-2])
        ax.errorbar(time, volume_mean_array, yerr=volume_error, color=color[1])

        # Plot of the smoothed values for each window value
        for i, smooth in enumerate(smoothing):
            ax.plot(smooth, color=color[i], label=f"Window size : {window_size[i]}")
        ax.legend()
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(r"Volume (nm$^{3}$)")
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.show()


if __name__ == "__main__":
    # Definition
    #   Pathes
    path = "/home/ccattin/Documents/EAU/HSP90_simulation"
    water_file = ["ES.txt", "GS.txt"]
    file_name = "/home/ccattin/Documents/EAU/Elisa_parameters/computed_by_Come/TIP3P/300K/analysis/density.xvg"
    error_file = "/home/ccattin/Documents/EAU/Elisa_parameters/computed_by_Come/TIP3P/300K/analysis/errest.xvg"
    code = "/home/ccattin/Documents/Code/GMX/analysis_water_protein"

    # Get the volume and its error of water
    TIP3P_300K = eau.eau(file_name, error_file)
    water_volume, water_volume_error = TIP3P_300K.density_to_volume()[2:4]
    # water_volume = 30.72669350142569

    # Initialize a protein object
    HSP90 = protein(
        path=path,
        water_file=water_file,
        water_volume=water_volume,
        water_volume_error=water_volume_error,
    )

    # Run script to extract all the volume
    # HSP90.extract_volume_command(code=code)

    # Get the number of water molecule
    ES, GS = HSP90.number_water()

    # Remove the water volume
    HSP90.remove_water()

    # Plot the result
    output_name = (
        path
        + "/R6OA_GS01_2021_11_19_Amber19SB_OPC_NaCl170mM_GMX_JeanZay/no_water_md_concatenated.txt"
    )
    output_plot = "/home/ccattin/Documents/Code/outputs/no_water.pdf"
    window_size = [100, 1000, 10000]
    #   Color
    color = sns.color_palette("cool", len(window_size) + 2)
    # Initialize a configuration
    GS_1 = configuration(number=1, state="GS")
    # Get its results
    (time, no_water, volume_mean, volume_error) = GS_1.load_txt(output_name=output_name)
    # Smooth the results
    smooth = GS_1.smoothing(time, no_water, window_size)
    # Plot the results
    GS_1.plot(
        time=time,
        no_water=no_water,
        volume_mean=volume_mean,
        volume_error=volume_error,
        color=color,
        smoothing=smooth,
        window_size=window_size,
        output_path=output_plot,
    )

    # Mean on the ground states
    output_plot = "/home/ccattin/Documents/Code/outputs/mean_GS.pdf"
    # Get the means and the associated errors
    mean_list, error_list = HSP90.mean_GS_ES("GS")
    color = sns.color_palette("cool", 12)
    # Plot the result
    HSP90.plot_state_mean(
        mean_list=mean_list,
        error_list=error_list,
        color=color[6],
        output_path=output_plot,
    )

    # Mean on the excited states
    output_plot = "/home/ccattin/Documents/Code/outputs/mean_ES.pdf"
    # Get the means and the associated errors
    mean_list, error_list = HSP90.mean_GS_ES("ES")
    color = sns.color_palette("cool", 12)
    # Plot the result
    HSP90.plot_state_mean(
        mean_list=mean_list,
        error_list=error_list,
        color=color[6],
        output_path=output_plot,
    )
