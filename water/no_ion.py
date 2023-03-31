#! /usr/bin/env python3
"""
Extract the volume of the water + ions and get the protein volume
"""
import numpy as np
import matplotlib.pyplot as plt
import extract_gmx_energy as xtract
import os
import protein_volume
import seaborn as sns


class no_ion:
    """Class for the extraction of volume without ions"""

    def __init__(self, path_water, path_protein):
        """Init function of the class

        Parameters
        ----------
        path_water : str
            Path to the data of the water simulation without ions
        path_protein : str
            Path to the data of the protein simulation
        """
        self.path_water = path_water
        self.path_protein = path_protein
        # Initialyze a dictionnary that will contain all the properties
        self.state_dict = {}

    def mean_volume(self):
        """Compute the mean volume of the water + ions

        Update the state_dict with the mean and the error on the volume
        Update the water dictionnary with all the data on the volume of the system
        """
        # Init
        state_dict = {}
        water = {}
        # Every directory in the data path
        for name in os.listdir(self.path_water):
            state = name[5:9]
            # File path
            volume_file = os.path.join(self.path_water, name, "volume_md_11.xvg")
            error_file = os.path.join(self.path_water, name, "errest_md_11.xvg")
            ouput = xtract.gromacs_output(file_name=volume_file, error_file=error_file)
            # Get the mean and error and append
            mean_error = ouput.analyze()
            state_dict[state] = list(mean_error)
            # Get all the data
            water[state] = ouput.extract()
        self.state_dict = state_dict
        self.water = water

    def get_total_volume_simulation(self, concatenate=False):
        """Get the total volume of the protein simulation (protein + water + ions)

        Parameters
        ----------
        concatenate : bool, optional
            Concatenate or not the simulation results in a .txt file,
                by default False

        Update the state_dict with a tupple containig :
        time : np.array
            Contain the time of the simulation
        Vprot : np.array
            Volume of the system
        mean : float
            Time average of the volume
        error : float
            Error on the time average of the volume
        """
        # Every directories in the data path
        for name in os.listdir(path_protein):
            if name[5:9] in self.state_dict:
                config = protein_volume.configuration(int(name[7:9]), name[5:7])
                # Concatenate the result
                if concatenate:
                    config.concatenate(self.path_protein)
                concatenated = os.path.join(
                    self.path_protein, name, "md_concatenated.txt"
                )
                # Extract the result
                t_Vprot_mean_error = config.load_txt(concatenated)
                self.state_dict[name[5:9]].append(t_Vprot_mean_error)

    def no_water(self):
        """Get the volume of the protein without water and ions

        Create and update the waterless dict with a list containnig:
        time : np.array
            Contain the time of the simulation
        Vprot : np.array
            Volume of the protein
        """
        self.waterless = {}
        # Every configuration in the state_dict
        for state in self.state_dict:
            # Extract data
            # in nm^3
            water_volume = self.state_dict[state][0]
            water_error = self.state_dict[state][1]

            time = self.state_dict[state][-1][0]
            protein_volume = self.state_dict[state][-1][1]

            # Remove the water
            waterless = protein_volume - water_volume
            # Append to the dictionnary
            self.waterless[state] = [time, waterless]

    def plot_volume_state(self, output):
        """Plot the results from the no_water function

        Volume of the protein as a function of time

        Parameters
        ----------
        output : str
            Path to the output .pdf file
        """
        # Init
        fig, axs = plt.subplots(len(self.waterless), sharex=True, figsize=(5, 7))
        color = sns.color_palette("cool", 12)

        # Do a plot per state
        for i, state in enumerate(self.waterless):
            # Time
            x = self.waterless[state][0]
            # Volume
            y = self.waterless[state][1]

            # Smooth the result
            config = protein_volume.configuration(None, None)
            smooth = config.smoothing(x, y, 100)

            # Volume vs time
            axs[i].plot(x, y, color=color[6])
            # Smoothed vs time
            axs[i].plot(smooth, label="Smoothed")

            # Put the state name in the plot
            axs[i].text(0, 32.5, state)
            axs[i].set_ylim(24, 34)
            axs[i].set_xlim(0, 100000)
        axs[0].legend()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.text(0.5, 0.04, "Time (ps)", ha="center", va="center")
        fig.text(
            0.03, 0.5, r"Volume (nm$^3$)", ha="center", va="center", rotation="vertical"
        )
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.show()

    def results(self):
        """Print the results"""
        # Header
        print(
            "State    Water    Error Water    Total    Total error    Protein    Protein error"
        )
        for state in self.state_dict:
            # Everything in nm^3
            # Water volume
            water = self.state_dict[state][0]
            # Error on the water volume
            water_error = self.state_dict[state][1]
            # Total volume of the system (protein + water + ions)
            total = self.state_dict[state][-1][2]
            # Error on the total system volume
            total_error = self.state_dict[state][-1][3]
            # Protein only volume
            protein = total - water
            # Associated error
            protein_error = total_error + water_error
            print(
                f"{state}: {water} {water_error} {total} {total_error} {protein} {protein_error}"
            )
        # Difference between the ES and GS
        diff, error_diff = self.error_state
        print(f"ES - GS = {diff} +/- {error_diff}")

    def confidence_intervals(self, n1, s1, n2, s2):
        """Compute the confidence interval of the difference

        Parameters
        ----------
        n1 : int
            Number of data points for the first object
        s1 : float
            Variance of the first object
        n2 : int
            Number of data points for the second object
        s2 : float
            Variance of the second object

        Returns
        -------
        error : float
            Associated error to the difference
        """

        sp2 = ((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2)
        error = np.sqrt(sp2 * (1 / n1 + 1 / n2)) * 3.182

        return error

    def error_using_variance(self):
        """Compute the error using the variance

        Create and update the error and error_state dictionnary
        """
        # Init
        self.error = {}
        self.error_state = {}
        excited = []
        ground = []

        # Loop over all the state
        for state in self.state_dict:

            # Get data on the solvent
            water = self.water[state][1]
            n_water = len(water)
            var_water = np.var(water)

            # Get data on the total system
            total = self.state_dict[state][-1][1]
            n_total = len(total)
            var_total = np.var(total)

            # Compute the error and append
            error = self.confidence_intervals(n_total, var_total, n_water, var_water)
            self.error[state] = error

            # Get the protein-only volume
            water_mean = self.state_dict[state][0]
            total_mean = total = self.state_dict[state][-1][2]
            protein = total_mean - water_mean
            # Append to the correct list
            if "ES" in state:
                excited.append(protein)
            if "GS" in state:
                ground.append(protein)

        # Gather data for all the GS and ES
        mean_ES = np.mean(excited)
        mean_GS = np.mean(ground)
        n_ES = len(excited)
        n_GS = len(ground)
        var_ES = np.var(excited)
        var_GS = np.var(ground)
        # Compute the difference of the mean and the associated error
        diff = mean_ES - mean_GS
        error_diff = self.confidence_intervals(n_ES, var_ES, n_GS, var_GS)
        # Append
        self.error_state = [diff, error_diff]


if __name__ == "__main__":
    path_water = "/home/ccattin/Documents/EAU/NO_COUNTERIONS"
    path_protein = "/home/ccattin/Documents/EAU/HSP90_simulation"
    output = "/home/ccattin/Documents/Code/outputs/no_ion.pdf"
    NoIon = no_ion(path_water, path_protein)
    NoIon.mean_volume()
    NoIon.get_total_volume_simulation()
    NoIon.no_water()
    # NoIon.plot_volume_state(output=output)
    NoIon.error_using_variance()
    NoIon.results()
