#! /usr/bin/env python3
"""
Script to extract the molecular volume of water.
"""
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import extract_gmx_energy as xtract

class eau(xtract.gromacs_output):
    """
    Class of the different water model (OPC and TIP3P studied)
    """
    def __init__(self, file_name, molar_mass):
        """
        Initialize function

        Parameters
        ----------
        file_name : str
            path to the .xvg gromacs output file to extract data from 
        molar_mass : float
            molar mass of the water
        """
        self.file_name = file_name
        self.molar_mass = molar_mass


    def density_to_volume(self):
        """Convert the density calculated with gromacs to molecular volume

        Returns
        -------
        x : np.array
            Simulation time
        volume : np.array
            Molecular volume as a function of the simulation time
            Unit : Angstrom^3 / molecule
        """
        #Avogadro constant
        Avogadro = 6.0221408e23
        x, density = self.extract()
        #Mecular volume in m^3/molecule
        volume = self.molar_mass/(density*Avogadro)
        #Convert to Angstrom^3
        volume = volume*1e30
        #Compute the mean
        mean_volume = np.mean(volume)
        return x, volume, mean_volume
    
    def plot_volume(self, xlabel, ylabel, output_name,color="b"):
        """
        Plot the molecular volume as a function of time

        Parameters
        ----------
        xlabel : str
            Label of the x-axis
        ylabel : str
            Label of the y-axis
        output_name : str
            Path to save the file
        color : str, optional
            color of the plot, by default "b"
        """
        x, volume, mean_volume = self.density_to_volume()
        fig, ax = plt.subplots()
        ax.plot(x, volume, color=color)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.savefig(output_name, format="pdf", dpi=300, bbox_inches='tight')
        plt.show()


if __name__ == "__main__":
    file_name = "/home/ccattin/Documents/EAU/OPC/production/density.xvg"

    xlabel = "Time (ps)"
    ylabel = r"Molecular volume ($\AA^{3}$.molec$^{-1}$)"
    output_name = "/home/ccattin/Documents/Python/molar_volume.pdf"
    color = sns.color_palette("cool", 12)[6]

    molar_mass = 18.01528e-3

    OPC = eau(file_name, molar_mass)
    x, volume, mean_volume = OPC.density_to_volume()

    print(
        "\nMean molecular volume : {} Angstrom/molec\n".format(mean_volume))
    
    OPC.plot_volume(xlabel=xlabel, ylabel=ylabel, 
                    output_name=output_name, color=color)