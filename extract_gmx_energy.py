#! /usr/bin/env python3
"""
Extract output from gmx energy command
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class gromacs_output:
    """
    Class to extract data from gromacs ernergy outputs
    """
    def __init__(self,file_name):
        """
        Initialize the class

        Parameters
        ----------
        file_name : str
            Path to the .xvg gromacs output file to extract data from
        """
        self.file_name = file_name

    def extract(self):
        """
        Extract gromacs data

        Returns
        -------
        x : np.array
            Array containing the extracted time simulation
        y : np.array
            Array containing the extracted observable 
        """
        x,y = np.loadtxt(self.file_name, comments=["@", "#"], unpack=True)
        return x, y

    def plot(self, xlabel, ylabel, output_name,color="b"):
        """
        Plot the extracted gromacs output data

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
        x,y = self.extract()
        fig, ax = plt.subplots()
        ax.plot(x,y,color=color)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        plt.savefig(output_name, format="pdf", dpi=300, bbox_inches='tight')
        plt.show()


if __name__ == "__main__":

    file_name = "/home/ccattin/Documents/EAU/OPC/production/density.xvg"

    xlabel = "Time (ps)"
    ylabel = r"Density (kg.m$^{-3}$)"
    output_name = "density.pdf"
    color = sns.color_palette("cool", 12)[6]

    density = gromacs_output(file_name)
    x, y = density.extract()
    density.plot(xlabel=xlabel, ylabel=ylabel, output_name=output_name, color=color)