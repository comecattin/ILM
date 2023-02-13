#! /usr/bin/env python3
"""
Extract output from gmx energy command
"""
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns

class gromacs_output:
    """
    Class to extract data from gromacs ernergy outputs
    """
    def __init__(self, file_name, error_file="errest.xvg"):
        """
        Initialize the class

        Parameters
        ----------
        file_name : str
            Path to the .xvg gromacs output file to extract data from
        error_file : str
            Path to the .xvg gromacs analyze output file estimate of the error,
            by default "errest.xvg"
        """
        self.file_name = file_name
        self.error_file = error_file

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

    def analyze(self):
        x,y = self.extract()
        #Estimate the mean
        mean_val = np.mean(y)
        #Estimate the error
        estimate_error_line = subprocess.getoutput(
            """grep "e " {}""".format(self.error_file)
            )
        error = float(estimate_error_line[-10:-1])
        return mean_val, error



if __name__ == "__main__":

    file_name = "/home/ccattin/Documents/EAU/OPC/production/density.xvg"
    error_file = "/home/ccattin/Documents/EAU/OPC/production/errest.xvg"

    xlabel = "Time (ps)"
    ylabel = r"Density (kg.m$^{-3}$)"
    output_name = "/home/ccattin/Documents/Python/outputs/density.pdf"
    color = sns.color_palette("cool", 12)[6]

    density = gromacs_output(file_name,error_file)
    x, y = density.extract()
    density.plot(xlabel=xlabel, ylabel=ylabel, output_name=output_name, color=color)
    
    mean_val, error = density.analyze()
    print(
        """
        Mean value: {}
        Estimated error: {}
        """.format(mean_val,error)
        )
