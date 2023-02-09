#! /usr/bin/env python3
"""
Extract output from gmx energy command
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class gromacs_output:
    
    def __init__(self,file_name):
        self.file_name = file_name

    def extract(self):
        x,y = np.loadtxt(self.file_name, comments=["@", "#"], unpack=True)
        return x, y

    def plot(self, xlabel, ylabel, output_name,color="b"):
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