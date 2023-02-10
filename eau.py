#! /usr/bin/env python3
import seaborn as sns
import matplotlib.pyplot as plt
import extract_gmx_energy as xtract

class eau(xtract.gromacs_output):

    def __init__(self, file_name, molar_mass):
        self.file_name = file_name
        self.molar_mass = molar_mass


    def density_to_volume(self):
        Avogadro = 6.0221408e23
        x, density = self.extract()
        volume = self.molar_mass/(density*Avogadro)
        volume = volume*1e30
        return x, volume
    
    def plot_volume(self, xlabel, ylabel, output_name,color="b"):
        x, volume = self.density_to_volume()
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
    x, volume = OPC.density_to_volume()
    OPC.plot_volume(xlabel=xlabel, ylabel=ylabel, output_name=output_name, color=color)