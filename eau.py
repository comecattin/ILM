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

    def density_to_volume(self):
        """Convert the density calculated with gromacs to molecular volume

        Returns
        -------
        x : np.array
            Simulation time
        volume : np.array
            Molecular volume as a function of the simulation time
            Unit : Angstrom^3 / molecule
        volume_mean : float
            Time average of the volume
        volume_error : float
            Volume associated error
        """
        #Avogadro constant
        Avogadro = 6.0221408e23
        #Molar mass of water
        molar_mass = 18.01528e-3

        x, density = self.extract()
        
        #Mecular volume in m^3/molecule
        volume = molar_mass/(density*Avogadro)
        #Convert to Angstrom^3
        volume = volume*1e30
        
        #Compute the mean and the error
        density_mean_error = np.array(self.analyze())
        density_mean = density_mean_error[0]
        density_error  = density_mean_error[1]

        volume_mean = molar_mass/(density_mean*Avogadro)*1e30
        volume_error = molar_mass/(density_mean**2 *Avogadro)*1e30

        return x, volume, volume_mean, volume_error
    
    def plot_volume(self, xlabel, ylabel, output_name,color="b",show=True):
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
            Color of the plot, by default "b"
        show : bool, optional
            Show the plot in an external windown or not, by default True
        """
        x, volume, volume_mean, volume_error = self.density_to_volume()
        fig, ax = plt.subplots()
        ax.plot(x, volume, color=color, label="Molecular volume")
        ax.plot(x, np.ones(len(x))*volume_mean,linewidth=5, label="Mean")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        plt.savefig(output_name, format="pdf", dpi=300, bbox_inches='tight')
        if show:
            plt.show()


if __name__ == "__main__":
    file_name = "/home/ccattin/Documents/EAU/OPC/production/density.xvg"
    error_file = "/home/ccattin/Documents/EAU/OPC/production/errest.xvg"

    xlabel = "Time (ps)"
    ylabel = r"Molecular volume ($\AA^{3}$.molec$^{-1}$)"
    output_name = "/home/ccattin/Documents/Python/molar_volume.pdf"
    color = sns.color_palette("cool", 12)[6]

    #######
    # OPC #
    #######
    OPC = eau(file_name, error_file)
    x, volume, volume_mean, volume_error = OPC.density_to_volume()

    print(
        """
        Mean molecular volume : {} +/- {} Angstrom/molec
        """.format(volume_mean, volume_error))
    
    OPC.plot_volume(xlabel=xlabel, ylabel=ylabel, 
                    output_name=output_name, color=color)