#! /usr/bin/env python3
"""
Script to extract the molecular volume of water.
"""
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import extract_gmx_energy as xtract
import os


class eau(xtract.gromacs_output):
    """
    Class of the different water model (OPC and TIP3P studied)
    """

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #       Density to volume       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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
        # Avogadro constant
        Avogadro = 6.0221408e23
        # Molar mass of water
        molar_mass = 18.01528e-3

        x, density = self.extract()

        # Mecular volume in m^3/molecule
        volume = molar_mass / (density * Avogadro)
        # Convert to Angstrom^3
        volume = volume * 1e30

        # Compute the mean and the error
        density_mean_error = np.array(self.analyze())
        density_mean = density_mean_error[0]
        density_error = density_mean_error[1]

        volume_mean = molar_mass / (density_mean * Avogadro) * 1e30
        volume_error = molar_mass / (density_mean**2 * Avogadro) * 1e30

        return x, volume, volume_mean, volume_error

    def plot_volume(self, xlabel, ylabel, output_name, color="b", show=True):
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
        ax.plot(x, np.ones(len(x)) * volume_mean, linewidth=5, label="Mean")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        plt.savefig(output_name, format="pdf", dpi=300, bbox_inches="tight")
        if show:
            plt.show()

    def write(self, file_name):
        """Write the results in an .txt file

        Parameters
        ----------
        file_name : str
            Path where to save the file
        """
        # Compute results
        x, volume, volume_mean, volume_error = self.density_to_volume()
        # Save the results
        # The first line contain the mean volume and the estimated error
        # The line after that contains the time and the volume
        np.savetxt(
            file_name,
            np.array([x, volume]),
            header="{} {}".format(volume_mean, volume_error),
        )

    def read(self, file_name):
        """Read the data from .txt file

        Parameters
        ----------
        file_name : str
            Path where the file has been saved

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
        x, volume = np.loadtxt(file_name)
        with open(file_name) as f:
            line = f.readline()
        line = line.replace("#", "").replace("\n", "").split(" ")[1:]
        volume_mean = float(line[0])
        volume_error = float(line[1])

        return x, volume, volume_mean, volume_error

    #~~~~~~~~~~~~~~~~~~~~#
    #       Cutoff       #
    #~~~~~~~~~~~~~~~~~~~~#

    def cutoff_extract_volume(
        self, dir="/home/ccattin/Documents/EAU/Change_cutoff/TIP3P/300K/"
    ):
        """Extract the volume and the associated cut-off from simulations

        Parameters
        ----------
        dir : str, optional
            Directory where the pipeline has been run,
            by default "/home/ccattin/Documents/EAU/Change_cutoff/TIP3P/300K/"

        Returns
        -------
        cutoff_list : np.array
            Array containing the cutoff
        volume_mean_list : np.array
            Array containing the average volumes
        volume_error_list : np.array
            Array containing the error associated
        """
        # Initialize
        volume_mean_list = []
        volume_error_list = []
        cutoff_list = []
        # Get the cutoff directories name
        for name in os.listdir(dir):
            # Ignore other directories
            if "cutoff_" in name:
                cutoff = float(name.replace("cutoff_", ""))
                cutoff_list.append(cutoff)
                # Extract values
                (x, volume, volume_mean, volume_error) = self.read(
                    dir + name + "/analysis/time_volume_mean_error.txt"
                )
                volume_mean_list.append(volume_mean)
                volume_error_list.append(volume_error)

        return (
            np.array(cutoff_list),
            np.array(volume_mean_list),
            np.array(volume_error_list),
        )

    def cutoff_plot(self, cutoff_dir, color="b", show=True, output_path="cutoff.pdf"):
        """Plot the volume as a function of the used cut-off.

        Parameters
        ----------
        cutoff_dir : str
            Path where the pipeline has been run
        color : str, optional
            Color of the plot, by default "b"
        show : bool, optional
            Show the plot in an external windown or not, by default True
        output_path : str, optional
            Path to save the file, by default "cutoff.pdf"
        """
        # Extract data
        (cutoff_list, volume_mean_list, volume_error_list) = self.cutoff_extract_volume(
            cutoff_dir
        )
        fig, ax = plt.subplots()
        ax.errorbar(
            cutoff_list, volume_mean_list, volume_error_list, fmt=".", color=color
        )
        ax.grid()
        ax.set_xlabel("Cut-off (nm)")
        ax.set_ylabel(r"Molecular volume ($\AA^{3}$.molec$^{-1}$)")

        plt.savefig(output_path, dpi=300, bbox_inches="tight")

        # Plot in an external window
        if show:
            plt.show()

    def cutoff_write(self, dir, file_name):

        (cutoff_list,
            volume_mean_list,
            volume_error_list) = self.cutoff_extract_volume(dir)

        np.savetxt(file_name,
                    np.array([cutoff_list,
                            volume_mean_list,
                            volume_error_list]))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #       Simulation length       #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    def simulation_length_extract_volume(
        self, dir="/home/ccattin/Documents/EAU/Change_simulatio_length/TIP3P/300K/"
    ):
        """
        Extract the volume and the associated simulation length.

        Parameters
        ----------
        dir : str, optional
            Directory where the pipeline has been run,
            by default
            "/home/ccattin/Documents/EAU/Change_simulatio_length/TIP3P/300K/"

        Returns
        -------
        nsteps_list : np.array
            Array containing the number of steps done
        volume_mean_list : np.array
            Array containing the average volumes
        volume_error_list : np.array
            Array containing the error associated
        """
        # Initialize
        volume_mean_list = []
        volume_error_list = []
        nsteps_list = []
        # Get the nsteps directories name
        for name in os.listdir(dir):
            # Ignore other directories
            if "nsteps_" in name:
                nsteps = float(name.replace("nsteps_", ""))
                nsteps_list.append(nsteps)
                # Extract values
                (x, volume, volume_mean, volume_error) = self.read(
                    dir + name + "/analysis/time_volume_mean_error.txt"
                )
                volume_mean_list.append(volume_mean)
                volume_error_list.append(volume_error)

        return (
            np.array(nsteps_list),
            np.array(volume_mean_list),
            np.array(volume_error_list),
        )

    def simulation_length_plot(
        self,
        simulation_length_dir,
        color="b",
        show=True,
        output_path="simulation_length.pdf",
    ):
        """Plot the volume as a function of the used nsteps.

        Parameters
        ----------
        cutoff_list : str
            Path where the pipeline has been run
        color : str, optional
            Color of the plot, by default "b"
        show : bool, optional
            Show the plot in an external windown or not, by default True
        output_path : str, optional
            Path to save the file, by default "simulation_length.pdf"
        """
        # Extract data
        (
            simulation_length_list,
            volume_mean_list,
            volume_error_list,
        ) = self.simulation_length_extract_volume(simulation_length_dir)
        fig, ax = plt.subplots()
        ax.errorbar(
            simulation_length_list,
            volume_mean_list,
            volume_error_list,
            fmt=".",
            color=color,
        )
        ax.grid()
        ax.set_xlabel("nsteps")
        ax.set_ylabel(r"Molecular volume ($\AA^{3}$.molec$^{-1}$)")

        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        # Plot in an external window
        if show:
            plt.show()


if __name__ == "__main__":
    file_name = "/home/ccattin/Documents/EAU/Ludovic_parameters/computed_by_Come/OPC/300K/production/density.xvg"
    error_file = "/home/ccattin/Documents/EAU/Ludovic_parameters/computed_by_Come/OPC/300K/production/errest.xvg"

    xlabel = "Time (ps)"
    ylabel = r"Molecular volume ($\AA^{3}$.molec$^{-1}$)"
    output_name = "/home/ccattin/Documents/Code/outputs/molar_volume.pdf"
    output_result_name = (
        "/home/ccattin/Documents/Code/outputs/time_volume_mean_error.txt"
    )
    color = sns.color_palette("cool", 12)[6]
    show = True
    #######
    # OPC #
    #######
    OPC = eau(file_name, error_file)
    x, volume, volume_mean, volume_error = OPC.density_to_volume()

    print(
        """
        Mean molecular volume : {} +/- {} Angstrom/molec
        """.format(
            volume_mean, volume_error
        )
    )

    OPC.plot_volume(
        xlabel=xlabel, ylabel=ylabel, output_name=output_name, color=color, show=show
    )

    OPC.write(output_result_name)
    (x, volume, volume_mean, volume_error) = OPC.read(output_result_name)

    ##########
    # Cutoff #
    ##########
    cutoff_dir = "/home/ccattin/Documents/EAU/Change_cutoff/TIP3P/300K/"
    output_path = "/home/ccattin/Documents/Code/outputs/cutoff.pdf"
    file_name = "/home/ccattin/Documents/Code/outputs/cutoff.txt"
    (cutoff_list, volume_mean_list, volume_error_list) = OPC.cutoff_extract_volume(
        cutoff_dir
    )
    OPC.cutoff_plot(cutoff_dir, color, show=show, output_path=output_path)
    OPC.cutoff_write(cutoff_dir, file_name)

    #####################
    # Simulation length #
    #####################
    simulation_length_dir = (
        "/home/ccattin/Documents/EAU/Change_simulation_length/TIP3P/300K/"
    )
    output_path = "/home/ccattin/Documents/Code/outputs/simulation_length.pdf"
    (
        cutoff_list,
        volume_mean_list,
        volume_error_list,
    ) = OPC.simulation_length_extract_volume(simulation_length_dir)
    OPC.simulation_length_plot(
        simulation_length_dir, color, show=show, output_path=output_path
    )
