#! /usr/bin/env python3
"""
Command Line Interface for the water simulations.
"""
import sys
import os
import argparse
import seaborn as sns
import eau
import extract_gmx_energy

def main():
    """
    CLI main function
    """
    parser = argparse.ArgumentParser(
        "CLI for water simulations", description="Compute the molecular volume of water."
    )
    #File Path
    parser.add_argument(
        "file_name",
        help="Path to the .xvg gromacs output file to extract data from",
        type=str
        )
    #Error File path
    parser.add_argument(
        "error_file",
        help="Path to the .xvg gromacs analyze output file estimate of the error",
        type=str
        )
    #Output the graph in an external windown or not
    parser.add_argument(
        "-p","--plot",
        help="Plot the molecular volume as a function of time",
        action="store_true"
    )
    args = parser.parse_args()

    #Define the output and plot properties
    output_name = "/home/ccattin/Documents/Python/molar_volume.pdf"
    color = sns.color_palette("cool", 12)[6]
    xlabel = "Time (ps)"
    ylabel = r"Molecular volume ($\AA^{3}$.molec$^{-1}$)"
    output_name = "/home/ccattin/Documents/Python/molar_volume.pdf"

    #Open a water model
    water_model = eau.eau(args.file_name, args.error_file)
    #Compute its volume from its density
    x, volume, volume_mean, volume_error = water_model.density_to_volume()

    #Print out the result
    print(
        """
        Mean molecular volume : {} +/- {} Angstrom/molec
        """.format(volume_mean, volume_error)
        )
    
    #Plot the result
    water_model.plot_volume(xlabel=xlabel, ylabel=ylabel, 
                    output_name=output_name, color=color,
                    show=args.plot)


if __name__ == "__main__":
    main()