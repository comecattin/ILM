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
    #Get the path to the output pdf file
    parser.add_argument(
        "output_pdf",
        help="Path where to output the .pdf file of the molar volume plot",
        type=str,
        default="/home/ccattin/Documents/Code/outputs/molar_volume.pdf"
    )
    #Get the path to the output txt file
    parser.add_argument(
        "output_txt",
        help="Path where to output the .tx file of the results",
        type=str,
        default="/home/ccattin/Documents/Code/outputs/time_volume_mean_error.txt"
    )
    #Extract or not the cutoff
    parser.add_argument(
        "-c","--cutoff",
        help="Extract and plot the volume as a funtion of the cutoff",
        action="store_true"
    )
    parser.add_argument(
        "--cutoff-pdf",
        help="Path to save the pdf file of the plot of the volume as a function of the cutoff",
        default="cutoff.pdf",
        dest="cutoff_pdf"
    )
    parser.add_argument(
        "--cutoff-dir",
        help="Directory where the pipeline has been run.",
        default="./",
        dest="cutoff_dir"
    )
    args = parser.parse_args()

    #Define the output and plot properties
    color = sns.color_palette("cool", 12)[6]
    xlabel = "Time (ps)"
    ylabel = r"Molecular volume ($\AA^{3}$.molec$^{-1}$)"

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
                    output_name=args.output_pdf, color=color,
                    show=args.plot)
    #Save the result
    water_model.write(args.output_txt)

    #Cutoff part
    if args.cutoff:
        (cutoff_list,
            volume_mean_list,
            volume_error_list) = water_model.cutoff_extract_volume(
                                                args.cutoff_dir
                                                )
        print(args.cutoff_pdf)
        water_model.cutoff_plot(args.cutoff_dir,
                                color,
                                show=args.plot,
                                output_path=args.cutoff_pdf)


if __name__ == "__main__":
    main()