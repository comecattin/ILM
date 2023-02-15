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
    #Analysis or not the molecular volume
    parser.add_argument(
        "--volume",
        dest="output_volume",
        help="Analysis of the molecular volume",
        action="store_true"
    )
    #File Path
    parser.add_argument(
        "-f",
        "--file-name",
        help="Path to the .xvg gromacs output file to extract data from",
        type=str,
        dest="file_name"
        )
    #Error File path
    parser.add_argument(
        "-e",
        "--error-file",
        help="Path to the .xvg gromacs analyze output file estimate of the error",
        type=str,
        dest="error_file",
        )
    #Output the graph in an external windown or not
    parser.add_argument(
        "-p","--plot",
        help="Plot the molecular volume as a function of time",
        action="store_true"
    )
    #Get the path to the output pdf file
    parser.add_argument(
        "--pdf",
        dest="output_pdf",
        help="Path where to output the .pdf file of the molar volume plot",
        type=str,
        default="molar_volume.pdf"
    )
    #Get the path to the output txt file
    parser.add_argument(
        "--txt",
        dest="output_txt",
        help="Path where to output the .txt file of the results",
        type=str,
        default="time_volume_mean_error.txt"
    )
    #Extract or not the cutoff
    parser.add_argument(
        "-c","--cutoff",
        help="Extract and plot the volume as a funtion of the cutoff",
        action="store_true"
    )
    #Save the .pdf from the cutoff analysis
    parser.add_argument(
        "--cutoff-pdf",
        help="Path to save the pdf file of the plot of the volume as a function of the cutoff",
        default="cutoff.pdf",
        dest="cutoff_pdf"
    )
    #Cutoff directory
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

    #Raise error if no options are provided
    if args.output_volume == False and args.cutoff == False:
        raise Exception("No options provided. Please use the -h option.")

    else:

        #Molecular Volume part
        if args.output_volume:
            
            #Raise error if some file are missing
            if args.file_name == None or args.error_file == None:
                raise Exception("No file or no error file provided. Please use the -h option.")
            
            else:
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
            water_model = eau.eau()
            (cutoff_list,
                volume_mean_list,
                volume_error_list) = water_model.cutoff_extract_volume(
                                                    args.cutoff_dir
                                                    )
            water_model.cutoff_plot(args.cutoff_dir,
                                    color,
                                    show=args.plot,
                                    output_path=args.cutoff_pdf)


if __name__ == "__main__":
    main()