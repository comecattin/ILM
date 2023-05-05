#! /usr/bin/env python3
"""LETHE's parser"""

import argparse

def parsing():
    parser = argparse.ArgumentParser(
        "CLI LETHE program",
        description="Automatize MSM building using PyEmma.",
    )
    # List of files
    parser.add_argument(
        '-f',
        '--files',
        required=True,
        nargs='+',
        help='List of files'
        )
    # Distances to add to the feat
    parser.add_argument(
        '-d',
        '--distances',
        nargs='+',
        help='List of the pair names'
    )
    # Topology .pdb file
    parser.add_argument(
        '-t',
        '--topology',
        required=True,
        help='Reference .pdb file'
    )
    # What properties to plot
    parser.add_argument(
        '-p',
        '--plot',
        nargs='+',
        help='Plot wanted (feat_hist, density_energy)'
    )
    # Do not display the plot 
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Do not display plots'
    )
    # Output directory to save plots
    parser.add_argument(
        '-o',
        '--outdir',
        help='Path to save the plots'
    )
    # Temperature in K of the system
    parser.add_argument(
        '--T',
        help='Temperature of the system'
    )
    # PCA dimension reduction
    parser.add_argument(
        '--pca',
        help='Do a PCA dimension reduction',
        action='store_true'
    )
    # TICA dimension reduction
    parser.add_argument(
        '--tica',
        help='Do a TICA dimension reduction',
        action='store_true'
    )
    # Lag time
    parser.add_argument(
        '--lag',
        help='Lag time for the MSM',
        type=int
    )

    args = parser.parse_args()

    return parser, args