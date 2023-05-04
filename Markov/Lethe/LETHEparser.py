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
    parser.add_argument(
        '-d',
        '--distances',
        nargs='+',
        help='List of the pair names'
    )
    parser.add_argument(
        '-t',
        '--topology',
        required=True,
        help='Reference .pdb file'
    )
    parser.add_argument(
        '-p',
        '--plot',
        nargs='+',
        help='Plot wanted (feat_hist, density_energy)'
    )
    parser.add_argument(
        '--no-plot',
        action='store_true',
        help='Do not display plots'
    )
    parser.add_argument(
        '-o',
        '--outdir',
        help='Path to save the plots'
    )
    parser.add_argument(
        '--T',
        help='Temperature of the system'
    )
    parser.add_argument(
        '--pca',
        help='Do a PCA dimension reduction',
        action='store_true'
    )
    parser.add_argument(
        '--tica',
        help='Do a TICA dimension reduction',
        action='store_true'
    )
    parser.add_argument(
        '--lag',
        help='Lag time for the MSM',
        type=int
    )

    args = parser.parse_args()

    return parser, args