#! /usr/bin/env python3
"""LETHE automatizes the MSM analysis"""

import argparse
import markov
import tools

def main():
    """
    CLI main function
    """
    
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
        help='Plot wanted (feat_hist)'
    )
    
    args = parser.parse_args()

    #Handles error
    if not args.files:
        parser.error("No trajectories file given")
    if not args.distances:
        parser.error("No distances given")

    file_list = args.files
    pairNames = args.distances
    pdb = args.topology
    
    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)

    feat = markov.create_feat(pdb,pair_indices)

    data = markov.load_data(traj=file_list,feat=feat)

    if args.plot:
        if 'feat_hist' in args.plot:
            markov.plot_feat_hist(data,feat)

if __name__ == '__main__':
    main()