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

    #Handle plots
    if args.plot:
        # Parameters
        display = not args.no_plot
        if args.outdir:
            save=True
            outdir = args.outdir
        if not args.outdir:
            save=False
            outdir=''
        
        #Plots
        if 'feat_hist' in args.plot:
            markov.plot_feat_hist(data,feat,
                                  display=display,
                                  save=save,
                                  outdir=outdir)
        if 'density_energy' in args.plot:
            if not args.T:
                parser.error('Please provide a temperature')
            T = float(args.T)
            markov.plot_density_energy(data=data,
                                       T=T,
                                       pairNames=pairNames,
                                       save=save,
                                       display=display,
                                       outdir=outdir)
            
    # Dimension reduction
    if args.pca:
        if args.no_plot:
            raise NotImplementedError('Please use the plot option')
        pca = markov.pca_reduction(data=data,
                                   T=T,
                                   save=save,
                                   display=display,
                                   outdir=outdir)

if __name__ == '__main__':
    main()