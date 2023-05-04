#! /usr/bin/env python3
"""LETHE automatizes the MSM analysis"""

import markov
import tools
import LETHEparser

def main():
    """
    CLI main function
    """
    
    # Parsing
    parser, args = LETHEparser.parsing()

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
            parser.error('Please use the plot option')
        pca = markov.pca_reduction(data=data,
                                   T=T,
                                   save=save,
                                   display=display,
                                   outdir=outdir)
    
    if args.tica:
        if args.no_plot:
            parser.error('Please use the plot option')
        if not args.lag:
            parser.error('Please provide a lag time')
        lag = args.lag
        tica = markov.tica_reduction(data=data,
                                     T=T,
                                     lag=lag,
                                     save=save,
                                     display=display,
                                     outdir=outdir)

if __name__ == '__main__':
    main()