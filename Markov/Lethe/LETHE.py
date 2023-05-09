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

    LETHEparser.LETHE_handle_error(parser,args)

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
            T = float(args.T)
            markov.plot_density_energy(data=data,
                                       T=T,
                                       pairNames=pairNames,
                                       save=save,
                                       display=display,
                                       outdir=outdir)
            
    # Dimension reduction
    if args.pca:
        red = markov.pca_reduction(data=data,
                                   T=T,
                                   save=save,
                                   display=display,
                                   outdir=outdir)
    
    if args.tica:
        lag = args.lag
        red = markov.tica_reduction(data=data,
                                     T=T,
                                     lag=lag,
                                     save=save,
                                     display=display,
                                     outdir=outdir)
    # Clustering
    if args.cluster:

        if args.stride:
            stride = args.stride
        elif not args.stride:
            stride = 1
        
        cluster = markov.clustering(reduction=red,
                                    method=args.cluster,
                                    k=args.cluster_number,
                                    stride=stride,
                                    save=save,
                                    display=display,
                                    outdir=outdir)


if __name__ == '__main__':
    main()