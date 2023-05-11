#! /usr/bin/env python3
"""LETHE automatizes the MSM analysis"""

import load_feat
import dimension_reduction
import tools
import LETHEparser
import markov_analysis
import validation
import pcca

def main():
    """
    CLI main function
    """
    
    #====PARSING====#
    parser, args = LETHEparser.parsing()

    LETHEparser.LETHE_handle_error(parser,args)

    file_list = args.files
    pairNames = args.distances
    pdb = args.topology
    
    pair_indices = tools.create_pairIndices_from_pairNames(
        pdb,pairNames
        )

    feat = load_feat.create_feat(
        pdb,pair_indices
        )

    data = load_feat.load_data(
        traj=file_list,feat=feat
        )

    #====HANDLE INITIAL PLOTS====#
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
            load_feat.plot_feat_hist(
                data,feat,
                display=display,
                save=save,
                outdir=outdir
                )
        if 'density_energy' in args.plot:
            T = float(args.T)
            load_feat.plot_density_energy(
                data=data,
                T=T,
                pairNames=pairNames,
                save=save,
                display=display,
                outdir=outdir
                )
            
    #====DIMENSION REDUCTION====#
    if args.reduction == 'pca':
        red = dimension_reduction.pca_reduction(
            data=data,
            T=args.T,
            save=save,
            display=display,
            outdir=outdir
            )
    
    if args.reduction == 'tica':
        lag = args.lag
        red = dimension_reduction.tica_reduction(
            data=data,
            T=args.T,
            lag=lag,
            save=save,
            display=display,
            outdir=outdir
            )
        
    if args.reduction == 'none' :
        red = data


    #====LOAD PREVIOUS MODEL====#
    if args.load:
        msm, cluster = tools.load_model(
            outdir=outdir,
            filename=args.load[0],
            model_name=args.load[1])

    #=====DO NOT LOAD PREVIOUS MODEL====#    
    if not args.load:
        
        #====CLUSTERING====#
        if args.cluster:

            if args.stride:
                stride = args.stride
            elif not args.stride:
                stride = 1
            
            cluster = dimension_reduction.clustering(
                reduction=red,
                method=args.cluster,
                k=args.cluster_number,
                stride=stride
                )
            
            if 'cluster' in args.plot:
                dimension_reduction.clustering_plot(
                    reduction=red,
                    cluster=cluster,
                    save=save,
                    outdir=outdir,
                    display=display
                    )
                
        #====CREATE MSM====#
        msm = markov_analysis.create_msm(
            cluster=cluster,
            lag=args.lag,
            error=args.confidence
            )
        

    #====MSM VALIDATION====#
    if args.its:
        its = validation.implied_time_scale(
            cluster=cluster,
            lags=args.its,
            nits=args.nits)
        
        if 'its' in args.plot:
            validation.plot_its(
                its=its,
                data=red,
                cluster=cluster,
                save=save,
                display=display,
                outdir=outdir
                )
    
    if args.its_cluster:
        validation.cluster_its(
            data=red,
            lags=args.its,
            nits=args.nits,
            k_list=args.its_cluster,
            save=save,
            display=display,
            outdir=outdir
            )
    
    if 'cktest' in args.plot:
        validation.cktest(
            msm=msm,
            stable_state=args.state,
            display=display,
            outdir=outdir,
            save=save
            )
    
    #====MSM ANALYSIS====#
    print(
        'fraction of states used = {:f}'.format(
        msm.active_state_fraction
        )
        )
    print(
        'fraction of counts used = {:f}'.format(
        msm.active_count_fraction
        )
        )
    
    if 'stationary' in args.plot:
        markov_analysis.plot_stationary(
            msm=msm,
            cluster=cluster,
            data=red,
            display=display,
            save=save,
            outdir=outdir
            )
    
    if 'eigenvectors' in args.plot:
        markov_analysis.plot_eigenvect(
            msm=msm,
            cluster=cluster,
            data=red,
            display=display,
            save=save,
            outdir=outdir
        )

    #====PCCA AND TPT====#
    if args.stationary_prob:
        pcca.stationary_prob(
            msm=msm,
            nstate=args.state
            )
    
    #====SAVE MSM====#
    if args.save:
        # Save the model
        tools.save_model(
            cluster=cluster,
            msm=msm,
            outdir=outdir,
            filename=args.save[0],
            model_name=args.save[1]
        )




if __name__ == '__main__':
    main()