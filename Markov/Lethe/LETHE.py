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
    # Get the args
    parser, args = LETHEparser.parsing()
    # Handle commun errors
    LETHEparser.LETHE_handle_error(parser,args)

    file_list = args.files
    pairNames = args.distances
    pdb = args.topology
    
    # Convert name in indices
    pair_indices = tools.create_pairIndices_from_pairNames(
        pdb,pairNames
        )
    # Create a feat
    feat = load_feat.create_feat(
        pdb,pair_indices
        )
    # Load the data
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
        # Histogram plot of the feat
        if 'feat_hist' in args.plot:
            load_feat.plot_feat_hist(
                data,feat,
                display=display,
                save=save,
                outdir=outdir
                )
        # Density energy map of the feat
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
    # PCA reduction
    if args.reduction == 'pca':
        red = dimension_reduction.pca_reduction(
            data=data,
            T=args.T,
            save=save,
            display=display,
            outdir=outdir
            )
    
    # TICA reduction
    if args.reduction == 'tica':
        red = dimension_reduction.tica_reduction(
            data=data,
            T=args.T,
            lag=args.tica_lag,
            save=save,
            display=display,
            outdir=outdir
            )
    
    # No reduction, raw data    
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
            
            # Plot the clustering result
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
    # ITS analysis
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
    
    # ITS analysis as a function of the number of cluster
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
    
    # CK test
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
    
    # Plot the stationary distribution
    if 'stationary' in args.plot:
        markov_analysis.plot_stationary(
            msm=msm,
            cluster=cluster,
            data=red,
            display=display,
            save=save,
            outdir=outdir
            )
    # Plot the first 6 eigenvectors
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
    
    # Display the stationary probabilities
    if args.pcca:
        pcca.stationary_prob(
            msm=msm,
            nstate=args.state
            )
        
        # Plot the metastable membership
        if 'metastable_membership' in args.plot:
            pcca.plot_metastable_membership(
                msm=msm,
                nstate=args.state,
                data=red,
                cluster=cluster,
                display=display,
                save=save,
                outdir=outdir
            )
        
        # Compute PCCA and TPT
        (
            metastable_traj,
            highest_membership,
            coarse_state_centers
        ) = pcca.concatenate(
            msm=msm,
            cluster=cluster
            )
        
        mfpt, inverse_mfpt = pcca.get_mfpt(
            msm=msm,
            nstates=args.state
            )
        
        # Plot MFPT
        if 'mfpt' in args.plot:
            pcca.plot_mftp(
                data=red,
                nstates=args.state,
                mfpt=mfpt,
                inverse_mfpt=inverse_mfpt,
                metastable_traj=metastable_traj,
                coarse_state_centers=coarse_state_centers,
                display=display,
                save=save,
                outdir=outdir
                )

        # Compute the flux for the committor    
        flux, cgflux = pcca.tpt(
            msm=msm,
            state=args.state_path,

        )

        # Plot the committor
        if 'committor' in args.plot:
            pcca.plot_committor_tpt(
                data=red,
                cluster=cluster,
                flux=flux,
                state=args.state_path,
                cgflux=cgflux,
                coarse_state_centers=coarse_state_centers,
                nstates=args.state,
                display=display,
                outdir=outdir,
                save=save
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