#! /usr/bin/env python3
"""LETHE automatizes the MSM analysis"""

import aesthetic
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

    aesthetic.header()

    # ====PARSING====#
    # Get the args
    parser, args = LETHEparser.parsing()
    print("Option given")
    print(vars(args))

    # Handle commun errors
    LETHEparser.LETHE_handle_error(parser, args)

    # ====FEAT====#
    aesthetic.feat()
    # Create a feat
    feat = load_feat.create_feat(args.topology)

    if args.distances:
        # Convert name in indices
        pair_indices = tools.create_pairIndices_from_pairNames(
            args.topology, args.distances
        )
        feat = load_feat.feat_atom_distances(feat,pair_indices)

    if args.residue:
        # Get the indices
        pair_indices = tools.create_pairIndices_from_indices(
            args.residue
        )
        feat = load_feat.feat_residue_midist(feat,pair_indices)
        
    if args.feat_txt:
        pair_indices = tools.read_feat_from_txt(
            file_path = args.feat_txt[0],
            quality_max = int(args.feat_txt[1])
            )
        feat = load_feat.feat_atom_distances(feat,pair_indices)

    # Load the data
    data = load_feat.load_data(
        traj=args.files, feat=feat, stride=args.stride, ram=args.ram
    )

    # Compute the VAMP2 score
    if args.vamp_score:
        vamp = load_feat.vamp_score(data=data, dim=args.dim)

    # ====HANDLE INITIAL PLOTS====#
    if args.plot:
        # Parameters
        display = not args.no_plot
        if args.outdir:
            save = True
            outdir = args.outdir
        if not args.outdir:
            save = False
            outdir = ""

        # Plots
        # Histogram plot of the feat
        if "feat_hist" in args.plot:
            print("Rendering feat histogram...")
            load_feat.plot_feat_hist(
                data, feat, display=display, save=save, outdir=outdir
            )
        # Density energy map of the feat
        if "density_energy" in args.plot:
            print("Rendering density and energy map...")
            T = float(args.T)
            load_feat.plot_density_energy(
                data=data,
                T=T,
                pairNames=args.distances,
                save=save,
                display=display,
                outdir=outdir,
            )

    # ====DIMENSION REDUCTION====#
    aesthetic.dimension_reduction()

    # Get the lag time and the number of dimension
    if 'vamp_lag_dim' in args.plot:
        print("Rendering VAMP2 score as a function of lag time and dimension...")
        dimension_reduction.plot_lag_dim_vamp(
            lags=args.vamp_lags,
            data=data,
            dim=args.vamp_dim,
            save=save,
            outdir=outdir,
            display=display
        )

    # PCA reduction
    if args.reduction == "pca":
        print("PCA reduction...")
        red = dimension_reduction.pca_reduction(data=data, dim=args.dim)
        if "pca" in args.plot:
            print("Rendering PCA reduction plot...")
            dimension_reduction.plot_pca(
                pca=red,
                T=args.T,
                dim=args.dim,
                save=save,
                outdir=outdir,
                display=display,
            )

    # TICA reduction
    if args.reduction == "tica":
        print("TICA reduction...")
        red = dimension_reduction.tica_reduction(
            data=data, lag=args.tica_lag, dim=args.dim
        )
        if "tica" in args.plot:
            print("Rendering TICA reduction plot...")
            dimension_reduction.plot_tica(
                tica=red,
                T=args.T,
                dim=args.dim,
                display=display,
                save=save,
                outdir=outdir,
            )

    # VAMP reduction
    if args.reduction == "vamp":
        print("VAMP reduction...")
        red = dimension_reduction.vamp_reduction(
            data=data, lag=args.vamp_lag, dim=args.dim
        )
        if "vamp" in args.plot:
            print("Rendering VAMP reduction plot...")
            dimension_reduction.plot_vamp(
                vamp=red,
                T=args.T,
                dim=args.dim,
                save=save,
                display=display,
                outdir=outdir,
            )

    # No reduction, raw data
    if args.reduction == "none":
        print("No reduction, raw data taken")
        red = data

    # VAMP2 score as a function of the number of cluster
    if args.vamp_cluster:
        dimension_reduction.plot_vamp_cluster(
            n_clustercenters=args.vamp_cluster,
            lag=args.lag,
            data=red,
            save=save,
            display=display,
            outdir=outdir
        )


    # ====LOAD PREVIOUS MODEL====#
    if args.load:
        msm, cluster = tools.load_model(
            outdir=outdir, filename=args.load[0], model_name=args.load[1]
        )

    # =====DO NOT LOAD PREVIOUS MODEL====#
    if not args.load:
        # ====CLUSTERING====#

        if args.cluster:
            aesthetic.cluster()
            print("Clustering...")
            cluster = dimension_reduction.clustering(
                reduction=red, method=args.cluster, k=args.cluster_number, stride=1
            )

            # Plot the clustering result
            if "cluster" in args.plot:
                print("Rendering clustering plot...")
                dimension_reduction.clustering_plot(
                    reduction=red,
                    cluster=cluster,
                    save=save,
                    outdir=outdir,
                    display=display,
                )

        # ====CREATE MSM====#
        print("Creating MSM...")
        msm = markov_analysis.create_msm(
            cluster=cluster, lag=args.lag, error=args.confidence
        )

    
    # ====SAVE MSM====#
    if args.save:
        # Save the model
        tools.save_model(
            cluster=cluster,
            msm=msm,
            outdir=outdir,
            filename=args.save[0],
            model_name=args.save[1],
        )

    # ====MSM VALIDATION====#
    aesthetic.validation()
    # ITS analysis
    if args.its:
        print("ITS analysis...")
        its = validation.implied_time_scale(
            cluster=cluster, lags=args.its, nits=args.nits
        )

        if "its" in args.plot:
            validation.plot_its(
                its=its,
                data=red,
                cluster=cluster,
                save=save,
                display=display,
                outdir=outdir,
            )

    # ITS analysis as a function of the number of cluster
    if args.its_cluster:
        print("ITS cluster analysis...")
        validation.cluster_its(
            data=red,
            lags=args.its,
            nits=args.nits,
            k_list=args.its_cluster,
            stride=1,
            save=save,
            display=display,
            outdir=outdir,
        )

    # CK test
    if "cktest" in args.plot:
        print("CK Test...")
        validation.cktest(
            msm=msm, stable_state=args.state, display=display, outdir=outdir, save=save
        )

    # ====MSM ANALYSIS====#
    aesthetic.analysis()
    print("fraction of states used = {:f}".format(msm.active_state_fraction))
    print("fraction of counts used = {:f}".format(msm.active_count_fraction))

    # Plot the stationary distribution
    if "stationary" in args.plot:
        print("Rendering stationary plot...")
        markov_analysis.plot_stationary(
            msm=msm,
            cluster=cluster,
            data=red,
            display=display,
            save=save,
            outdir=outdir,
        )
    
    # Plot the eigenvalues
    if 'eigenvalues' in args.plot:
        print("Rendering eigenvalues plot...")
        markov_analysis.plot_eigenvalues(
            msm=msm,
            nvalues=args.nits,
            save=save,
            display=display,
            outdir=outdir
        )

    # Plot the first 6 eigenvectors
    if "eigenvectors" in args.plot:
        print("Rendering eigenvectors plot...")
        markov_analysis.plot_eigenvect(
            msm=msm,
            cluster=cluster,
            data=red,
            display=display,
            save=save,
            outdir=outdir,
        )

    # ====PCCA AND TPT====#

    # Display the stationary probabilities
    if args.pcca:
        aesthetic.pcca()
        print("Computing stationary probabilities...")
        pcca.stationary_prob(msm=msm, nstate=args.state)

        # Plot the metastable membership
        if "metastable_membership" in args.plot:
            print("Rendering metastable membership...")
            pcca.plot_metastable_membership(
                msm=msm,
                nstate=args.state,
                data=red,
                display=display,
                save=save,
                outdir=outdir,
            )

        # Compute PCCA and TPT
        print("Concatenating results...")
        (metastable_traj, highest_membership, coarse_state_centers) = pcca.concatenate(
            msm=msm, cluster=cluster
        )

        print("Computing MFPT")

        mfpt, inverse_mfpt = pcca.get_mfpt(msm=msm, nstates=args.state)

        # Plot MFPT
        if "mfpt" in args.plot:
            print("Rendering MFPT")
            pcca.plot_mftp(
                data=red,
                nstates=args.state,
                mfpt=mfpt,
                inverse_mfpt=inverse_mfpt,
                metastable_traj=metastable_traj,
                coarse_state_centers=coarse_state_centers,
                display=display,
                save=save,
                outdir=outdir,
            )

        # Compute the flux for the committor
        print("Computing flux...")
        flux, cgflux = pcca.tpt(
            msm=msm,
            state=args.state_path,
        )

        # Plot the committor
        if "committor" in args.plot:
            print("Rendering committor plot...")
            pcca.plot_committor_tpt(
                data=red,
                msm=msm,
                flux=flux,
                state=args.state_path,
                cgflux=cgflux,
                coarse_state_centers=coarse_state_centers,
                nstates=args.state,
                display=display,
                outdir=outdir,
                save=save,
            )

        # Render pdb sample files
        if args.pdb:
            print("Saving .pdb sample files...")
            pcca.sample_structures(
                msm=msm,
                number_of_sample=args.pdb,
                feat=feat,
                files=args.files,
                outdir=outdir
            )


if __name__ == "__main__":
    main()
