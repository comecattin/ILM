#! /usr/bin/env python3
"""
Markov building using PyEmma
"""


if __name__ == "__main__":
    # Path
    pdb = "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb"
    traj = [
        "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc",
        "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc",
    ]
    outdir = "/home/ccattin/Documents/Code/outputs"
    # Feat
    pairNames = ["64_CA-130_CA", "119_CA-24_CA"]
    # Parameters
    save = True
    display = True
    T = 300
    lag = 1000

    pair_indices = tools.create_pairIndices_from_pairNames(pdb, pairNames)
    feat = create_feat(pdb, pair_indices)
    data = load_data(traj, feat)

    plot_feat_hist(data, feat, display=display, save=save, outdir=outdir)
    plot_density_energy(
        data=data, T=T, pairNames=pairNames, display=display, save=save, outdir=outdir
    )
    pca = pca_reduction(data=data, T=T, save=save, display=display, outdir=outdir)
    tica = tica_reduction(
        data=data, lag=lag, T=T, save=save, display=display, outdir=outdir
    )
    cluster = clustering(
        reduction=tica,
        method="kmeans",
        k=200,
        stride=1,
        save=save,
        display=display,
        outdir=outdir,
    )
