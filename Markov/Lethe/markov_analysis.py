#! /usr/bin/env python3
"""Markov State Model analysis"""

import pyemma
import tools
from load_feat import *
from dimension_reduction import *
from validation import *


def create_msm(cluster, lag, dt_traj="1 ps", error=False):
    """Create a MSM

    Parameters
    ----------
    cluster : pyemma.cluster
        Cluster to make the analysis
    lag : int
        Lag time for the MSM
    dt_traj : str, optional
        Axe legend, by default '1 ps'
    error : bool, optional
        Bayesian MSM or not, by default False

    Returns
    -------
    msm : pyemma.msm
        MSM or bayesian MSM
    """

    if error:
        msm = pyemma.msm.bayesian_markov_model(
            cluster.dtrajs, lag=lag, dt_traj=dt_traj, conf=0.95
        )

    else:
        msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=lag, dt_traj=dt_traj)

    return msm


def plot_stationary(msm, cluster, data, display=False, save=False, outdir=""):
    """Stationary plot

    Parameters
    ----------
    msm : pyemma.msm
        MSM or bayesian MSM
    cluster : pyemma.cluster
        Cluster to make the analysis
    data : pyemma.load
        Data loaded from pyemma loader
    save : bool, optional
        Save or not the plot, by default False
    display : bool, optional
        Display or not the plot, by default False
    outdir : str, optional
        Output directory to save the plot, by default ''

    Raises
    ------
    Exception
        Provide a directory to save the file
    """

    # Data without dimension reduction
    if type(data) == list:
        data_concatenated = np.concatenate(data)
    # Dimension reduction
    else:
        data_concatenated = np.concatenate(data.get_output())

    dtrajs_concatenated = np.concatenate(msm.dtrajs_active)

    fig, ax, misc = pyemma.plots.plot_contour(
        *data_concatenated.T[0:2],
        msm.pi[dtrajs_concatenated],
        cbar_label="Stationary distribution",
        method="nearest",
        mask=True,
    )
    ax.scatter(*cluster.clustercenters.T[0:2], s=15, c="C1")
    ax.set_xlabel("Feat 1")
    ax.set_ylabel("Feat 2")

    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(
                f"{outdir}/stationary_distribution.pdf", dpi=300, bbox_inches="tight"
            )

    if display:
        plt.show()


def plot_eigenvect(msm, data, cluster, display=False, save=False, outdir=""):
    """Plot the first 6 eigen vectors

    Parameters
    ----------
    msm : pyemma.msm
        MSM or bayesian MSM
    cluster : pyemma.cluster
        Cluster to make the analysis
    data : pyemma.load
        Data loaded from pyemma loader
    save : bool, optional
        Save or not the plot, by default False
    display : bool, optional
        Display or not the plot, by default False
    outdir : str, optional
        Output directory to save the plot, by default ''

    Raises
    ------
    Exception
        Provide a directory to save the file
    """

    # Data without dimension reduction
    if type(data) == list:
        data_concatenated = np.concatenate(data)
    # Dimension reduced
    else:
        data_concatenated = np.concatenate(data.get_output())

    dtrajs_concatenated = np.concatenate(msm.dtrajs_active)

    eigvec = msm.eigenvectors_right()

    # Confirmation print
    print(
        "first eigenvector is one: {} (min={}, max={})".format(
            np.allclose(eigvec[:, 0], 1, atol=1e-15),
            eigvec[:, 0].min(),
            eigvec[:, 0].max(),
        )
    )

    fig, axes = plt.subplots(2, 3, figsize=(17, 6))

    for i, ax in enumerate(axes.flat):
        pyemma.plots.plot_contour(
            *data_concatenated.T[0:2],
            eigvec[dtrajs_concatenated, i + 1],
            ax=ax,
            cmap="PiYG",
            cbar_label="{}. right eigenvector".format(i + 2),
            mask=True,
        )
        #ax.scatter(*cluster.clustercenters.T[0:2], s=15, c="C1")
        ax.set_xlabel("Feat 1")
        ax.set_ylabel("Feat 2")

    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(f"{outdir}/eigen_vect.pdf", dpi=300, bbox_inches="tight")

    if display:
        plt.show()

def plot_eigenvalues(msm, nvalues,save=False,display=False,outdir=''):
    
    eigenvalues = np.abs(msm.eigenvalues(nvalues))

    fig, ax = plt.subplots()

    ax.bar(
        np.arange(nvalues),
        eigenvalues,
        width=0.05,
        color='grey',
        zorder=0
        )
    ax.scatter(
        np.arange(nvalues),
        eigenvalues,
        zorder=1,
        s=50
        )
    ax.plot(
        np.arange(nvalues),
        np.ones(nvalues),
        '--',
        color='lightgrey'
    )

    ax.set_xlabel("Index")
    ax.set_ylabel("Eigenvalue")
    ax.set_ylim(0, 1.1)
    
    if save:
        if outdir == "":
            raise Exception("Please provide a directory to save the file")
        else:
            plt.savefig(f"{outdir}/eigen_values.pdf", dpi=300, bbox_inches="tight")

    if display:
        plt.show()



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
    save = False
    display = False
    T = 300
    dim = 2
    lag = 1000
    nits = 4
    lags = [1, 2, 5, 10, 20, 50]
    stable_state = 4

    # Save
    filename = "test.pyemma"
    model_name = "GS01_GS02"

    pair_indices = tools.create_pairIndices_from_pairNames(pdb, pairNames)
    feat = create_feat(pdb)
    feat = feat_atom_distances(feat,pair_indices)
    data = load_data(traj, feat,stride=2,ram=True)

    tica = tica_reduction(data=data, lag=lag, dim=dim)

    cluster = clustering(reduction=tica, method="kmeans", k=200, stride=1)

    msm = create_msm(cluster=cluster, lag=lag, error=True)

    plot_stationary(
        msm=msm, data=tica, cluster=cluster, display=display, save=save, outdir=outdir
    )

    plot_eigenvect(
        msm=msm, data=tica, cluster=cluster, display=display, outdir=outdir, save=save
    )

    tools.save_model(
        cluster=cluster,
        msm=msm,
        outdir=outdir,
        filename=filename,
        model_name=model_name,
    )

    plot_eigenvalues(
        msm,
        6,
        save=save,
        display=display,
        outdir=outdir)