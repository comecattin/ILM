#! /usr/bin/env python3
"""Markov State Model validation"""

import pyemma
from dimension_reduction import *
from load_feat import *
from markov_analysis import *

def implied_time_scale(cluster,lags,nits):
    """Compute the implied time scale (ITS)

    Parameters
    ----------
    cluster : pyemma.cluster
        Cluster to make the analysis
    lags : list
        List if lag times for the MSM
    nits : int
        Number of iteration

    Returns
    -------
    its : pyemma.its
        Implied time scale
    """

    its = pyemma.msm.its(cluster.dtrajs,lags=lags,nits=nits,errors='bayes')

    return its

def plot_its(its,data, cluster, save=False, display=False,outdir=''):
    """Plot the ITS validation

    Parameters
    ----------
    its : pyemma.its
        Implied time scale
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

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    
    # Data were not reduced using TICA or PCA
    if type(data) == list:
        data_concatenated = np.concatenate(data)
    # Data reduced
    else:
        data_concatenated = np.concatenate(data.get_output())
    
    # Histogram plot
    pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=['Feat 1', 'Feat 2'], ax=axes[0])
    # Density plot
    pyemma.plots.plot_density(*data_concatenated.T[0:2], ax=axes[1], cbar=False, alpha=0.1)
    axes[1].scatter(*cluster.clustercenters.T[0:2], s=15, c='C1')
    axes[1].set_xlabel('Feat 1')
    axes[1].set_ylabel('Feat 2')
    # ITS plot
    pyemma.plots.plot_implied_timescales(its, ax=axes[2], units='steps')
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/its_validation.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()


def cluster_its(data,lags,nits, k_list,save=False,display=False,outdir=''):
    """ITS validation as a function of the cluster numbers

    Parameters
    ----------
    data : pyemma.load
        Data loaded from pyemma loader
    lag : int
        Lag time for the MSM
    nits : int
        Number of iteration
    k_list : list
        List of the different number of cluster size to test
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

    # Data were not reduced
    if type(data) == list:
        data_concatenated = np.concatenate(data)
    # Data were reduced 
    else:
        data_concatenated = np.concatenate(data.get_output())
    
    fig, axes = plt.subplots(2, len(k_list))
    
    for i, k in enumerate(k_list):
        # Loop over the number of cluster provided    
        cluster = pyemma.coordinates.cluster_kmeans(data, k=k, max_iter=500, stride=10)
        
        # Density plot
        pyemma.plots.plot_density(*data_concatenated.T[0:2], ax=axes[0, i], cbar=False, alpha=0.1)
        
        axes[0, i].scatter(*cluster.clustercenters.T[0:2], s=15, c='C1')
        axes[0, i].set_xlabel('Feat 1')
        axes[0, i].set_ylabel('Feat 2')
        axes[0, i].set_title('k = {} centers'.format(k))
        
        # ITS plot
        pyemma.plots.plot_implied_timescales(
            pyemma.msm.its(cluster.dtrajs, lags=lags, nits=nits, errors='bayes'),
            ax=axes[1, i], units='step')
        
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/its_cluster.pdf',dpi=300,bbox_inches='tight')

    if display:
        plt.show()

    
def cktest(msm,
           stable_state,
           save=False,
           display=False,
           outdir=''):
    """CK test

    Parameters
    ----------
    msm : pyemma.msm
        MSM or bayesian MSM
    stable_state : int
        Number of meta stable states to consider
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
    
    # Plot of the CK test
    pyemma.plots.plot_cktest(msm.cktest(stable_state), units='step')

    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/cktest.pdf',dpi=300,bbox_inches='tight')

    if display:
        plt.show()



if __name__ == '__main__':

     # Path
    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    outdir='/home/ccattin/Documents/Code/outputs'
    # Feat
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    # Parameters
    save = False
    display = False
    T = 300
    dim=2
    lag=1000
    nits = 4
    lags=[1, 2, 5, 10, 20, 50]
    k_list=[20,50,100]
    stable_state = 4

    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

    tica = tica_reduction(
        data=data,
        lag=lag,
        dim=dim
        )
    cluster = clustering(
        reduction=tica,
        method='kmeans',
        k=200,
        stride=1,
        )
    
    its = implied_time_scale(cluster=cluster,
                             lags=lags,
                             nits=nits)
    
    plot_its(its=its,
             data=tica,
             cluster=cluster,
             save=save,
             display=display,
             outdir=outdir)
    

    cluster_its(data=tica,
                lags=lags,
                nits=nits,
                k_list=k_list,
                save=save,
                display=display,
                outdir=outdir)
    
    msm = create_msm(cluster=cluster,
                     lag=lag,
                     error=True)
    
    cktest(msm=msm,
           stable_state=stable_state,
           display=display,
           outdir=outdir,
           save=save)
    
    