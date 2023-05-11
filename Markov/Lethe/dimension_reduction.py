#! /usr/bin/env python3
"""Dimension reduction for faster analysis MSM"""

import pyemma
import tools
from load_feat import *
import numpy as np
import matplotlib.pyplot as plt

def pca_reduction(data,T,save=False,display=False,outdir=''):
    """Do a PCA dimension reduction and plot it

    Parameters
    ----------
    data : pyemma.load
        Data loaded from pyemma loader
    T : float
        Temperature
    save : bool, optional
        Save or not the plot, by default False
    display : bool, optional
        Display or not the plot, by default False
    outdir : str, optional
        Output directory to save the plot, by default ''

    Returns
    -------
    pca : pyemma.pca
        Data reduced using PCA

    Raises
    ------
    Exception
        Provide a directory to save the file
    """
    # Dimension reduction
    pca = pyemma.coordinates.pca(data,dim=2)
    # Concatenate
    pca_concatenated = np.concatenate(pca.get_output())
    
    # Histogramm plot of the new dimension
    fig,ax  =pyemma.plots.plot_feature_histograms(pca_concatenated,ignore_dim_warning=True)
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/pca_histogram.pdf', dpi=300,bbox_inches="tight")
    if display:
        plt.show()
    
    # Free energy plot of the new dimension
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    kT = tools.get_kT(T)
    pyemma.plots.plot_free_energy(*pca_concatenated.T[0:2],
                                            ax=axes,
                                            kT=kT,
                                            cbar_label='free energy / kJ.mol-1') 
    axes.set_xlabel('PC 1')
    axes.set_ylabel('PC 2')
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/pca_free_energy_direct_from_density.pdf', dpi=300,bbox_inches="tight")
    if display:
        plt.show()
    
    return pca

def tica_reduction(data,lag,T,save=False,display=False,outdir=''):
    """Do a TICA dimension reduction and plot it

    Parameters
    ----------
    data : pyemma.load
        Data loaded from pyemma loader
    lag : int
        Lag time
    T : float
        Temperature of the system
    save : bool, optional
        Save or not the plot, by default False
    display : bool, optional
        Display or not the plot, by default False
    outdir : str, optional
        Output directory to save the plot, by default ''

    Returns
    -------
    pca : pyemma.tica
        Data reduced using TICA

    Raises
    ------
    Exception
        Provide a directory to save the file
    """

    # TICA reduction
    tica = pyemma.coordinates.tica(data, dim=2, lag=lag)
    tica_concatenated = np.concatenate(tica.get_output())

    # Histogramm plot
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    pyemma.plots.plot_feature_histograms(
        tica_concatenated,
        ['TICA {}'.format(i + 1) for i in range(tica.dimension())],
        ax=axes,
        ignore_dim_warning=True)
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/tica_histogram.pdf',dpi=300,bbox_inches='tight')
    if display:
        plt.show()

    # Free energy plot
    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    kT=tools.get_kT(T)
    pyemma.plots.plot_free_energy(*tica_concatenated.T[0:2],
                                            ax=axes,
                                            kT=kT,cbar_label='free energy / kJ.mol-1') 
    axes.set_xlabel('TIC 1')
    axes.set_ylabel('TIC 2')
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/tica_free_energy_direct_from_density.pdf',dpi=300,bbox_inches='tight')
    if display:
        plt.show()
    
    return tica

def clustering(reduction,method,k,stride):
    """Clustering of the data

    Parameters
    ----------
    reduction : pyemma.data
        Data from pyemma loader or from PCA, TICA reduction
    method : str
        Method to use, 'kmeans' and 'regspace' are implemented
    k : int
        Number of cluster
    stride : int
        Stride

    Returns
    -------
    cluster : pyemma.cluster
        Data clustered
    """
    if method == 'kmeans':
        cluster = pyemma.coordinates.cluster_kmeans(reduction, k=k, stride=stride,max_iter=200)
    if method == 'regspace':
        cluster = pyemma.coordinates.cluster_regspace(reduction,k=k,stride=stride)
    
    return cluster

def clustering_plot(reduction,cluster,save=False,outdir='',display=False):
    """Plot of the clustering

    Parameters
    ----------
    reduction : pyemma.data
        Data from pyemma loader or from PCA, TICA reduction
    cluster : pyemma.cluster
        Data clustered
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

    # If no reduction has been performed
    if type(reduction) == list:
        fig, axe = plt.subplots(1,1,figsize=(10, 4))
        reduction_concatenated = np.concatenate(reduction)
        pyemma.plots.plot_density(*reduction_concatenated.T, cbar=False, alpha=0.1, logscale=True,ax=axe)
        axe.scatter(*cluster.clustercenters.T, s=15, c='C1')
        axe.set_xlabel('IC 1')
        axe.set_ylabel('IC 2')
    
    # Reduction has been performed
    else:
        fig, axes = plt.subplots(1,2,figsize=(10, 4))
        reduction_concatenated = np.concatenate(reduction.get_output())
        # Histogram plot
        pyemma.plots.plot_feature_histograms(
            reduction_concatenated,
            ['IC {}'.format(i + 1) for i in range(reduction.dimension())],
            ax=axes[0])
        # Density plot
        pyemma.plots.plot_density(*reduction_concatenated.T, ax=axes[1], cbar=False, alpha=0.1, logscale=True)
        axes[1].scatter(*cluster.clustercenters.T, s=15, c='C1')
        axes[1].set_xlabel('IC 1')
        axes[1].set_ylabel('IC 2')
        
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/cluster.pdf',dpi=300,bbox_inches='tight')
    if display:
        plt.show()


if __name__ == '__main__' :

    # Path
    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    outdir='/home/ccattin/Documents/Code/outputs'
    # Feat
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    # Parameters
    save = True
    display = True
    T = 300
    lag=1000

    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

    plot_feat_hist(data,feat,
                   display=display,
                   save=save,
                   outdir=outdir)
    plot_density_energy(data=data,
                        T=T,
                        pairNames=pairNames,
                        display=display,
                        save=save,
                        outdir=outdir
                        )
    pca = pca_reduction(data=data,
                        T=T,
                        save=save,
                        display=display,
                        outdir=outdir)
    tica = tica_reduction(data=data,
                          lag=lag,
                          T=T,
                          save=save,
                          display=display,
                          outdir=outdir)
    cluster = clustering(reduction=tica,
                         method='kmeans',
                         k=200,
                         stride=1,
                         save=save,
                         display=display,
                         outdir=outdir)