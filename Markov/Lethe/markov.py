#! /usr/bin/env python3
"""
Markov building using PyEmma
"""

import pyemma
import tools
import numpy as np
import matplotlib.pyplot as plt

def create_feat(pdb,pair_indices):
    """Create a PyEmma featurizer

    Parameters
    ----------
    pdb : str
        Path to a reference topology .pdb file
    pair_indices : list
        List containing the indices of the pair to compute distances

    Returns
    -------
    feat : PyEmma object
        PyEmma featurizer
    """
    #Create the feat and add the right coordinates
    feat = pyemma.coordinates.featurizer(pdb)
    feat.add_distances(indices = pair_indices, periodic=True, indices2=None)
    print(f'PyEmma feat description:\n{feat.describe()}')
    return feat

def load_data(traj, feat):
    """Load all the given trajectories in PyEmma

    Parameters
    ----------
    traj : list
        List that contain the path of all the trajectories to analyze
    feat : PyEmma object
        PyEmma featurizer

    Returns
    -------
    data : List
        List of all the value of the feat for each snapshot
    """
    data = pyemma.coordinates.load(traj, features=feat,stride=1)
    
    print('Lengths (number of trajectories ):', len(data))
    print('Shape of elements (for each trajectories, number of timestep, number of features):', data[0].shape)
    
    return data

def plot_feat_hist(data, feat,display=False,save=False,outdir=''):
    """Make a histogram plot of the feat values

    Parameters
    ----------
    data : List
        List of all the value of the feat for each snapshot
    feat : PyEmma object
        PyEmma featurizer
    display : bool, optional
        Display or not the plots, by default False
    save : bool, optional
        Save or not the plots, by default False
    outdir : str, optional
        Path where to save the files, by default ''
    """
    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=feat,ignore_dim_warning=True)
    
    fig.set_figwidth(10)
    fig.set_figheight(5)
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/feat_hist.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()

def plot_density_energy(data, T, pairNames, save=False, display=False, outdir=''):
    """Plot the density map and the energy map

    Parameters
    ----------
    data : List
        List of all the value of the feat for each snapshot
    T : float
        Temperature
    pairNames : list
        List containing the name of the pairs
    save : bool, optional
        Save or not the plots, by default False
    display : bool, optional
        Display or not the plots, by default False
    outdir : str, optional
        Path where to save the files, by default ''
    """
    data_concatenated = np.concatenate(data)
    
    # Density plot
    fig, axes = plt.subplots(1, 1, figsize=(6, 4), sharex=True, sharey=True)    

    pyemma.plots.plot_density(*data_concatenated.T[0:2],ax=axes)
    
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')

    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/data_density.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()

    # Energy plot
    fig, axes = plt.subplots(1, 1, figsize=(6, 4), sharex=True, sharey=True)
    kT = tools.get_kT(T)
    fig,axes = pyemma.plots.plot_free_energy(*data_concatenated.T[0:2],
                                            ax=axes,
                                            kT=kT,
                                            cbar_label='free energy / kJ.mol-1')  
    
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')
    
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/data_free_energy_direct_from_density.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()

def pca_reduction(data,T,save=False,display=False,outdir=''):
    pca = pyemma.coordinates.pca(data,dim=2)
    pca_concatenated = np.concatenate(pca.get_output())
    fig,ax  =pyemma.plots.plot_feature_histograms(pca_concatenated,ignore_dim_warning=True)
    if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/pca_histogram.pdf', dpi=300,bbox_inches="tight")
    if display:
        plt.show()
    

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

    tica = pyemma.coordinates.tica(data, dim=2, lag=lag)
    tica_concatenated = np.concatenate(tica.get_output())

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

def clustering(reduction,method,k,stride,save=False,display=False,outdir=''):
    if method == 'kmeans':
        cluster = pyemma.coordinates.cluster_kmeans(reduction, k=k, stride=stride)
    if method == 'regspace':
        cluster = pyemma.coordinates.cluster_regspace(reduction,k=k,stride=stride)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    reduction_concatenated = np.concatenate(reduction)
    pyemma.plots.plot_feature_histograms(
        reduction_concatenated,
        ['IC {}'.format(i + 1) for i in range(reduction.dimension())],
        ax=axes[0])
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

    return cluster


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
    