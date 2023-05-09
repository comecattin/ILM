import pyemma
import tools
from load_feat import *
import numpy as np
import matplotlib.pyplot as plt

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
        cluster = pyemma.coordinates.cluster_kmeans(reduction, k=k, stride=stride,max_iter=200)
    if method == 'regspace':
        cluster = pyemma.coordinates.cluster_regspace(reduction,k=k,stride=stride)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    reduction_concatenated = np.concatenate(reduction.get_output())
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