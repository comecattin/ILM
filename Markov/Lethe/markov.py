#! /usr/bin/env python3
"""
Markov building using PyEmma
"""

import pyemma
import tools
import numpy as np
import matplotlib.pyplot as plt

def create_feat(pdb,pair_indices):
    feat = pyemma.coordinates.featurizer(pdb)
    feat.add_distances(indices = pair_indices, periodic=True, indices2=None)
    print(f'PyEmma feat description:\n{feat.describe()}')
    return feat

def load_data(traj, feat):
    data = pyemma.coordinates.load(traj, features=feat,stride=1)
    print('Lengths (number of trajectories ):', len(data))
    print('Shape of elements (for each trajectories, number of timestep, number of features):', data[0].shape)
    return data

def plot_feat_hist(data, feat,display=False,save=False,outdir=''):
    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=feat,ignore_dim_warning=True)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    if save:
        plt.savefig(f'{outdir}/feat_hist.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()

def plot_density_energy(data, kT, pairNames, save=False, display=False, outdir=''):

    fig, axes = plt.subplots(1, 1, figsize=(6, 4), sharex=True, sharey=True)
    
    data_concatenated = np.concatenate(data)

    pyemma.plots.plot_density(*data_concatenated.T[0:2],ax=axes)
    
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')

    if save:
        plt.savefig(f'{outdir}/data_density.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()

    fig, axes = plt.subplots(1, 1, figsize=(6, 4), sharex=True, sharey=True)

    fig,axes = pyemma.plots.plot_free_energy(*data_concatenated.T[0:2],
                                            ax=axes,
                                            kT=kT,
                                            cbar_label='free energy / kJ.mol-1')  
    
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')
    
    if save:
        plt.savefig(f'{outdir}/data_free_energy_direct_from_density.pdf', dpi=300, bbox_inches="tight")
    if display:
        plt.show()






if __name__ == '__main__' :

    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    outdir='/home/ccattin/Documents/Code/outputs'
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    
    save = True
    display = True

    T = 300
    kT = tools.get_kT(T)

    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

    plot_feat_hist(data,feat,
                   display=display,
                   save=save,
                   outdir=outdir)
    plot_density_energy(data=data,
                        kT=kT,
                        pairNames=pairNames,
                        display=display,
                        save=save,
                        outdir=outdir
                        )