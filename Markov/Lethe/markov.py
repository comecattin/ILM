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

def plot_feat_hist(data, feat):
    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=feat,ignore_dim_warning=True)
    fig.set_figwidth(10)
    fig.set_figheight(5)
    fig.tight_layout()
    plt.show()

if __name__ == '__main__' :

    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    
    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

    plot_feat_hist(data,feat)