#! /usr/bin/env python3
"""
Markov building using PyEmma
"""

import pyemma
import tools

def create_feat(pdb,pair_indices):
    feat = pyemma.coordinates.featurizer(pdb)
    feat.add_distances(indices = pair_indices, periodic=True, indices2=None)
    print(f'PyEmma feat description:\n{feat.describe()}')
    return feat

if __name__ == '__main__' :

    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)