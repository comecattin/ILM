#! /usr/bin/env python3
"""Markov State Model analysis"""

import pyemma
import tools
from load_feat import *
from dimension_reduction import *
from validation import *

def create_msm(cluster,lag,dt_traj='1 ps',error=False):

    if error:
        msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs,
                                               lag=lag,
                                               dt_traj=dt_traj,
                                               conf=0.95)
        
    else:
        msm = pyemma.msm.estimate_markov_model(cluster.dtrajs,
                                           lag=lag,
                                           dt_traj=dt_traj)

    return msm




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
    lag=1000
    nits = 4
    lags=[1, 2, 5, 10, 20, 50]
    stable_state = 4

    pair_indices = tools.create_pairIndices_from_pairNames(pdb,pairNames)
    feat = create_feat(pdb,pair_indices)
    data = load_data(traj,feat)

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
    

    msm = create_msm(cluster=cluster,
           lag=lag,
           error=True)