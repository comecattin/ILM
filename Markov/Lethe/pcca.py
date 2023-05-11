#! /usr/bin/env python3
"""Markov State Model PCCA++ and TPT analysis"""

import pyemma
import tools
from load_feat import *
from dimension_reduction import *
from validation import *
from markov_analysis import *

def stationary_prob(msm,nstate):
     msm.pcca(nstate)
     print('Sationary probabilities of the metastable sets')
     for i, s in enumerate(msm.metastable_sets):
          print('π_{} = {:f}'.format(i + 1, msm.pi[s].sum()))
























if __name__ == '__main__':

     # Path
    pdb = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    traj = ['/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS01_md_all_fitBB_protonly.xtc','/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS02_md_all_fitBB_protonly.xtc']
    outdir='/home/ccattin/Documents/Code/outputs'
    filename='test.pyemma'
    model_name = 'GS01_GS02'
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
    
#     cluster = clustering(reduction=tica,
#                          method='kmeans',
#                          k=200,
#                          stride=1)
    
#     msm = create_msm(cluster=cluster,
#            lag=lag,
#            error=False)
    
#     tools.save_model(
#         cluster=cluster,
#         msm=msm,
#         outdir=outdir,
#         model_name=model_name,
#         filename=filename
#         )
    msm, cluster = tools.load_model(
         outdir=outdir,
         filename=filename,
         model_name=model_name
    )

    stationary_prob(msm=msm,
                    nstate=stable_state)
