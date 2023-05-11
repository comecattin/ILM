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

def plot_metastable_membership(
          msm,
          nstate,
          data,
          cluster,
          display=False,
          save=False,
          outdir=''
          ):

     # Data without dimension reduction
     if type(data) == list:
          data_concatenated = np.concatenate(data)
     # Dimension reduction
     else:
          data_concatenated = np.concatenate(data.get_output())
          
     dtrajs_concatenated = np.concatenate(cluster.dtrajs)

     fig, axes = plt.subplots(1, nstate, figsize=(15, 3))
     for i, ax in enumerate(axes.flat):
          pyemma.plots.plot_contour(
               *data_concatenated.T,
               msm.metastable_distributions[i][dtrajs_concatenated],
               ax=ax,
               cmap='afmhot_r', 
               mask=True,
               cbar_label='Metastable distribution {}'.format(i + 1)
               )
          ax.set_xlabel('Feat 1')
     axes[0].set_ylabel('Feat 2')
     
     if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/metastable_membership.pdf',dpi=300,bbox_inches='tight')

     if display:
        plt.show()


def concatenate(msm,cluster):
     dtrajs_concatenated = np.concatenate(cluster.dtrajs)

     metastable_traj = msm.metastable_assignments[
         dtrajs_concatenated
         ]
     
     highest_membership = msm.metastable_distributions.argmax(1)
     
     coarse_state_centers = cluster.clustercenters[
         msm.active_set[highest_membership]
         ]
     
     return metastable_traj, highest_membership, coarse_state_centers


def get_mfpt(msm,nstates):
     mfpt = np.zeros((nstates, nstates))
     for i in range(nstates):
          for j in range(nstates):
               mfpt[i, j] = msm.mfpt(
                    msm.metastable_sets[i],
                    msm.metastable_sets[j])

     inverse_mfpt = np.zeros_like(mfpt)
     nz = mfpt.nonzero()
     inverse_mfpt[nz] = 1.0 / mfpt[nz]
     return mfpt, inverse_mfpt

def plot_mftp(
          data,
          nstates,
          mfpt,
          inverse_mfpt,
          display=False,
          save=False,
          outdir=''
          ):

     # Data without dimension reduction
     if type(data) == list:
          data_concatenated = np.concatenate(data)
     # Dimension reduction
     else:
          data_concatenated = np.concatenate(data.get_output())

     fig, ax = plt.subplots(figsize=(10, 7))

     _, _, misc = pyemma.plots.plot_state_map(
          *data_concatenated.T,
          metastable_traj,
          ax=ax,
          zorder=-1)
     # set state numbers 1 ... nstates
     misc['cbar'].set_ticklabels(range(1, nstates + 1))  

     pyemma.plots.plot_network(
          inverse_mfpt,
          pos=coarse_state_centers,
          figpadding=0,
          arrow_label_format='%.1f ps',
          arrow_labels=mfpt,
          size=12,
          show_frame=True,
          ax=ax)

     ax.set_xlabel('Feat 1')
     ax.set_ylabel('Feat 2')
     ax.set_xlim(
         min(data_concatenated[:,0]),
         max(data_concatenated[:,0])
         )
     ax.set_ylim(
         min(data_concatenated[:,1]),
         max(data_concatenated[:,1])
         )
     
     if save:
        if outdir=='':
            raise Exception('Please provide a directory to save the file')
        else:
            plt.savefig(f'{outdir}/metastable_membership.pdf',dpi=300,bbox_inches='tight')

     if display:
        plt.show()


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
    stable_state = 2

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

    stationary_prob(
        msm=msm,
        nstate=stable_state
        )
    plot_metastable_membership(
        msm=msm,
        nstate=stable_state,
        data=tica,
        cluster=cluster,
        display=display,
        save=save,
        outdir=outdir
    )

    (
        metastable_traj,
        highest_membership,
        coarse_state_centers
        ) = concatenate(
        msm=msm,
        cluster=cluster
    )
    print(
        metastable_traj,
        highest_membership,
        coarse_state_centers
        )
    
    mfpt, inverse_mfpt = get_mfpt(
         msm=msm,
         nstates=stable_state
    )
    print(mfpt,inverse_mfpt)


    plot_mftp(
        data=tica,
        nstates=stable_state,
        mfpt=mfpt,
        inverse_mfpt=inverse_mfpt,
        display=display,
        save=save,
        outdir=outdir
    )
