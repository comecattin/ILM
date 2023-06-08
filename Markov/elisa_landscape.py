#! /usr/bin/env python3
"""Plot of the Elisa's free energy landscape"""
import pyemma
import numpy as np
import matplotlib.pyplot as plt
import glob
from Lethe import tools

def create_feat(pdb,pair_indices,offset=11):
    feat = pyemma.coordinates.featurizer(pdb)
    pair_indices = pair_indices - offset
    feat.add_residue_mindist(residue_pairs=pair_indices)

    return feat

def load_data(feat,traj):
    data = pyemma.coordinates.load(traj, features=feat, stride=1)
    return data

def clustering(data,k):
    cluster = pyemma.coordinates.cluster_kmeans(
            data, k=k, stride=1, max_iter=200
        )

    return cluster

def create_msm(cluster,lag):
    msm = pyemma.msm.bayesian_markov_model(
            cluster.dtrajs, lag=lag, conf=0.95
        )
    return msm

def plot_reweighted_free_energy(data,msm,T,path):
    data_concatenated = np.concatenate(data)
    dtrajs_concatenated = np.concatenate(msm.dtrajs_active)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)

    # Stationary distribution
    pyemma.plots.plot_contour(
        *data_concatenated[:, [0,1]].T,
        msm.pi[dtrajs_concatenated],
        ax=axes[0],
        mask=True,
        cbar_label='Stationary distribution')
    
    # Free energy
    kT = tools.get_kT(T)
    pyemma.plots.plot_free_energy(
        *data_concatenated[:, [0,1]].T,
        weights=np.concatenate(msm.trajectory_weights()),
        ax=axes[1],
        legacy=False,
        kT=kT,
        cbar_label="free energy / kJ.mol-1"
    )
    for ax in axes.flat:
        ax.set_xlabel(f'IC {1}')
    axes[0].set_ylabel(f'IC {2}')
    #axes[0].set_title('Stationary distribution', fontweight='bold')
    #axes[1].set_title('Reweighted free energy surface', fontweight='bold')
    
    axes[1].scatter(path[:,0],path[:,1])
   
    plt.savefig(f"reweighted_free_energy.pdf", dpi=300, bbox_inches="tight")

    plt.show()

if __name__ == '__main__':

    pdb = "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb"
    GS = glob.glob("/home/ccattin/Documents/data_elisa/MD_Gromacs_NMR_Structures_GS_ER-AF-CNSrefined-19nov21_mutant-R60A/*/md_1-10_every1ns_fitBB_protonly.xtc")
    ES = glob.glob("/home/ccattin/Documents/data_elisa/MD_Gromacs_NMR_Structures_ES_ER-AF-CNSrefined-8nove21_R46A/*/md_1-10_evert1ns_fitBB_protonly.xtc")
    files = GS + ES
    path_amber14_TIP3P = np.array(
        [
            [0.34484E+02,0.78440E+01],
            [0.32540E+02,0.82295E+01],
            [0.30562E+02,0.82573E+01],
            [0.28942E+02,0.93670E+01],
            [0.28564E+02,0.11309E+02],
            [0.27423E+02,0.12925E+02],
            [0.25537E+02,0.13517E+02],
            [0.24185E+02,0.14965E+02],
            [0.23469E+02,0.16810E+02],
            [0.22020E+02,0.18162E+02],
            [0.20859E+02,0.19763E+02],
            [0.20506E+02,0.21714E+02],
            [0.20157E+02,0.23642E+02],
            [0.18208E+02,0.24003E+02],
            [0.16372E+02,0.24739E+02],
            [0.14529E+02,0.24008E+02],
        ]
    )
    path_amber19_OPC = np.array(
        [
            [0.34490E+02,0.92100E+01],
            [0.32441E+02,0.97030E+01],
            [0.30536E+02,0.10668E+02],
            [0.31938E+02,0.12173E+02],
            [0.31585E+02,0.14251E+02],
            [0.30464E+02,0.15931E+02],
            [0.28657E+02,0.17055E+02],
            [0.26567E+02,0.17434E+02],
            [0.24577E+02,0.18186E+02],
            [0.22665E+02,0.19133E+02],
            [0.20542E+02,0.19219E+02],
            [0.18585E+02,0.20002E+02],
            [0.17280E+02,0.21652E+02],
            [0.17133E+02,0.23777E+02],
            [0.15782E+02,0.25395E+02],
            [0.14210E+02,0.26840E+02]
        ]
    )
    pair_indices = np.array([[119,24],[130,64]])
    offset = 11
    T = 300
    k = 100
    lag = 200

    feat = create_feat(
        pdb=pdb,
        pair_indices=pair_indices,
        offset=offset
        )
    
    data = load_data(feat,files)
    cluster = clustering(data,k)
    msm = create_msm(cluster,lag)

    plot_reweighted_free_energy(data,msm,T,path=path_amber19_OPC*1e-1)