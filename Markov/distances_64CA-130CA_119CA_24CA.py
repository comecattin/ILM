#! /usr/bin/env python3
#%%
import matplotlib.pyplot as plt
import MDAnalysis as mda
import scipy as sp
import numpy as np
import pyemma
import os
# for visualization of molecular structures:

# Importing the library to get the memory usages 
#%%
def get_dir():
    WD = os.getcwd()
    PICTDIR = WD + '/Figures'
    #DATADIR = WD + '/Data/'
    DATADIR = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC'
    OUTDIR =  WD + f'/Output_{DataType}_{ClusterType}{ClusterNumber}_lag{LagSteps}/'

    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    pdb  = DATADIR + "/GS_cluster1.pdb"
    allfiles =  os.listdir(DATADIR)
    GSfiles = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("GS") and fn.endswith(".xtc"))])
    ESfiles = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("ES") and fn.endswith(".xtc"))])
    TRJ15files = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("trj15") and fn.endswith(".xtc"))])

    allxtc = GSfiles + ESfiles + TRJ15files
    #allxtc = TRJ15files

    refGS = DATADIR+ '/GS_cluster1.pdb'
    refES = DATADIR+ '/ES_cluster1.pdb'

    return pdb, refGS, refES, allxtc, OUTDIR

def create_pairIndices_from_pairNames(pdbfilename, pairNames):
    refu = mda.Universe(pdbfilename)
    pairsListIndices = []
    for pairname in pairNames:
        print(f'search indices for the pair {pairname}')
        atomName1 =  pairname.split('-')[0]
        res1 = atomName1.split('_')[0]
        name1 = atomName1.split('_')[1]

        atomName2 =  pairname.split('-')[1]
        res2 = atomName2.split('_')[0]
        name2 = atomName2.split('_')[1]
        
        print(name1, res1)
        index1 =   refu.select_atoms(f"name {name1} and resid {res1}").indices
        print(name2, res2)
        index2 =   refu.select_atoms(f"name {name2} and resid {res2}").indices
        
        if len(index1) == 1 and len(index2) == 1 :
            pairsListIndices.append([index1[0],index2[0]])
        else :
            print(f'WARNING : the pair defined by the name {pairname} do not lead to a pair on indices')
            
    return pairsListIndices

def pyemma_feat(pdb,pairNames,refGS):
    feat = pyemma.coordinates.featurizer(pdb)
    pair_indices = create_pairIndices_from_pairNames(pdbfilename = refGS, pairNames = pairNames)

    feat.add_distances(indices = pair_indices, periodic=True, indices2=None)
    return feat

def load_allxtc(allxtc,feat):
    data = pyemma.coordinates.load(allxtc,
                                   features=feat,stride=1
                                   )
    return data

def plot_density_free_energy(data_concatenated,
                             allxtc,
                             timestep,
                             pairNames,
                             OUTDIR,
                             kT):

    fig, axes = plt.subplots(1, 
                             1, 
                             figsize=(6, 4), 
                             sharex=True, 
                             sharey=True)
    
    pyemma.plots.plot_density(
        *data_concatenated.T[0:2],
        ax=axes
        )
    
    plt.title(
        f'All {len(allxtc)}  simulations, stored every {timestep} ns',
        fontweight='bold')
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')
    fig.tight_layout()
    
    fig.savefig(f'{OUTDIR}/data_density.png')

    
    fig, axes = plt.subplots(
        1,
        1,
        figsize=(6, 4),
        sharex=True,
        sharey=True
        )
    fig,axes = pyemma.plots.plot_free_energy(
        *data_concatenated.T[0:2],
        ax=axes,
        kT=kT,cbar_label='free energy / kJ.mol-1'
        )  

    plt.title(
        f'All {len(allxtc)} simulations, stored every {timestep} ns',
        fontweight='bold')
    axes.set_xlabel(f' {pairNames[0]} (nm)')
    axes.set_ylabel(f' {pairNames[1]} (nm)')
    fig.tight_layout()
    fig.savefig(f'{OUTDIR}/data_free_energy_direct_from_density.png')

def write_distance(allxtc,outdir,data):
    for i, name in enumerate(allxtc):
        
        if "/GS" in name:
            state = name.split("/")[-1][:4]
            conf = int(state[-2:])
            output_name = os.path.join(
                OUTDIR,
                f"GS_{conf}_distance.txt"
            )
            
        if "/ES" in name:
            state = name.split("/")[-1][:4]
            conf = int(state[-2:])
            output_name = os.path.join(
                OUTDIR,
                f"ES_{conf}_distance.txt"
            )
        else:
            continue
        
        np.savetxt(output_name,
                   data[i])

write_distance(allxtc,OUTDIR,data)
    #%%

if __name__ == '__main__' :
    

    T = 300
    kT = sp.constants.R *T/1000

    ClusterType = "kmeans"
    DataType = "distances_64CA-130CA_119CA-24CA"
    LagSteps = 2000
    ClusterNumber = 500
    timestep = 0.1
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    
    #%%
    pdb, refGS, refES, allxtc, OUTDIR = get_dir()
    indices = create_pairIndices_from_pairNames(pdbfilename = refGS, 
                                                pairNames = pairNames )
    feat = pyemma_feat(pdb,pairNames,refGS)
    
    #%%
    data = load_allxtc(allxtc,feat)
    
    #%%
    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(data_concatenated, feature_labels=feat,ignore_dim_warning=True)
    # %%
    plot_density_free_energy(
        data_concatenated,
        allxtc,
        timestep,
        pairNames,
        OUTDIR,
        kT
    )

# %%
