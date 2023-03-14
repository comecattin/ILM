#! /usr/bin/env python3
"""Extract and write distance"""
import matplotlib.pyplot as plt
import MDAnalysis as mda
import scipy as sp
import numpy as np
import pyemma
import os

#%%
def get_dir(DATADIR,OUTDIR):
    """Get and store in variables the different path
    
    Parameters
    ----------
    DATADIR : str
            Path to the data directory
    OUTDIR : str
            Path to the output directory
    Returns
    -------
    pdb : str
        Path to the .pdb file
    refGS : str
        Path to the .pdb reference ground state
    refES : str
        Path to the .pdb reference excited state
    allxtc : list
        All the .xtc files 
    """
    
    #If the output directory doesn't exist -> create one
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)

    #PDB file
    pdb  = DATADIR + "/GS_cluster1.pdb"
    #XTC files
    allfiles =  os.listdir(DATADIR)
    GSfiles = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("GS") and fn.endswith(".xtc"))])
    ESfiles = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("ES") and fn.endswith(".xtc"))])
    TRJ15files = sorted([ f"{DATADIR}/{fn}" for fn in allfiles if (fn.startswith("trj15") and fn.endswith(".xtc"))])
    allxtc = GSfiles + ESfiles + TRJ15files

    #Reference pdb files
    refGS = DATADIR+ '/GS_cluster1.pdb'
    refES = DATADIR+ '/ES_cluster1.pdb'

    return pdb, refGS, refES, allxtc

def create_pairIndices_from_pairNames(pdbfilename, pairNames):
    """Get the indices using pair names

    Parameters
    ----------
    pdbfilename : str
        Path to the pdb file
    pairNames : list of str
        Pair names to analyze

    Returns
    -------
    pairsListIndices : list
        Indices of the pair, can be opened with pyEmma
    """

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
    """Get the features

    Parameters
    ----------
    pdb : str
        Path to the pdb file to load
    pairNames : list of str
        Pair names to analyze
    refGS : str
        Path to the .pdb reference ground state

    Returns
    -------
    feat : PyEmma feature object
        Features
    """
    #Initialyze
    feat = pyemma.coordinates.featurizer(pdb)
    #Get the indices
    pair_indices = create_pairIndices_from_pairNames(pdbfilename = refGS, pairNames = pairNames)
    #Add the feature
    feat.add_distances(indices = pair_indices, periodic=True, indices2=None)
    
    return feat

def load_allxtc(allxtc,feat):
    """Load all the .xtc file. Can be long

    Parameters
    ----------
    allxtc : list of str
        contain all the path to the .xtc files
    feat : PyEmma feature object
        Features to analyze

    Returns
    -------
    data : list of np.array
        List containing all the data extracted
    """
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
    """Plot the density and the associated free energy

    Parameters
    ----------
    data_concatenated : np.array
        All the data of all the trajectories
    allxtc : list of str
        contain all the path to the .xtc files
    timestep : float
        Time step between every file save in the MD simulation
    pairNames : list of str
        Pair names to analyze
    OUTDIR : str
            Path to the output directory
    kT : float
        Value of kT
    """
    #Dentisty map
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

    #Free energy map
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
    """Write for every time step and every trajectories the value of the distance

    Parameters
    ----------
    allxtc : list of str
        contain all the path to the .xtc files
    outdir : str
        Path where to save the file
    data : list of np.array
        List containing all the data extracted
    """
    #Explore every trajectory
    for i, name in enumerate(allxtc):
        #Ground state
        if "/GS" in name:
            state = name.split("/")[-1][:4]
            conf = int(state[-2:])
            output_name = os.path.join(
                OUTDIR,
                f"GS_{conf}_distance.txt"
            )
        #Exited state 
        if "/ES" in name:
            state = name.split("/")[-1][:4]
            conf = int(state[-2:])
            output_name = os.path.join(
                OUTDIR,
                f"ES_{conf}_distance.txt"
            )
        else:
            continue
        #Save
        np.savetxt(output_name,
                   data[i])

if __name__ == '__main__' :
    #Temperature
    T = 300
    kT = sp.constants.R *T/1000

    #Input
    ClusterType = "kmeans"
    DataType = "distances_64CA-130CA_119CA-24CA"
    LagSteps = 2000
    ClusterNumber = 500
    timestep = 0.1
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    #Path definition
    DATA = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC'
    OUTDIR = '/home/ccattin/Documents/Markov/volume_pressure/Output_distances_64CA-130CA_119CA-24CA/'
    
    #%%
    #Get the directories, 
    #       the indices
    #       and the features
    pdb, refGS, refES, allxtc = get_dir(DATA,OUTDIR)
    indices = create_pairIndices_from_pairNames(pdbfilename = refGS, 
                                                pairNames = pairNames )
    feat = pyemma_feat(pdb,pairNames,refGS)
    
    #%%
    #Load the date
    #Warning : this can be really long
    data = load_allxtc(allxtc,feat)
    
    #%%
    #Plot the histogramm of the features
    data_concatenated = np.concatenate(data)
    fig, ax = pyemma.plots.plot_feature_histograms(
        data_concatenated,
        feature_labels=feat,
        ignore_dim_warning=True
        )
    # %%
    #Plot the density and the free energy landscape associated
    plot_density_free_energy(
        data_concatenated,
        allxtc,
        timestep,
        pairNames,
        OUTDIR,
        kT
    )
    #%%
    #Write the results
    write_distance(allxtc,
                   OUTDIR,
                   data)
# %%
