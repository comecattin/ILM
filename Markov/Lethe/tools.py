#! /usr/bin/env python3
"""
Tools used in LETHE
"""

import MDAnalysis as mda
import scipy
import pyemma
import numpy as np


def create_pairIndices_from_pairNames(pdbfilename, pairNames):
    """Get the indices from the names

    Parameters
    ----------
    pdbfilename : str
        Path to the .pdb file
    pairNames : list
        List of string containing the different pair name

    Returns
    -------
    pairsListIndices : list
        List containing the different indices
    """
    refu = mda.Universe(pdbfilename)
    pairsListIndices = []
    for pairname in pairNames:
        print(f"Search indices for the pair {pairname}")
        atomName1 = pairname.split("-")[0]
        res1 = atomName1.split("_")[0]
        name1 = atomName1.split("_")[1]

        atomName2 = pairname.split("-")[1]
        res2 = atomName2.split("_")[0]
        name2 = atomName2.split("_")[1]

        # print(name1, res1)
        index1 = refu.select_atoms(f"name {name1} and resid {res1}").indices
        # print(name2, res2)
        index2 = refu.select_atoms(f"name {name2} and resid {res2}").indices

        if len(index1) == 1 and len(index2) == 1:
            pairsListIndices.append([index1[0], index2[0]])

        else:
            print(
                f"WARNING : the pair defined by the name {pairname} do not lead to a pair on indices"
            )

    print(f"Found indices: {pairsListIndices}")

    return np.array(pairsListIndices)

def create_pairIndices_from_indices(pairNames):
    """Get the indices list from the indices input

    Parameters
    ----------
    pairNames : list
        List containing the indices. Pairs are separated by a space. Atoms inside a pair are separated by '-'

    Returns
    -------
    pairsListIndices : list
        List containing the different indices
    """
    # Initialization
    pairsListIndices = []
    # Loop over all the pair and convert
    for i, pair in enumerate(pairNames):
        pair = pair.split('-')
        pair = [int(i) for i in pair]
        pairsListIndices.append(pair)

    print(f"Found indices: {pairsListIndices}")
    
    return np.array(pairsListIndices)


def get_kT(T):
    """Compute kT

    Parameters
    ----------
    T : float
        Temperature of the system

    Returns
    -------
    kT : float
        kT
    """
    return scipy.constants.R * T / 1000


def save_model(cluster, msm, outdir, filename, model_name):
    """Save the given model in a .pyemma file

    Parameters
    ----------
    cluster : pyemma.cluster type
        PyEmma cluster
    msm : pyemma.msm
        PyEmma estimation MSM
    outdir : str
        Output directory
    filename : str
        Name of the file in the output directory
    model_name : str
        Name to give to the model
    """
    # Save the cluster
    cluster.save(
        f"{outdir}/{filename}", model_name=f"{model_name}_cluster", overwrite=True
    )
    # Save the MSM
    msm.save(f"{outdir}/{filename}", model_name=f"{model_name}_msm", overwrite=True)
    # Confirmation print
    print(f"Cluster and MSM saved in {outdir}/{filename} with model name {model_name}")


def load_model(outdir, filename, model_name):
    """Load previous PyEmma model

    Parameters
    ----------
    outdir : str
        Output directory
    filename : str
        Name of the file in the output directory
    model_name : str
        Name to give to the model

    Returns
    -------
    cluster : pyemma.cluster type
        PyEmma cluster
    msm : pyemma.msm
        PyEmma estimation MSM
    """
    # Load MSM
    msm = pyemma.load(f"{outdir}/{filename}", model_name=f"{model_name}_msm")
    # Load cluster
    cluster = pyemma.load(f"{outdir}/{filename}", model_name=f"{model_name}_cluster")
    # Confirmation print
    print(
        f"Cluster and MSM loaded from {outdir}/{filename} with model name {model_name}"
    )

    return msm, cluster

def read_feat_from_txt(file_path,quality_max):
    """Read features from .txt files

    Parameters
    ----------
    file_path : str
        Path to the .txt file. This file contain all the interaction we want to consider. Each line is an interaction. The first line is a header.
        - Column 1 : index of the interaction (not useful)
        - Column 2 : index of the i residue
        - Column 3 : index of the j residue
        - Column 4 : quality of the interaction (1: good quality then decreasing, must be int)
    
    quality_max : int
        Maximum quality to consider (1: good quality then decreasing, must be int)

    Returns
    -------
    selection : np.array
        Array containing the selected atom indices
    """

    data = np.loadtxt(
        file_path,
        skiprows=1,
        usecols=(0,1,2,3),
        dtype=int
    )

    selection = data[data[:,-1] == quality_max]
    pair_indices = selection[:,1:3]
    unique = remove_double(pair_indices)
    print(f"Found indices {unique}")
    
    return unique

def remove_double(pair_indices):
    sorted = np.sort(pair_indices,axis=1)
    unique = np.unique(sorted,axis=0)
    return unique



if __name__ == "__main__":
    refGS = "/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb"
    pairNames = ["64_CA-130_CA", "119_CA-24_CA"]
    print(create_pairIndices_from_pairNames(pdbfilename=refGS, pairNames=pairNames))
    pairNames = ['16-109','17-109','18-109']
    print(create_pairIndices_from_indices(pairNames))
    file_path = '/home/ccattin/Documents/Markov/HSP90/Amber19_OPC_300K/elisa_feat/2022_10_26_Liste_interactions_simulations.txt'
    pair_indices = read_feat_from_txt(file_path,1)
    pair_indices = remove_double(pair_indices=pair_indices)

