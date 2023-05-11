#! /usr/bin/env python3
"""
Tools used in LETHE
"""

import MDAnalysis as mda
import scipy
import pyemma


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
        print(f'Search indices for the pair {pairname}')
        atomName1 =  pairname.split('-')[0]
        res1 = atomName1.split('_')[0]
        name1 = atomName1.split('_')[1]

        atomName2 =  pairname.split('-')[1]
        res2 = atomName2.split('_')[0]
        name2 = atomName2.split('_')[1]
        
        #print(name1, res1)
        index1 =   refu.select_atoms(f"name {name1} and resid {res1}").indices
        #print(name2, res2)
        index2 =   refu.select_atoms(f"name {name2} and resid {res2}").indices
        
        if len(index1) == 1 and len(index2) == 1 :
            pairsListIndices.append([index1[0],index2[0]])
            
        else :
            print(f'WARNING : the pair defined by the name {pairname} do not lead to a pair on indices')

    print(f'Found indices: {pairsListIndices}')    
    
    return pairsListIndices

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
    return scipy.constants.R *T/1000

def save_model(cluster,msm,outdir,filename,model_name):
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
        f'{outdir}/{filename}',
        model_name=f'{model_name}_cluster',
        overwrite=True
        )
    # Save the MSM
    msm.save(
        f'{outdir}/{filename}',
        model_name=f'{model_name}_msm',
        overwrite=True
    )
    # Confirmation print
    print(
        f'Cluster and MSM saved in {outdir}/{filename} with model name {model_name}'
        )

def load_model(outdir,filename,model_name):
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
    msm = pyemma.load(
        f'{outdir}/{filename}',
        model_name=f'{model_name}_msm'
        )
    # Load cluster
    cluster = pyemma.load(
        f'{outdir}/{filename}',
        model_name=f'{model_name}_cluster'
        )
    # Confirmation print
    print(
        f'Cluster and MSM loaded from {outdir}/{filename} with model name {model_name}'
    )
    
    return msm, cluster

if __name__ == '__main__':

    refGS = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    print(create_pairIndices_from_pairNames(pdbfilename = refGS, pairNames = pairNames ))