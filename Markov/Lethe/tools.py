#! /usr/bin/env python3
"""
Tools used in LETHE
"""

import MDAnalysis as mda


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

if __name__ == '__main__':

    refGS = '/data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS_cluster1.pdb'
    pairNames = ['64_CA-130_CA', '119_CA-24_CA']
    print(create_pairIndices_from_pairNames(pdbfilename = refGS, pairNames = pairNames ))