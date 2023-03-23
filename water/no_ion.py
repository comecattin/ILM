#! /usr/bin/env python3
"""
Extract the volume of the water + ions and get the protein volume
"""
import numpy as np
import extract_gmx_energy as xtract
import os
import protein_volume

class no_ion:

    def __init__(self,path_water,path_protein):
        
        self.path_water = path_water
        self.path_protein = path_protein
        
        self.state_dict = {}

    def mean_volume(self):
        state_dict = {}
        for name in os.listdir(self.path_water):
            state = name[5:9]
            volume_file = os.path.join(
                self.path_water,name,'volume_md_11.xvg'
            )
            error_file = os.path.join(
                self.path_water,name,'errest_md_11.xvg'
            )
            ouput = xtract.gromacs_output(
                file_name=volume_file,
                error_file=error_file
            )
            mean_error = ouput.analyze()
            state_dict[state] = list(mean_error)
        self.state_dict = state_dict
        
    def get_total_volume_simulation(self,concatenate=False):
        for name in os.listdir(path_protein):
            if name[5:9] in self.state_dict:
                config = protein_volume.configuration(int(name[7:9]),name[5:7])
                if concatenate:
                    config.concatenate(self.path_protein)
                concatenated = os.path.join(
                    self.path_protein,name,'md_concatenated.txt'
                    )
                t_Vprot_mean_error = config.load_txt(concatenated)
                self.state_dict[name[5:9]].append(t_Vprot_mean_error)

    def no_water(self):
        self.waterless = {}
        for state in self.state_dict:

            #in nm^3
            water_volume = self.state_dict[state][0]
            water_error = self.state_dict[state][1]

            time = self.state_dict[state][-1][0]
            protein_volume = self.state_dict[state][-1][1]
            
            waterless = protein_volume - water_volume

            self.waterless[state] = [time,waterless]
        












if __name__ == '__main__':
    path_water = '/home/ccattin/Documents/EAU/NO_COUNTERIONS'
    path_protein = '/home/ccattin/Documents/EAU/HSP90_simulation'
    NoIon = no_ion(path_water,path_protein)
    NoIon.mean_volume()
    NoIon.get_total_volume_simulation()
    NoIon.no_water()