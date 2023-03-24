#! /usr/bin/env python3
"""
Extract the volume of the water + ions and get the protein volume
"""
import numpy as np
import matplotlib.pyplot as plt
import extract_gmx_energy as xtract
import os
import protein_volume
import seaborn as sns

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
        
    def plot_volume_state(self,output):
        fig, axs = plt.subplots(len(self.waterless),
                                sharex=True,
                                figsize=(5,7))
        color = sns.color_palette('cool',12)
        for i, state in enumerate(self.waterless):
            x = self.waterless[state][0]
            y = self.waterless[state][1]
            
            config = protein_volume.configuration(None,None)
            smooth = config.smoothing(x,y,100)

            axs[i].plot(x,y,color=color[6])
            axs[i].plot(smooth,label='Smoothed')
            axs[i].text(0,32.5,state)
            axs[i].set_ylim(24,34)
            axs[i].set_xlim(0,100000)
        axs[0].legend()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.text(0.5, 0.04, 'Time (ps)', ha="center", va="center")
        fig.text(0.03, 0.5, r'Volume (nm$^3$)', ha="center", va="center", rotation="vertical")
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.show()  

    def results(self):
        print("State    Water    Error Water    Total    Total error    Protein    Protein error")
        for state in self.state_dict:
            water = self.state_dict[state][0]
            water_error = self.state_dict[state][1]
            total = self.state_dict[state][-1][2]
            total_error = self.state_dict[state][-1][3]
            protein = total - water
            protein_error = total_error + water_error
            print(f"{state}: {water} {water_error} {total} {total_error} {protein} {protein_error}")







if __name__ == '__main__':
    path_water = '/home/ccattin/Documents/EAU/NO_COUNTERIONS'
    path_protein = '/home/ccattin/Documents/EAU/HSP90_simulation'
    output = '/home/ccattin/Documents/Code/outputs/no_ion.pdf'
    NoIon = no_ion(path_water,path_protein)
    NoIon.mean_volume()
    NoIon.get_total_volume_simulation()
    NoIon.no_water()
    #NoIon.plot_volume_state(output=output)
    NoIon.results()