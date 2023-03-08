#! /usr/bin/env python3
"""
Extract the molecular volume of protein of multiple simualtion.
"""
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import extract_gmx_energy as xtract
import eau
import os
import subprocess

class protein():
    
    def __init__(self,path,water_file,water_volume):
        self.path=path
        self.water_file=[]
        self.water_file.append(os.path.join(path,water_file[0]))
        self.water_file.append(os.path.join(path,water_file[1]))
        self.water_volume=water_volume

    def extract_volume_command(self,code):
        subprocess.call([code,self.path])
    
    def number_water(self):
        water=[[],[]]
        for i in (0,1):
            with open(self.water_file[i],'r') as f:
                for line in f:
                    water[i].append(int(line.split()[-1]))
        return water
    


    def remove_water(self,water_number):
        
        for name in os.listdir(self.path):
            if ".txt" in name:
                continue
            for traj in range(1,11):
                volume_file=os.path.join(
                    self.path,
                    name,
                    "volume_md_"+str(traj)+".xvg"
                )
                error_file=os.path.join(
                    self.path,
                    name,
                    "errest_md_"+str(traj)+".xvg"
                )
                output = os.path.join(
                    self.path,
                    name,
                    "no_water_md_"+str(traj)+".txt"
                )
                error_file=os.path.join(
                    self.path,
                    name,
                    "errest_md_"+str(traj)+".xvg"
                )
                volume = xtract.gromacs_output(
                    file_name=volume_file,
                    error_file=error_file
                    )
                time, volume_array = volume.extract()
                no_water=volume_array-self.water_volume*1e-3*water_number

                np.savetxt(
                    output,
                    time,
                    no_water
                )

        return no_water

            

if __name__ == "__main__":
    path="/home/ccattin/Documents/EAU/HSP90_simulation"
    water_file=["ES.txt","GS.txt"]
    water_volume=30.72669350142569
    HSP90 = protein(
        path=path,
        water_file=water_file,
        water_volume=water_volume
        )
    code="/home/ccattin/Documents/Code/GMX/analysis_water_protein"
    #HSP90.extract_volume_command(code=code)
    ES, GS = HSP90.number_water()
    no_water = HSP90.remove_water(GS[-1])