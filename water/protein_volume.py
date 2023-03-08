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
    
    def __init__(self,path,water_file):
        self.path=path
        self.water_file=[]
        self.water_file.append(os.path.join(path,water_file[0]))
        self.water_file.append(os.path.join(path,water_file[1]))

    def extract_volume_command(self,code):
        subprocess.call([code,self.path])
    
    def number_water(self):
        water=[[],[]]
        for i in (0,1):
            with open(self.water_file[i],'r') as f:
                for line in f:
                    water[i].append(int(line.split()[-1]))
        return water


            

if __name__ == "__main__":
    path="/home/ccattin/Documents/EAU/HSP90_simulation"
    water_file=["ES.txt","GS.txt"]
    HSP90 = protein(
        path=path,
        water_file=water_file
        )
    code="/home/ccattin/Documents/Code/GMX/analysis_water_protein"
    #HSP90.extract_volume_command(code=code)
    ES = HSP90.number_water()