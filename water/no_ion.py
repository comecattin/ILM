#! /usr/bin/env python3
"""
Extract the volume of the water + ions and get the protein volume
"""
import numpy as np
import extract_gmx_energy as xtract
import os

class no_ion:

    def __init__(self,path):
        self.path = path

    def mean_volume(self):
        state_dict = {}
        for name in os.listdir(self.path):
            state = name[5:9]
            volume_file = os.path.join(
                self.path,name,'volume_md_11.xvg'
            )
            error_file = os.path.join(
                self.path,name,'errest_md_11.xvg'
            )
            ouput = xtract.gromacs_output(
                file_name=volume_file,
                error_file=error_file
            )
            mean_error = ouput.analyze()
            state_dict[state] = mean_error
        















if __name__ == '__main__':
    path = '/home/ccattin/Documents/EAU/NO_COUNTERIONS'
    NoIon = no_ion(path)
    NoIon.mean_volume()