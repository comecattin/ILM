#! /bin/bash

# Bash script to automatize the alignment of protein

pdb_path=$1
vmd_exe=vmd
alignment_script=/home/ccattin/Documents/Code/VMD/alignment_script.tcl
# Launch VMD
$vmd_exe -f $1 -e $alignment_script