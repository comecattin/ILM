#! /bin/bash

run_dir=/home/ccattin/Documents/HSP90_278K/multiple_input
cluster_dir=/home/ccattin/Documents/Cluster

cd $run_dir
~/Documents/Code/GMX/concatenate

cd $cluster_dir
for state in ES GS all
do
    dir=$cluster_dir/$state
    mkdir -p $dir
    cd $dir
    ln -s $run_dir/analysis/$state.xtc .
    ln -s $run_dir/ES01/step_10/analysis/prod_prot_fit.gro $state.gro
    ttclust -f $state.xtc -t $state.gro -sr "backbone and (residue 98 to 136)" -sa "backbone and (residue 40 to 97 or residue 137 to 220)" -axis frame
    cd $cluster_dir 
done