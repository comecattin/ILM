#! /bin/bash
# Executable
./LETHE.py \
`# Files to load` \
-f /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/*_md_all_fitBB_protonly.xtc \
`# Topology .pdb file` \
-t /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/ES_cluster1.pdb \
`# Distances to analyse` \
-d 64_CA-130_CA 119_CA-24_CA \
`# Temperature of the system` \
--T 300 \
`# Plot to draw` \
-p feat_hist density_energy tica its cluster cktest stationary eigenvectors metastable_membership mfpt committor \
`# Do not display the plots` \
--no-plot \
`# Do a reduction (tica pca or none)` \
--reduction none \
`#--tica-lag 1000` \
--dim 2 \
`# Number of stride` \
--stride 10 \
`# Algorithm for clustering` \
--cluster kmeans \
`# ITS validation on lagtime list` \
--its 200 500 1000 2000 4000 8000 10000\
`# Number of iteeration for the ITS validation` \
--nits 4 \
`# ITS validation as a function of the number of clusters` \
--its-cluster 200 400 800 1000 2000 4000\
`# Number of clusters` \
-k 200 \
`# Lag time` \
--lag 500 \
`# Bayesian MSM` \
--confidence \
`# Number of metastable states to consider` \
--state 3 \
`# PCCA analysis` \
--pcca \
`# Transition Path Theory analysis between two state (begin at 0)` \
--state-path 1 2 \
`# Output directories for plot and save` \
-o /home/ccattin/Documents/Code/outputs/LETHE \
`# Save model` \
--save all_lag500_k200_stride10_reductionNone.pyemma GSES \
`# Load previous model` \
`#--load all_lag1000_k200_stride4_reductionNone.pyemma GSES` \
> lethe.out \
