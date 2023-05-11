#! /bin/bash
# Executable
./LETHE.py \
`# Files to load` \
-f /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS0{1..2}_md_all_fitBB_protonly.xtc \
`# Topology .pdb file` \
-t /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/ES_cluster1.pdb \
`# Distances to analyse` \
-d 64_CA-130_CA 119_CA-24_CA \
`# Temperature of the system` \
--T 300 \
`# Plot to draw` \
-p feat_hist density_energy its cluster cktest stationary eigenvectors metastable_membership mfpt committor \
`# Do a reduction (tica pca or none)` \
--reduction none \
`# Lag time` \
--lag 1000 \
`# Algorithm for clustering` \
--cluster kmeans \
`# Number of clusters` \
-k 200 \
`# Number of stride` \
--stride 1 \
`# ITS validation on lagtime list` \
--its 1 2 5 10 20 50 \
`# Number of iteeration for the ITS validation` \
--nits 4 \
`# ITS validation as a function of the number of clusters` \
--its-cluster 20 50 100 \
`# Bayesian MSM` \
--confidence \
`# Number of metastable states to consider` \
--state 3 \
`# PCCA analysis` \
--pcca \
`# Transition Path Theory analysis between two state (begin at 0)` \
--state-path 1 2 \
`# Output directories for plot and save` \
-o /home/ccattin/Documents/Code/outputs \
`# Save model` \
--save test.pyemma GS01_GS02 \
`# Load previous model` \
--load test.pyemma GS01_GS02
