#! /bin/bash
# Executable
./LETHE.py \
`# Files to load` \
-f /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/GS0{1..4}_md_all_fitBB_protonly.xtc \
`# Topology .pdb file` \
-t /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/ES_cluster1.pdb \
`# Distances to analyse` \
-d 64_CA-130_CA 119_CA-24_CA 65_CA-130_CA 66_CA-130_CA \
--residue 16-109 17-109 18-109 18-110 19-109 19-110 20-109 20-110 21-109 21-110 109-113 109-169 109-170 109-171 109-172 110-114 \
--feat-txt /home/ccattin/Documents/Markov/HSP90/Amber19_OPC_300K/elisa_feat/2022_10_26_Liste_interactions_simulations.txt 1 \
`# Compute the VAMP2 score` \
--vamp-score \
`# Load into RAM` \
--ram \
`# Temperature of the system` \
--T 300 \
`# Plot to draw` \
-p feat_hist density_energy vamp_lag_dim tica vamp its cluster cktest stationary eigenvalues eigenvectors reweight_free_energy metastable_membership mfpt committor \
--axis 0 2 \
`# Do not display the plots` \
`#--no-plot` \
--vamp-lags 1 2 5 10 20 50 100 200 500 1000 \
--vamp-dim 5 \
`# Do a reduction (tica pca vamp or none)` \
--reduction vamp \
`#--tica-lag 1000` \
--vamp-lag 1000 \
--dim 3 \
`# Number of stride` \
--stride 10 \
`# Algorithm for clustering` \
--cluster kmeans \
`# ITS validation on lagtime list` \
--its 20 50 100\
`# Number of iteeration for the ITS validation` \
--nits 4 \
`# ITS validation as a function of the number of clusters` \
--its-cluster 200 400 \
`# Get the number of cluster` \
--vamp-cluster 50 100 200 400 \
`# Number of clusters` \
-k 200 \
`# Lag time` \
--lag 100 \
`# Bayesian MSM` \
--confidence \
`# Number of metastable states to consider` \
--state 3 \
`# PCCA analysis` \
--pcca \
`# Transition Path Theory analysis between two state (begin at 0)` \
--state-path 1 2 \
`# Writte .pdb file of the meta stable states` \
--pdb 10 \
`# Output directories for plot and save` \
-o /home/ccattin/Documents/Code/outputs/LETHE \
`# Save model` \
`#--save all_lag500_k200_stride10_reductionNone.pyemma GSES` \
`# Load previous model` \
`#--load all_lag1000_k200_stride4_reductionNone.pyemma GSES` \
#> lethe.out \
