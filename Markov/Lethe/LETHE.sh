#! /bin/bash
./LETHE.py \
-f /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/ES0{1..2}_md_all_fitBB_protonly.xtc \
-t /data/cloison/Simulations/HSP90-NT/SIMULATIONS_TRAJECTORIES/AMBER19SB_OPC/ES_cluster1.pdb \
-d 64_CA-130_CA 119_CA-24_CA \
--T 300 \
-p feat_hist density_energy its \
--pca \
--tica \
--lag 1000 \
--cluster kmeans \
-k 100 \
--stride 1 \
--its 1 2 5 10 20 50 \
--nits 4 \
-o /home/ccattin/Documents/Code/outputs