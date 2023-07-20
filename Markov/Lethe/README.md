# LETHE
<img width="472" alt="LETHE_logo" src="https://github.com/comecattin/ILM/assets/75748278/dba00363-eadb-4370-8411-e324f2e8e245">

Lethe is a Client line interface to automatize MSM building and analysis on PyEmma.

## Instalation
To install LETHE:
  - Clone the GitHub repo
  - Install the latest version of PyEmma: https://github.com/markovmodel/PyEMMA

## How to run
Lethe is called by the command `./LETHE.py [options]`.
Since the number of options can be large, it is recommanded to use a `.sh` file to submit the job.
An example of submitting `.sh` file is given

### Input parameters and files
  - `-h` Show the help message and exit
  - `-f [FILES]` Input the `.xtc` GROMACS trajectory files
  - `--ram` Load all the file in the RAM. More efficient but use a lot of RAM
  - `--no-plot` Do not display the plot every time. Usefull when using LETHE on HPC clusters
  - `--T TEMPERATURE` `--temperature TEMPERATURE` Temperature of the system.
  - `-o OUTDIR` `--outdir OUTDIR` Output directory to save all the produced data and plots
  - `--stride STRIDE` Number of stride to consider. Data will be read at regular stride interval
  - `--save FILENAME MODELNAME` File name and model name to save the clustering and the MSM built
  - `--load FILENAME MODELNAME` Load previous model. File name (`.pyemma`) and model name
  - `-p [PLOTS]` Plot the different arguments. Please refer to the next section.

### Features selection
  - `--vamp_feat_type {torsion, distance, txt}` Test the VAMP2 score on different feat type.
  - `--lags_vamp_feat LAGS` The different lag time to perform the VAMP2 score as a function of the feat type.
  - `-d DISTANCES` `--distances DISTANCE` List of the pair names to consider their CA distances.
  - `--residue [RESIDUE]` Use the shortest distance between the residue indices given.
    Pair of residue are separated with a space and '-' separate the two residue inside a pair.
  - `--feat-txt FILE MAXQUALITY` Load the feature from a `.txt` file with a maximum quality to consider.
  - `--vamp-score` Compute the VAMP2 score of the selected features.
  - `--axis AXISX AXISY` Axis to project the plots. Axis begin at 0. By default the two first dimensions (0, 1)

### Dimension reduction
  - `--vamp-lags [LAGS]` List of the different dimension reduction lags to test.
  - `--vamp-dim DIM` Maximum number of dimension to test the dimension reduction.
  - `--reduction {pca, tica, vamp, none}` Which dimension reduction to perform. PCA, TICA or VAMP dimension reduction of no dimension reduction.
  - `--tica-lag LAG` Lag time for the TICA dimension reduction.
  - `--vamp-lag LAG` Lag time for the VAMP dimension reduction.
  - `--dim DIM` Number of dimension to project the dimension reduction.

### Clustering
  - `--cluster {kmeans, regspace}` Clustering method to use. 
  - `--vamp-cluster [CLUSTER_CENTERS]` Compute the VAMP2 score as a function of the number of cluster center.
    This give the optimal number of cluster for after.
  - `-k CLUSTERS` `--cluster-number CLUSTERS` Number of micro-clusters to consider

### MSM lag
  - `--its [ITS]` List of the lag times to perform the ITS convergence validation.
  - `--nits N` Number of eigenvalues to consider in the ITS convergence validation.
  - `--its-cluster` ITS convergence as a function of the number of cluster. Not really usefull...
  - `--lag LAG` Lag time to use for the MSM.

### PCCA++ Analysis
  - `--pcca` perform a pcca++ and TPT analysis to cluster the phase space into macro-clusters. Display the stationnary probabilities.
  - `--state STATE` Number of macro-clusters to consider.
  - `--pdb CONFORMATIONS` Extract the structures of the macro-clusters. Generate `.pdb` files. The number of conformation to consider is given in the option.
  - `--state-path START END` Get the TPT between the two given states.

### MSM Validation
  - `--cktest {True, False}` Perform a CK Test on the given number of macro-state. Option `True` for having the error, `False` for not.
  - `--confidence` Create a Bayesian MSM to get the error. Recommended. Compute the 95% confidence interval.

## Possible plots
LETHE output plots. The following options are to put after the `-p` flag.

### Features selection
- `feat_hist` Plot the feature histogram of the selected features.
![pca_histogram](https://github.com/comecattin/ILM/assets/75748278/17b422b0-c511-48d2-b082-7e89e834cc25)

- `vamp_feat_type` Plot the VAMP2 score as a function of the feature type
![vamp_compare_feat](https://github.com/comecattin/ILM/assets/75748278/8e270df3-969c-4bba-879d-66ebbb72a38d)

- `density_energy` Extract the free energy landscape directly from the sample density.
  ![pca_free_energy_direct_from_density](https://github.com/comecattin/ILM/assets/75748278/5ac73581-7abb-4978-b0ca-e6f7740b318b)

### Dimension reduction

- `vamp_lag_dim` Plot the VAMP2 score as a function of the number of considered dimension and the TICA lag time.
  ![vamp_lag_dim](https://github.com/comecattin/ILM/assets/75748278/f4f021d7-3aba-45c7-950f-2e4925ebc8de)

- `pca` ; `tica` ; `vamp` Plot the histogram of the reduced dimensions.
  ![tica_histogram](https://github.com/comecattin/ILM/assets/75748278/22aedb93-b4ce-40ae-87e8-6eaa888594b0)

### Clustering

- `cluster` Plot the micro-clustered phase space.
  ![cluster](https://github.com/comecattin/ILM/assets/75748278/938b5e6e-14e0-468d-9b94-25d00084382b)

### MSM lag

- `its` Plot the ITS convergence test.
  ![its_validation](https://github.com/comecattin/ILM/assets/75748278/306c3ef1-d28d-43de-9264-625664bfb427)

### PCCA++ analysis

- `stationary` Plot the stationary distribution
  ![stationary_distribution](https://github.com/comecattin/ILM/assets/75748278/a0e37c84-f7eb-43e3-8a5a-f337bbd68a80)

- `eigenvalues` ; `eigenvectors` Plot the spectrum analysis of the MSM.
  ![eigen_values](https://github.com/comecattin/ILM/assets/75748278/f58fa591-6a51-44b0-bf82-3e7e8dbfe4e1)
  ![eigen_vect](https://github.com/comecattin/ILM/assets/75748278/44c0769b-7d03-46be-80b1-7b00d85b8509)

- `reweight_free_energy` Plot the free energy landscape reweighted by the stationary distribution.
  ![reweighted_free_energy](https://github.com/comecattin/ILM/assets/75748278/2a442b03-2365-46f1-903d-5cdeec9a3433)

- `metastable_membership` Extract the metastable membership used by the PCCA++ algorithm.
  ![metastable_membership](https://github.com/comecattin/ILM/assets/75748278/32b6adcd-9a4c-4925-8235-68fddf5100b4)

- `mfpt` Plot the macro-state map with all the associated mean first passage time.
  ![mfpt](https://github.com/comecattin/ILM/assets/75748278/6833f383-b51b-4b69-aaa7-feee213d3559)

- `state_map` Plot the macro-state map. Without the MFPT.
  ![state_map](https://github.com/comecattin/ILM/assets/75748278/5e2690fe-3c2a-4144-ad0f-b45a54a98d0d)

- `committor` Plot the commitor between the chosen two macro-states.
  ![mfpt_committor](https://github.com/comecattin/ILM/assets/75748278/4b7a41db-d832-4f3d-93ec-43ea70b9bdb3)

  ### MSM validation

  - `cktest` Plot the CK test on the chosen number of macro-states.
    ![cktest](https://github.com/comecattin/ILM/assets/75748278/cc3da527-49ad-45b2-9aad-c4362d152046)

