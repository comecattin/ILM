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


