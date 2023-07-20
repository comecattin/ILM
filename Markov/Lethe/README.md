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

### Features selection
  - `--vamp_feat_type {torsion, distance, txt}` Test the VAMP2 score on different feat type.
  - `--lags_vamp_feat LAGS` The different lag time to perform the VAMP2 score as a function of the feat type.
  - `-d DISTANCES` `--distances DISTANCE` List of the pair names to consider their CA distances.
  - `--residue [RESIDUE]` Use the shortest distance between the residue indices given.
    Pair of residue are separated with a space and '-' separate the two residue inside a pair.
  - `--feat-txt FILE MAXQUALITY` Load the feature from a `.txt` file with a maximum quality to consider.
  - `--vamp-score` Compute the VAMP2 score of the selected features.

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
