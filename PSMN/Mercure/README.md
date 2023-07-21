# Mercure
![mercure](https://github.com/comecattin/ILM/assets/75748278/7d877227-4057-4bf8-b9a2-9c886f34c962)

Mercure is a code to automatize multiple submission scripts on the PSMN clusters

## Installation

To install, please clone this repo.

## How to use

### One state to simulate
  - Into the executable `mercure` change and adapt the job name.
  - Adapt all the `input/step_i` files. `step_0` is the first step, it include solvatation of the system, equilibration and production. `step_1` include only production from previously equilibrated simulations.
  - Inside the `script/` directory adapt all the submission script. Can be done using `sed` I think. In each script, change the `$code_dir` variable by the directory of the cloned repo.
  - Change the `$SCRATCH` variable. This is where the computation is done.
  - Adapt all the names, node, time and memory.
  - Run with `sbatch mercure`.


### Multiple states to simulate
  - Into the working directory, create a directory named `GRO` with all the `.gro` file inside state directory. Ex: `WORKING_DIR/GRO/ES01/processed.gro`.
  - Run with `./multiple_input`. The programm then loop over all the `.gro` file and submit for you the job.
