#! /bin/bash
#SBATCH --job-name=HSP90_278K
#SBATCH -o ./%x.%j.%N.out           # output file
#SBATCH -e ./%x.%j.%N.err           # errors file
#
#SBATCH -p Cascade
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00           # day-hours:minutes:seconds
#SBATCH --mem=1G
#
shopt -s extglob

nMPI=1
nOMP=96
GMX_EXE="gmx_mpi"
MPRUN_EXE="mpirun"
code_dir=/home/ccattin/Stage/Code

echo "The job ${SLURM_JOB_ID} is running on these nodes:"
echo ${SLURM_NODELIST}
echo
#
cd $SLURM_SUBMIT_DIR    # go to the work / submission directory
#Module
module purge
module use /applis/PSMN/debian11/Cascade/modules/all
module load GROMACS/2021.5-foss-2021b
#Copy to scratch
SCRATCH=/scratch/Cascade/ccattin/${SLURM_JOB_ID}
mkdir -p $SCRATCH
cp -r !(*.err|*.out) $SCRATCH/
cd $SCRATCH/
#Script to run
$code_dir/PSMN/runGMX_PSMN $nMPI $nOMP $GMX_EXE $MPRUN_EXE

cd $SCRATCH/
mv $SCRATCH/* $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
rm -r $SCRATCH
