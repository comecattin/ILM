#! /bin/bash
#SBATCH --job-name=HSP90_4
#SBATCH -o ./%x.%j.%N.out           # output file
#SBATCH -e ./%x.%j.%N.err           # errors file
#
#SBATCH -p Cascade
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=1-16:00:00           # day-hours:minutes:seconds
#SBATCH --mem=64G
#
# Variables
code_dir=/home/ccattin/Stage/Code
gmx_code_dir=$code_dir/GMX
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
GMX_EXE="gmx"

run_dir=$(pwd)
step_3_dir=$run_dir/step_3
step_4_dir=$run_dir/step_4
input_dir=$step_4_dir/input
input_global_dir=$run_dir/input
sub_script_dir=$run_dir/script
log_dir=$run_dir/LOGS

printf "
4th Production step : 100 ns
"
# Just for the copy
shopt -s extglob

# Heading
echo "The job ${SLURM_JOB_ID} is running on these nodes:"
echo ${SLURM_NODELIST}
echo
#
cd $step_4_dir

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
mkdir -p $SCRATCH/production
OUTPUTDIR=$SCRATCH/OUTPUTDIR
mkdir -p $OUTPUTDIR
cd $SCRATCH/production
#Copy files
cp $input_dir/prod.mdp .
ln -s $step_3_dir/production/prod.gro prev.gro
ln -s $step_3_dir/production/prod.cpt prev.cpt
ln -s $input_dir/topol.top .
ln -s $input_dir/index.ndx .
#Running gmx
echo "        Running gmx grompp"
$GMX_EXE grompp -f prod.mdp -c prev.gro -t prev.cpt -p topol.top -o prod.tpr -n index.ndx &> $OUTPUTDIR/grompp_prod.out
echo "    Running gmx mdrun"
$GMX_EXE mdrun -deffnm prod -nt $SLURM_CPUS_PER_TASK &> $OUTPUTDIR/mdrun.out
echo "        Production run done."

cp -r $SCRATCH/* $step_4_dir
cd $run_dir
rm -r $SCRATCH

step_5_dir=$run_dir/step_5
mkdir -p $step_5_dir
cp -r $input_global_dir/step_1 $step_5_dir/input
mv *.err $log_dir
mv *.out $log_dir
sbatch $sub_script_dir/submit_5