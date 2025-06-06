#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=8mm_free
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=willhies,amarsden
#SBATCH --exclude=sh03-16n02
# #SBATCH --exclude=sh-107-[59-64]

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

TOTAL_TASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))

# executable is here 
SRCDIR=$PWD

# run in scratch, name with the job name
RUNDIR=$SCRATCH/aortic_${SLURM_JOBID/.*}_384_ce2eb5d_0mm_radial_8mm_circ_basic_rcr_extra_free_edge_redone
mkdir $RUNDIR

# set up run info 
BASE_NAME=aortic_384
INPUT_NAME=input_aortic_384
RUN_LINE="srun main3d"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"
SESSION_NAME="aortic_384_visit2pt9.session"
VIEW_CLIPPING="-0.2"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME*                               $RUNDIR
cp aorta_bdry.vertex                         $RUNDIR
cp atrium_bdry.vertex                        $RUNDIR
cp $INPUT_NAME                               $RUNDIR
cp *.cpp                                     $RUNDIR
cp main3d                                    $RUNDIR
cp fourier_coeffs*.txt                       $RUNDIR
cp watchdog_job_restart.py                   $RUNDIR
cp kill_all_mpi.sh                           $RUNDIR
cp run_aortic_N384.sbatch                    $RUNDIR
cp run_parallel_movie.py                     $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS" 

# load stuff for movie making 
source ~/.bash_profile
python run_parallel_movie.py $SESSION_NAME $TOTAL_TASKS $VIEW_CLIPPING

SESSION_NAME_PARAVIEW="~/heart_valves/scripts/velocity_slices_5.py"
visit -cli -nowin -s ~/heart_valves/scripts/run_parallel_convert_visit_to_paraview.py 48 $TOTAL_TASKS
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS
