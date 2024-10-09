#!/bin/bash -l
#
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=64
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=deg270C
#SBATCH --partition=standard
#SBATCH --time=0-01:00:00
#SBATCH --output=mylog.out
#SBATCH --error=myfail.out
#SBATCH --mail-user='fyshi@udel.edu'
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
vpkg_require openmpi
#UD_QUIET_JOB_SETUP=YES
#UD_USE_SRUN_LAUNCHER=YES
#UD_DISABLE_CPU_AFFINITY=YES
#UD_MPI_RANK_DISTRIB_BY=CORE
#UD_DISABLE_IB_INTERFACES=YES
#UD_SHOW_MPI_DEBUGGING=YES

. /opt/shared/slurm/templates/libexec/openmpi.sh

${UD_MPIRUN} /home/1047/Fully_dispersive_model/src_FD/f2w_vessel
mpi_rc=$?


#
exit $mpi_rc
