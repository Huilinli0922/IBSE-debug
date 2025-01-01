#!/bin/bash

#SBATCH --job-name=Matlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --time=100:00:00

module purge
module load matlab/2020b

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)
export MATLAB_LOG_DIR=$SLURM_JOBTMP

cat<<EOF | matlab -nodisplay

pll = gcp('nocreate');
%%parpool('local', $SLURM_CPUS_PER_TASK);
if isempty(pll)
    parobj = parpool('local',$SLURM_CPUS_PER_TASK);
end

main_ice_boat_moving

exit

EOF





