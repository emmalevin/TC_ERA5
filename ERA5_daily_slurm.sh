#!/bin/bash
#SBATCH --job-name=myjob            # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks=1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=32G            # memory per cpu-core (4G is default)
#SBATCH --time=05:00:00             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin,end,fail  # receive email notifications
#SBATCH --mail-user=el2358@princeton.edu

module purge
module load anaconda3/2025.6
conda activate /home/GEOCLIM/software/anaconda3/2024.10/envs/geoclim

cd /home/el2358/notebooks/2020/
python ERA5_daily.py