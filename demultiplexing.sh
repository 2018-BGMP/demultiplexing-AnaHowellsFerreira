
#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=demultiplexing.py
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14


module load prl
module load python/3.6.0
which python

#location of my python script on Talapas
/projects/bgmp/ahowells/demultiplexing/demultiplexing_part2/demultiplexing.py
