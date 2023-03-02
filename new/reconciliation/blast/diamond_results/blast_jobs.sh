#!/bin/bash
#SBATCH -p admintest
#SBATCH --mem-per-cpu=10G
#SBATCH -t 8:00:00
#SBATCH --array=0-33
#SBATCH --job-name=./blast_jobs.txt
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.evans@yale.edu
#SBATCH --ntasks=1

/ysm-gpfs/apps/software/dSQ/0.92/dSQBatch.py ./blast_jobs.txt
