#!/bin/bash
#MSUB -A p30190
#!/bin/bash
#SBATCH -A p30190             ## account
#SBATCH -p gengpu             ## "-p" instead of "-q"
#SBATCH -N 1                 ## number of nodes
#SBATCH --mem 40G
#SBATCH -t 2:00:00          ## walltime
#SBATCH --gres=gpu:k80:1
#SBATCH	--job-name="py_ivf"    ## name of job

module purge all
module load anaconda3           ## Load modules (unchanged)

cd /projects/p30190/py_mar
source activate my-numba

~/.conda/envs/my-numba/bin/python -u main.py
exit
