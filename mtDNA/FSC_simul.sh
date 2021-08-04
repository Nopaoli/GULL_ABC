#!/bin/bash -l

# author: momiglia
#SBATCH --account=project_2000465
#SBATCH -o FSC.ou
#SBATCH -e FSC.err
#SBATCH -p small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paolo.momigliano@helsinki.fi
#SBATCH --array=1-3
LOC=/scratch/project_2000465/BIRD/FSC/mtDNA
MOD=$(sed -n ${SLURM_ARRAY_TASK_ID}p $LOC/models)
/projappl/project_2000465/fsc26_linux64/fsc26 -t $LOC/$MOD/${MOD}.tpl -n1 -e $LOC/$MOD/${MOD}.est -E 10000 
cd  $LOC/$MOD 
sh LaunchArlSumStat.sh
sed -i "s,\tMaxEstLhood\tMaxObsLhood,,g" ${MOD}.params
paste $LOC/$MOD/outSumStats.txt $LOC/$MOD/${MOD}.params > $LOC/SIM_STATS/${MOD}_reftable.txt
for i in IM SI RI; do rm $LOC/$MOD/${MOD}*ar* -rf ; done


