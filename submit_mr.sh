#!/bin/bash 
#
module load R/3.4.0
module load plink2
module load zstd
today=`date +%m%d%y`
dirLog="$mr/log/"
[ -d "$dirLog" ] || mkdir -p "$dirLog"

# run multiple `job_mr.R' each with parameter rank[i], prev[i], weighted[i]
# (mmem[i] memory is used for each job)

# log files written to $dirLog/mr_${today}_${rank[$idx]}_${weighted[$idx]}.err, $dirLog/mr_${today}_${rank[$idx]}_${weighted[$idx]}.out

rank=(8 16 35)
prev=(0 0 0) # which iteration to start with; -1 means starting from the last intermediate results
mmem=(250000 250000 250000)
weighted=(0 0 0)
use_plink2=(1 1 1) # whether to use plink2 for matrix multiplication, 1 for 2.0, 0 for 1.9.

for idx in "${!rank[@]}"; do
  sbatch -p stat,mrivas --mail-type=END,FAIL --mail-user=junyangq@stanford.edu \
       -c 16 \
       --time=2:00:00 --mem=${mmem[$idx]} \
       --output=$dirLog/mr_${today}_${rank[$idx]}_${weighted[$idx]}.out\
       --error=$dirLog/mr_${today}_${rank[$idx]}_${weighted[$idx]}.err \
       --job-name="mr_${rank[$idx]}_${weighted[$idx]}" \
       --wrap="Rscript ./job_mr.R ${rank[$idx]} ${prev[$idx]} ${weighted[$idx]} ${use_plink2[$idx]}"
done
