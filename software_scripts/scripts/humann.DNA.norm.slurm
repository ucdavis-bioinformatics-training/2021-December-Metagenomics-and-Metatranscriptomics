#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-12
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/hmadn_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/hmadn_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/02-DNA-rmhost
export outP=$baseP/03-HUMANN-DNA
export databaseP=$baseP/References/databases

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `

echo $SAMPLE

if [ ! -e $outP/$SAMPLE ]; then
    mkdir -p $outP/$SAMPLE
fi


module load humann/3.0.1

source activate metaphlan-3.0.13


call="humann_renorm_table -i $outP/${SAMPLE}/${SAMPLE}_genefamilies.tsv 
      -o $outP/${SAMPLE}/${SAMPLE}_genefamilies.cpm.tsv -u cpm"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

