#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=0-12
#SBATCH --mem=80000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/mtp_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/mtp_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/02-mRNA-rRNArmvd
export outP=$baseP/03-Metaphlan-RNA-test
export databaseP=$baseP/References/databases

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `

echo $SAMPLE

if [ ! -e $outP/$SAMPLE ]; then
    mkdir -p $outP/$SAMPLE
fi


module load metaphlan/3.0.13

source activate metaphlan-3.0.13


call="metaphlan --nproc ${SLURM_CPUS_PER_TASK} \
      --tmp_dir $outP/${SAMPLE}/tmp \
      --input_type fastq $seqP/${SAMPLE}/${SAMPLE}.fq.gz $outP/${SAMPLE}/${SAMPLE}_profile.txt"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

