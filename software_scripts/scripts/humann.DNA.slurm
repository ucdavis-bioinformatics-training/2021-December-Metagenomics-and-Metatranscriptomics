#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-12
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/hmad_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/hmad_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/02-DNA-rmhost
export outP=$baseP/03-HUMANN-DNA
export databaseP=/software/humann/3.0.1/lssc0-linux/db

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `

echo $SAMPLE

if [ ! -e $outP/$SAMPLE ]; then
    mkdir -p $outP/$SAMPLE
fi


module load humann/3.0.1

source activate metaphlan-3.0.13

call="cat $seqP/${SAMPLE}/${SAMPLE}_hostrmvd_R1.fastq $seqP/${SAMPLE}/${SAMPLE}_hostrmvd_R2.fastq |gzip - > $outP/${SAMPLE}/${SAMPLE}.fastq.gz"

echo $call
eval $call


call="humann --threads ${SLURM_CPUS_PER_TASK} \
      --input $outP/${SAMPLE}/${SAMPLE}.fastq.gz \
      --pathways metacyc --protein-database $databaseP/uniref --nucleotide-database $databaseP/chocophlan \
      --output $outP/${SAMPLE} --output-basename ${SAMPLE} --o-log  $outP/${SAMPLE}/${SAMPLE}.log"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

