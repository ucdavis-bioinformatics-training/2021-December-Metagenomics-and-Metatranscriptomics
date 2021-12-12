#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-12
#SBATCH --mem=80000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/mgh_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/mgh_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/02-DNA-rmhost
export outP=$baseP/03-Megahit-DNA

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `

echo $SAMPLE

if [ ! -e $outP ]; then
    mkdir -p $outP
fi

# do not create $outP/${SAMPLE}, megahit fails because of existing dir

module load megahit/1.2.9

call="megahit -t ${SLURM_CPUS_PER_TASK} -m 0.4 \
      -1 $seqP/${SAMPLE}/${SAMPLE}_hostrmvd_R1.fastq -2 $seqP/${SAMPLE}/${SAMPLE}_hostrmvd_R2.fastq \
      -o $outP/${SAMPLE} --out-prefix ${SAMPLE}"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

