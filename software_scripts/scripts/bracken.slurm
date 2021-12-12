#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-12
#SBATCH --mem=80000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --output=slurmout/brk_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/brk_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/02-DNA-rmhost
export refP=$baseP/References
export outP=$baseP/03-Kraken


SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `

echo $SAMPLE

if [ ! -e "$outP/$SAMPLE" ]; then
    mkdir -p $outP/$SAMPLE
fi

module load kraken2/2.1.2
module load bracken/2.5

call="bracken -d $refP/krakendb -t 24 -i $outP/${SAMPLE}/${SAMPLE}.kraken_report.out -l S -o $outP/${SAMPLE}/${SAMPLE}_report_species.txt"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

