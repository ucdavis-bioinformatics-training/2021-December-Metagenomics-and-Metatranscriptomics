#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-10
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/btb_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/btb_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export refP=$baseP/References

module load bowtie2/2.4.2

call="bowtie2-build --threads 12 $refP/GCF_002263795.1_ARS-UCD1.2_genomic.fna $refP/GCF_002263795.1_ARS-UCD1.2_genomic"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

