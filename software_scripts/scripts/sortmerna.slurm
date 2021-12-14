#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-12
#SBATCH --mem=30000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/srnas_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/srnas_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export refP=$baseP/References
export outP=$baseP/02-mRNA-rRNArmvd-test

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `
TYPE=$1

echo $SAMPLE
echo $TYPE
export seqP=$baseP/01-HTS_Preproc/$TYPE

if [ ! -e $outP ]; then
    mkdir -p $outP
fi

if [ ! -e "$outP/$SAMPLE" ]; then
    mkdir -p $outP/$SAMPLE
fi

module load sortmerna/4.3.4

if [ ${SLURM_ARRAY_TASK_ID} == 2 ]
then
  call="sortmerna --ref $refP/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta \
	--ref $refP/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
	--ref $refP/RF00001.fa --ref $refP/RF00002.fa \
	--reads $seqP/${SAMPLE}/${SAMPLE}_mRNA_R1.fastq.gz \
	--reads $seqP/${SAMPLE}/${SAMPLE}_mRNA_R2.fastq.gz \
	--workdir $outP/${SAMPLE} --fastx --threads ${SLURM_CPUS_PER_TASK} \
	--other $outP/${SAMPLE}/${SAMPLE}"
else
  call="sortmerna --ref $refP/SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta \
	--ref $refP/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
	--ref $refP/RF00001.fa --ref $refP/RF00002.fa \
	--reads $seqP/${SAMPLE}/${SAMPLE}_mRNA_SE.fastq.gz \
	--workdir $outP/${SAMPLE} --fastx --threads ${SLURM_CPUS_PER_TASK} \
	--other $outP/${SAMPLE}/${SAMPLE}"
fi

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

