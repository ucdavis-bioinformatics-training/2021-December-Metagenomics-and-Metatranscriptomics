#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-12
#SBATCH --mem=30000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/btr_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/btr_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export refP=$baseP/References
export outP=$baseP/02-DNA-rmhost

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

module load bowtie2/2.4.2
module load bedtools2/2.29.2
module load samtools/1.11

nbt=$(( ${SLURM_CPUS_PER_TASK} - 4 ))
echo ${nbt}

## -f 12: extract only alignments with both reads unmapped
## -F 256: Do not extract alignments which are not primary alignment
call="bowtie2 -x $refP/GCF_002263795.1_ARS-UCD1.2_genomic \
	-1 <(zcat $seqP/${SAMPLE}/${SAMPLE}_DNA_R1.fastq.gz) -2 <(zcat $seqP/${SAMPLE}/${SAMPLE}_DNA_R2.fastq.gz) \
	-q --phred33 --sensitive --end-to-end \
	-p ${nbt} |samtools view -bh - |samtools view -bh -f 12 -F 256 - |samtools sort -n - | \
	bedtools bamtofastq -i - -fq $outP/${SAMPLE}/${SAMPLE}_hostrmvd_R1.fastq -fq2 $outP/${SAMPLE}/${SAMPLE}_hostrmvd_R2.fastq"


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

