#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=0-10
#SBATCH --mem=20000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --reservation=meta_workshop
#SBATCH --account=workshop
#SBATCH --partition=production
#SBATCH --output=slurmout/hts_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/hts_%A_%a.err # File to which STDERR will be written


start=`date +%s`
hostname

export baseP=/share/workshop/meta_workshop/$USER/meta_example
export seqP=$baseP/00-RawData
export cwd=$baseP/scripts

SAMPLE=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -1 `
TYPE=$1

echo $SAMPLE
echo $TYPE
export outdir=$baseP/01-HTS_Preproc/$TYPE

if [ ! -e $outdir ]; then
    mkdir -p $outdir
fi

if [ ! -e "$outdir/$SAMPLE" ]; then
    mkdir -p $outdir/$SAMPLE
fi

module load htstream/1.3.2

if [ $TYPE == "DNA" ]
then
  call="hts_Stats -L $outdir/$SAMPLE/$SAMPLE.stats.log \
      -1 $seqP/${SAMPLE}_${TYPE}_1.fastq.gz -2 $seqP/${SAMPLE}_${TYPE}_2.fastq.gz | \
      hts_SeqScreener -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_AdapterTrimmer -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_QWindowTrim -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_LengthFilter -m 50 -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_Stats -f $outdir/$SAMPLE/${SAMPLE}_${TYPE}.cleaned -A $outdir/$SAMPLE/$SAMPLE.stats.log"
else
  call="hts_Stats -L $outdir/$SAMPLE/$SAMPLE.stats.log \
      -1 $seqP/${SAMPLE}_${TYPE}_1.fastq.gz -2 $seqP/${SAMPLE}_${TYPE}_2.fastq.gz | \
      hts_SeqScreener -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_AdapterTrimmer -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_Overlapper -o 10 -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_QWindowTrim -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_LengthFilter -m 30 -A $outdir/$SAMPLE/$SAMPLE.stats.log | \
      hts_Stats -f $outdir/$SAMPLE/${SAMPLE}_${TYPE}.cleaned -A $outdir/$SAMPLE/$SAMPLE.stats.log"
fi


echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo Runtime: $runtime seconds

