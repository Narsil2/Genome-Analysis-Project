#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 06:00:00
#SBATCH -J trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load Trimmomatic

ADAPTER="/sw/generic/pixi-envs/shovill-1.4.2/.pixi/envs/default/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa"

INPUT_DIR=/proj/uppmax2026-1-61/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/raw
OUTDIR=~/Genome_Analysis/outputs/05.2-trimming/BH/

for SAMPLE in ERR1797972 ERR1797973 ERR1797974
do
  READ1=${INPUT_DIR}/${SAMPLE}_1.fastq.gz
  READ2=${INPUT_DIR}/${SAMPLE}_2.fastq.gz

  OUT_P1=${OUTDIR}/trim_paired_${SAMPLE}_1.fastq.gz
  OUT_S1=${OUTDIR}/trim_single_${SAMPLE}_1.fastq.gz
  OUT_P2=${OUTDIR}/trim_paired_${SAMPLE}_2.fastq.gz
  OUT_S2=${OUTDIR}/trim_single_${SAMPLE}_2.fastq.gz

  trimmomatic PE -threads 1 \
    $READ1 $READ2 \
    $OUT_P1 $OUT_S1 \
    $OUT_P2 $OUT_S2 \
    ILLUMINACLIP:${ADAPTER}:2:30:10\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


