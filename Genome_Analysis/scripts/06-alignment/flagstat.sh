#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -J bwa_mem
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load SAMtools

for bam in ~/Genome_Analysis/outputs/06-alignment/BH/*.sorted.bam
do
    samtools flagstat "$bam" > ~/Genome_Analysis/outputs/06-alignment/BH/$(basename "${bam%%.sorted.bam}").flagstat.txt
done    

for bam in ~/Genome_Analysis/outputs/06-alignment/Serum/*.sorted.bam
do
    samtools flagstat "$bam" > ~/Genome_Analysis/outputs/06-alignment/Serum/$(basename "${bam%%.sorted.bam}").flagstat.txt
done    
