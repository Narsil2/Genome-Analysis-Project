#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -J bwa_mem
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load BWA

bwa index -p Efaecium ~/Genome_Analysis/reference-genome/*.fna

for file in ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_BH/trimmed/*.fastq.gz
do
    bwa mem Efaecium "$file" \
    > ~/Genome_Analysis/outputs/06-alignment/BH/$(basename "${file%.fastq.gz}").sam
done

for file in ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_Serum/trimmed/*.fastq.gz
do
    bwa mem Efaecium "$file" \
    > ~/Genome_Analysis/outputs/06-alignment/Serum/$(basename "${file%.fastq.gz}").sam
done
