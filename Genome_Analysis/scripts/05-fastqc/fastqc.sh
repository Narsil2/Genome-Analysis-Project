#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -J fastqc
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load FastQC

fastqc -o ~/Genome_Analysis/outputs/05-fastqc/trimmed/BH/ /proj/uppmax2026-1-61/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trimmed/*.fastq.gz

fastqc -o ~/Genome_Analysis/outputs/05-fastqc/trimmed/Serum/ /proj/uppmax2026-1-61/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trimmed/*.fastq.gz

