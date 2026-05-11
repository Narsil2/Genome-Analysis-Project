#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -J fastqc
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load FastQC

fastqc -o ~/Genome_Analysis/outputs/05-fastqc/BH/ ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_BH/raw/*.fastq.gz

fastqc -o ~/Genome_Analysis/outputs/05-fastqc/serum/ ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_Serum/raw/*.fastq.gz

