#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 00:15:00
#SBATCH -J eval_quast
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load QUAST/5.2.0-foss-2024a-Python-2.7.18

python /sw/arch/eb/software/QUAST/5.3.0-gfbf-2024a/bin/quast ~/Genome_Anlaysis/outputs/01-genome-assembly/E.faecium.contigs.fasta \
	-r /home/yogesh22/Genome_Anlaysis/reference-genome/GCF_009734005.1_ASM973400v2_genomic.fna 
