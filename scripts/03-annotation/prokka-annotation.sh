#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 00:15:00
#SBATCH -J job_name
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

#load module
module load prokka

#run prokka
prokka --outdir ~/Genome_Analysis/outputs/03-annotation/ --prefix E.faecium-annotation \
        ~/Genome_Analysis/outputs/01-genome-assembly/E.faecium.contigs.fasta
