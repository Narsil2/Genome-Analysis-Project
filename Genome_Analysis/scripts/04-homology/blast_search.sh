#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 00:15:00
#SBATCH -J blastn_search
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load BLAST+

blastn -db 'placeholder' \ #fix issue with blast databases in uppmax
	-query /home/yogesh22/Genome_Analysis/outputs/01-genome-assembly/E.faecium.contigs.fasta \
	-out /home/yogesh22/Genome_Analysis/outputs/04-homology/blast.hits

