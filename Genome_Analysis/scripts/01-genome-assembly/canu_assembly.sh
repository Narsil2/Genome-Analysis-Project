#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 03:00:00
#SBATCH -J canu_pacbio_assembly
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load canu
module load SAMtools

canu \
	-p E.faecium \
	genomeSize=4m \
	useGrid=false \
	-pacbio-raw /home/yogesh22/Genome_Anlaysis/raw_data/genomics_data/PacBio/*.fastq.gz 1>canu-assembly.stdout 2>canu-assembly.stderr

