#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -J bwa_mem
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J bwa_mem
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load BWA
module load SAMtools

# Reference (index only once externally ideally)
REF=~/Genome_Analysis/reference-genome/*.fna
PREFIX=Efaecium

# If not already indexed, run once manually:
# bwa index -p $PREFIX $REF

# --- BH samples ---
for r1 in ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_BH/trimmed/*_paired_*_pass_1.fastq.gz
do
    r2=${r1/_pass_1/_pass_2}
    base=$(basename "${r1%%_pass_1.fastq.gz}")

    bwa mem -t 4 $PREFIX "$r1" "$r2" | \
    samtools view -bS | \
    samtools sort -o ~/Genome_Analysis/outputs/06-alignment/BH/${base}.sorted.bam
done


# --- Serum samples ---
for r1 in ~/Genome_Analysis/raw_data/transcriptomics_data/RNA-Seq_Serum/trimmed/*_paired_*_pass_1.fastq.gz
do
    r2=${r1/_pass_1/_pass_2}
    base=$(basename "${r1%%_pass_1.fastq.gz}")

    bwa mem -t 4 $PREFIX "$r1" "$r2" | \
    samtools view -bS | \
    samtools sort -o ~/Genome_Analysis/outputs/06-alignment/Serum/${base}.sorted.bam
done