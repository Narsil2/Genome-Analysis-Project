#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 15:00:00
#SBATCH -J htseq_counts
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load HTSeq

for file in ~/Genome_Analysis/outputs/06-alignment/BH/*.sorted.bam
do
    base=$(basename "$file" .sorted.bam)

    htseq-count -f bam -r pos -i ID -s yes -t CDS \
    "$file" ~/Genome_Analysis/outputs/03-annotation/cleaned.gff \
    > ~/Genome_Analysis/outputs/07-counts/BH/${base}_counts.txt
done

for file in ~/Genome_Analysis/outputs/06-alignment/Serum/*.sorted.bam
do

    base=$(basename "$file" .sorted.bam)

    htseq-count -f bam -r pos -i ID -s yes -t CDS \
    "$file" ~/Genome_Analysis/outputs/03-annotation/cleaned.gff \
    > ~/Genome_Analysis/outputs/07-counts/Serum/${base}_counts.txt
done
