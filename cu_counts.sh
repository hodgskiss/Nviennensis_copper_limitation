#!/bin/bash
#SBATCH --job-name=cu_counts
#SBATCH --cpus-per-task=8
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --output=log/counts-%j.out
#SBATCH --error=log/counts-%j.err

module load subread

featureCounts -a /scratch/hodgskiss/Cu_Project/featurecounts/GCA_000698785.1_ASM69878v1_genomic.gff -t gene -g ID -o ./featurecounts/cu_counts.txt /scratch/hodgskiss/Cu_Project/Bowtie2/*.sam
