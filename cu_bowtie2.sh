#!/bin/bash
#SBATCH --job-name=cu_bowtie2
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --output=log/indexgenome-%j.out
#SBATCH --error=log/indexgenome-%j.err

module load bowtie2

bowtie2-build ./Bowtie2/Genome_Index/GCA_000698785.1_ASM69878v1_genomic.fna ./Bowtie2/Genome_Index/n_viennensis

ls -1 /scratch/hodgskiss/Cu_Project/SortMeRNA/mRNA/*_mRNA.fastq.gz|while read filename;

do

pigz -dc $filename > temp

bin_name=`basename $filename|cut -f1 -d '.'`


bowtie2 -x /scratch/hodgskiss/Cu_Project/Bowtie2/Genome_Index/n_viennensis -U temp -S ./Bowtie2/$bin_name'.sam'

rm temp

done
