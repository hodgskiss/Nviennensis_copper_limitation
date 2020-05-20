#!/bin/bash
#SBATCH --job-name=sortmerna
#SBATCH --cpus-per-task=16
#SBATCH --mem=50000
#SBATCH --mail-type=ALL
#SBATCH --output=log/sortmerna-%j.out
#SBATCH --error=log/sortmerna-%j.err

module load sortmerna

indexdb_rna --ref /scratch/hodgskiss/Cu_Project/SortMeRNA/Index/Nviennensis_rRNAs.fasta,/scratch/hodgskiss/Cu_Project/SortMeRNA/Index/viennensisALL

ls -1 /scratch/hodgskiss/Cu_Project/PrinSeq/Good/*fastq.gz|while read filename;
do
pigz -d $filename

bin_name=`basename $filename|cut -f1 -d '.'`

sortmerna --ref /scratch/hodgskiss/Cu_Project/SortMeRNA/Index/Nviennensis_rRNAs.fasta,\
/scratch/hodgskiss/Cu_Project/SortMeRNA/Index/viennensisALL \
--reads /scratch/hodgskiss/Cu_Project/PrinSeq/Good/$bin_name'.fastq' \
--aligned ./rRNA/$bin_name'_rRNA' --other ./mRNA/$bin_name'_mRNA' --fastx --log -a 16

pigz -p 16 ./rRNA/$bin_name'_rRNA.fastq'
pigz -p 16 ./mRNA/$bin_name'_mRNA.fastq'

done
