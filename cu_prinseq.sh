#!/bin/bash
#SBATCH --job-name=cu_prinseq
#SBATCH --cpus-per-task=16
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --output=log/cu_prinseq-%j.out
#SBATCH --error=log/cu_prinseq-%j.err

module load prinseqlite
module load fastqc
mkdir PrinSeq
mkdir ./PrinSeq/Good
mkdir ./PrinSeq/Bad
ls -1 ./Trimmed/*fq.gz|while read filename;
do
bin_name=`basename $filename|cut -f1 -d '.'`
pigz -dc $filename > temp
prinseq-lite.pl -fastq temp -min_qual_mean 30 -out_bad ./PrinSeq/Bad/$bin_name'_bad' -out_good ./PrinSeq/Good/$bin_name'_filt'
fastqc -t 16 $bin_name'_filt.fastq'
pigz -p 16 $bin_name'_filt.fastq'
done
