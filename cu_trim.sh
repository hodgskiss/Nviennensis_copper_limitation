#!/bin/bash
#SBATCH --job-name=cu_trim
#SBATCH --cpus-per-task=16
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --output=log/cu_trim-%j.out
#SBATCH --error=log/cu_trim-%j.err

module load trimmomatic
mkdir Trimmed
ls -1 ./Renamed_Files/*fq|while read filename;
do
bin_name=`basename $filename|cut -f1 -d '.'`
echo $bin_name
trimmomatic SE -phred33 -threads 16 -trimlog R.trim.log $filename Trimmed/$bin_name'_trimmed.fq' ILLUMINACLIP:/apps/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:38 HEADCROP:13
pigz -p 16 Trimmed/$bin_name'_trimmed.fq'
done

