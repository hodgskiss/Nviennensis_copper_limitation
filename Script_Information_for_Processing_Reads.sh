#!/bin/bash

#The following includes programs and scripts used to prepare transcriptomic data for differential expression analysis of copper limiatation in Nitrososphaera viennensis.
#Results from this analysis are found in the paper "Genome wide transcriptomic analysis of the soil oxidizing archaeon Nitrososphaera viennensis upon exposure to copper limitation" Reyes et al. (2020).

#All raw sequence runs from the RNA-Sequencing analysis were deposited in the DNA Databank of Japan (DDBJ) under the accession number DRA008627.

###Renaming Files
#Renaming files.

#Files were renamed reflecting whether they were under copper replete (R) or copper limited (L) conditions, sample number (1-6, five replicates for each condition), and the day of harvest (Day7 or Day5).
# i.e. sample_name = L2Day7.bam (copper limited, culture number 2, harvested on day 7)

76784_TCCGCGAATATAGCCT_CCVFKANXX_1_20181204B_20181204.bam    ->   L2Day7.bam
76785_TCTCGCGCATAGAGGC_CCVFKANXX_1_20181204B_20181204.bam    ->   L3Day7.bam	
76786_AGCGATAGCCTATCCT_CCVFKANXX_1_20181204B_20181204.bam    ->   L5Day7.bam
76787_ATTACTCGGGCTCTGA_CCVFKANXX_1_20181204B_20181204.bam    ->   L6Day7.bam
76788_TCCGGAGAAGGCGAAG_CCVTGANXX_3_20181204B_20181204.bam    ->   R1Day5.bam
76789_GAGATTCCTAATCTTA_CCVTGANXX_3_20181204B_20181204.bam    ->   R2Day5.bam
76790_ATTCAGAACAGGACGT_CCVTGANXX_3_20181204B_20181204.bam    ->   R3Day5.bam
76791_GAATTCGTGTACTGAC_CCVTGANXX_3_20181204B_20181204.bam    ->   R4Day5.bam
76792_CTGAAGCTTATAGCCT_CCVTGANXX_3_20181204B_20181204.bam    ->   R5Day5.bam
76793_TAATGCGCATAGAGGC_CCVTGANXX_3_20181204B_20181204.bam    ->   L1Day7.bam


#Each command was done for each sample.
#Quality checking via FastQC was performed after Trimmotaic, Prinseq-Lite, and SortMeRNA.
#Note: It is possible to write scripts in a loop format to automatically analyze each sample for each respective step. 

###fastQC version 0.11.5
#Checking the quality of sequenced reads.

fastqc sample_name.bam

###samtools package (Version: 1.9 (using htslib 1.9))
#Converting .bam files to .fq files

samtools fastq --threads 16 sample_name.bam > sample_name.fq

###Trimmomatic Version 0.36
#Cropping reads and trimming out adapter sequences and low quality ends.

trimmomatic SE -phred33 -threads 16 -trimlog R.trim.log sample_name.fq sample_name_trimmed.fq ILLUMINACLIP:/path/to/adapters/adapters.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:38 HEADCROP:13

###PrinSeq-Lite version 0.20.4
#Sorting for a minimum quality score.

prinseq-lite.pl -fastq temp -min_qual_mean 30 -out_bad ./path/to/discarded/sequences -out_good ./path/to/retained/sequences/sample_name_trimmed_filt.fastq

###SortMeRNA version 2.1, 01/02/2016
#Separating rRNA and mRNA reads.  

indexdb_rna --ref /path/to/Nviennensis_rRNAs.fasta,/path/to/index/output

sortmerna --ref /path/to/Nviennensis_rRNAs.fasta,\
/path/to/index/output \
--reads /path/to/reads/sample_name_trimmed_filt.fastq \
--aligned ./rRNA/sample_name_trimmed_filt_rRNA.fastq --other ./mRNA/sample_name_trimmed_filt_mRNA.fastq --fastx --log -a 16

###Bowtie 2 version version 2.3.3.1
#Mapping mRNA reads to the genome of N. viennensis.

bowtie2-build ./path/to/genome_fna ./path/to/indexed/genome

bowtie2 -x /path/to/indexed/genome -U sample_name_trimmed_filt_mRNA.fastq -S ./path/to/aligned_reads/sample_name_trimmed_filt_mRNA.sam

###featurecounts v1.6.2
#Counting the number of mapped reads to each gene. 

featureCounts -a /path/to/genome_gff -t gene -g ID -o ./path/to/output/cu_counts.txt /path/to/aligned_reads/*sample_name_trimmed_filt_mRNA.sam

#The final output is a text file, "cu_counts.txt", that will be used in the DESeq2 differenital analysis software in R. 

###Citations

#fastQC version 0.11.5
#Andrews S. 2010. FastQC: a quality control tool for high throughput sequence data.

#samtools package (Version: 1.9 (using htslib 1.9))
#Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943].

#Trimmomatic Version 0.36
#Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30:2114–2120.

#PrinSeq-Lite version 0.20.4
#Schmieder R, Edwards R. 2011. Quality control and preprocessing of metagenomic datasets. Bioinformatics 27:863–864.

#SortMeRNA version 2.1, 01/02/2016
#Kopylova E, Noé L, Touzet H. 2012. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics 28:3211–3217.

#Bowtie 2 version version 2.3.3.1
#Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2. Nat Methods 9:357–9.

#featurecounts v1.6.2
#Liao Y, Smyth GK, Shi W. 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30:923–930.
