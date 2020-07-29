# Nviennensis_copper_limitation
Scripts for the analysis of transcriptomic data from Nitrososphaera viennensis under copper limitation. Information found here is associated with the paper "Genome wide transcriptomic analysis of the soil ammonia oxidizing archaeon *Nitrososphaera viennensis* upon exposure to copper limitation." (Reyes et al. 2020) https://doi.org/10.1038/s41396-020-0715-2.

##  Table of Contents  

1.  Transcriptomic Processing
    *  Workflow_Information_for_Processing_Reads
       *  Contains information on how to access transcriptomic data as well as programs and commands used to prepare data for DESeq2 (i.e. create the file cu_counts.txt).
    *  cu_trim.sh
       *  Script for trimming adapter sequences from raw reads using Trimmomatic.
    *  cu_prinseq.sh
       *  Script for filtering reads based on quality using PrinSeq-Lite.
    *  cu_sortmerna.sh
       *  Script for separating mRNA reads using SortMeRNA.
    *  Nviennensis_rRNAs.fasta
       *  Reference file for rRNA genes found in *N. viennensis*.  Used with SortMeRNA.
    *  cu_bowtie2.sh
       *  Script for aligning reads to a genome using Bowtie2.
    *  GCA_000698785.1_ASM69878v1_genomic.fna
       *  *N. viennensis* genome used for Bowtie2.
    *  cu_counts.sh
       *  Script for counting reads matched genes using featurecounts.
     * GCA_000698785.1_ASM69878v1_genomic.gff
       *  *N. viennensis* annotation file used for featurecounts.
    *  cu_counts.txt
       *  Text file used for input in DESeq2 after processing reads.
    *  DESeq2_Commands.R
       *  Commands used to process data in DESeq2 version 1.10.1.
    *  rlog_mat.txt
       *  Text file containig the rlog normalized reads of genes used to create the PCA plot of samples.
    *  PCA_Plot_Script.R
       *  R script used to create PCA plot of samples.
    
   
    
2.  arCOG Analysis
    *  Nvie_arCOGs.csv
       *  File containing the arCOG assignments of genes in *N. viennensis* (NVIE locus tags). 
    *  Up_Genes_NVIE.csv
       *  List of upregulated genes under copper limitation (adjusted p-value &lt; 0.01). 
    *  Down_Genes_NVIE.csv
       *  List of downregulated genes under copper limitation (adjusted p-value &lt; 0.01).
    *  Inputs_for_Hypergeometric_Text.R
       *  Script and information for creating the necessary inputs for the arCOG hypergeometric test analysis. 
    
3.  Motif Analysis
    *  Promoter_Regions_Top_25_Upregulated_Genes.txt
       *  Promoter regions of th top 25 upregulated genes under copper limitation in N. viennensis.  Used as input for MEME and FEME.
    *  Promoter_Regions_Top_25_Downregulted_Genes.txt
       *  Promoter regions of the top 25 downregulated genes under copper limitation in N. viennensis.  Used as input for MEME.
    *  Promoter_Regions_Selected_Upregulated_Genes.txt
       *  Promoter regions of selected upregulated genes under copper limitation in N. viennensis.  Used as input for MEME.
    *  Promoter_Regions_Downregulated_Putative_Homologs.txt
       *  Promoter regions of genes in closely related species that share high similarity and/or synteny with genes containing the downregulated motif in N. viennensis.  Used as input for MEME. 
    *  Motif_Upregulated_Genes.txt
       *  Motif found in the promoter region of selected upregulated genes.  Motif is in the MEME format.  Used as input for FEME.
    *  Motif_Downregulated_Genes.txt
       *  Motif found in the promoter region of the top 25 downregulated genes.	Motif is in the MEME format.
    
4.  Proteomic Precessing
    *  Proteomic_Data_Information.txt
       *  Information for accessing proteomic data.
    *  copper_mqpar.xml
       *  Parameter file for analyzing proteomic data using MaxQuant version 1.6.3.3.
    *  uniprot-proteome%3AUP000027093+AND+%28proteomecomponent%3AChromosome%29.fasta
       *  Proteome for *N. viennensis* used for MaxQuant analysis.
    
