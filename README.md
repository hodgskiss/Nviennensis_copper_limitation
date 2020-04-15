# Nviennensis_copper_limitation
Scripts for the analysis of transcriptomic data from Nitrososphaera viennensis under copper limitation. Information found here is associated with the paper "Genome wide transcriptomic analysis of the soil ammonia oxidizing archaeon *Nitrososphaera viennensis* upon exposure to copper limitation." (Reyes et al. 2020).

##  Table of Contents  

1.  Transcriptomic Processing
    *  Script_Information_for_Processing_Reads
       *  Contains information on how to access transcriptomic data as well as programs and commands used to prepare data for DESeq2 (i.e. create the file cu_counts.txt).
    *  cu_counts.txt
       *  Text file used for input in DESeq2 after processing reads.
    *  DESeq2_Commands.R
       *  Commands used to process data in DESeq2 version 1.10.1.
    *  rlog_mat.txt
       *  Text file containig the rlog normalized reads of genes used to create the PCA plot of samples.
    *  PCA_Plot_Script.R
       *  R script used to create PCA plot of samples.
    
2.  arGOC Analysis
    *  Nvie_arCOGs.csv
       *  File containing the arCOG assignments of genes in *N. viennensis* (NVIE locus tags). 
    *  Up_Genes_NVIE.csv
       *  List of upregulated genes under copper limitation (adjusted p-value &ge; 0.01). 
    *  Down_Genes_NVIE.csv
       *  List of downregulated genes under copper limitation (adjusted p-value &ge; 0.01).
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
       *  Motif found in the promoter region of selected upregulated genes.  Motif is in the MEME format.
    *  Motif_Downregulated_Genes.txt
       *  Motif found in the promoter region of the top 25 downregulated genes.	Motif is in the MEME format.
    
4.  Proteomic Precessing
    *  Proteomic_Data.txt
       *  Information for accessing proteomic data.
    *  copper_mqpar.xml
       *  Parameter file for analyzing proteomic data using MaxQuant version 1.6.3.3.
    
