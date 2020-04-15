#Commands used for DESeq2 version 1.10.1 in R version 3.2 with RStudio version 0.99.902.

#Set the pathway for the needed files using setwd(). 

#Load DESeq2.
library(DESeq2)

#Import the featurecounts output as a table in R.
cu_counts <- read.table("cu_counts.txt", header = TRUE, row.names=1)

#Remove unneccesary coloumns in the new table.
cu_counts$Chr <- NULL
cu_counts$Start <- NULL
cu_counts$End <- NULL
cu_counts$Strand <- NULL
cu_counts$Length <- NULL

#Rename the coloumns to something more recognizable (done 10x respectively).
names(cu_counts)[1]<-"L1Day7"

#Assign a condition (i.e. Limited or Replete) to each sample. 
cu_columns <- data.frame(row.names=c("L1Day7","L2Day7","L3Day7","L5Day7","L6Day7","R1Day5","R2Day5","R3Day5","R4Day5","R5Day5"),cu_conc=as.factor(c("Lim","Lim","Lim","Lim","Lim","Rep","Rep","Rep","Rep","Rep")))

#Create a matrix that can be used as an input for DESeq2.
ddsCopper <- DESeqDataSetFromMatrix(countData=cu_counts, colData = cu_columns, design = ~ cu_conc)

#Normalize the data in order to visulize it in a PCA plot.  This is NOT the normalization for the DESeq2 pipeline.  This is an extra step for visualization purposes. 
rlddsCopper <- rlog(ddsCopper)

#Make a text file of the rlog transformed values to use for a PCA plot (see PCA_plot_script.R).
write.table(rlddsCopper, file="rlog_mat.txt", sep = " ")

#Perform the DESeq2 analysis (uses a different normalization than the one used earlier for the PCA plot). 
DEanalysis_cu <- DESeq(ddsCopper)

#Compile results of the Limited condition vs. the Replete Condtion. 
res_Lim_Rep <- results(DEanalysis_cu, contrast=c("cu_conc", "Lim", "Rep"))

#Make a text file of the DESeq2 results. 
write.table(res_Lim_Rep, file="all_res_Lim_Rep.txt", sep = " ")

#Results can be found in Dataset_S1.xlsx in the tab "DESeq2_Output".

###Citations

#DESeq2
#Love MI, Huber W, Anders S. Moderated estimation of fold change and  dispersion for RNA-seq data with DESeq2. Genome Biol 2014; 15: 550. 



