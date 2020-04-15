#Set the pathway for the needed files using setwd(). 

arCOGs <- read.csv(file="Nvie_arCOGs.csv",header=TRUE,sep=",",na.strings="")
#Nvie_arCOGs.csv is a list of arCOG categories identifed in proteins from N. viennensis.

up_names <- read.csv(file="Up_Genes_NVIE.csv",header=TRUE,sep=",")
#Up_Genes_NVIE.csv is a list of locus tags for genes that were found to be differentially upregulated
#based on DESeq2 analysis and a p-value of 0.01.
up_name_vec = as.character(up_names[,1])


down_names <- read.csv(file="Down_Genes_NVIE.csv",header=TRUE,sep=",")
#Down_Genes_NVIE.csv is a list of locus tags for genes that were found to be differentially downregulated
#based on DESeq2 analysis and a p-value of 0.01.
down_name_vec = as.character(down_names[,1])

#These inputs are used for the script to perform the hypergeometric test. 
#Original script can be found at https://github.com/amyschmid/histone_arCOG.

#Required edits for N. viennensis:
#   In the original script, line 26 needs to read: resm <- matrix(0,22,3)
#     -The middle value must be 22 for N. viennensis. 
#   Delete lines 4-12.
#   Use source("name of the script") to activate the functions.

#For this study, namelist=up_name_vec or down_name_vec, cogfile=arCOGS, and pvalue=0.05.
#cogtest2(namelist, cogfile, pvalue) will give the functional categories that are enriched at a p-value of 0.05.

up_enriched = cogtest2(up_name_vec,arCOGs,0.05)
down_enriched = cogtest2(down_name_vec,arCOGs,0.05)

write.table(up_enriched, file="Upregulated_Enriched_Categories.txt", sep = " ")
write.table(down_enriched, file="Downregulated_Enriched_Categories.txt", sep = " ")


#To get the full list or arCOGs for up and downregulated genes, use the following commands.
#up_clust and down_clust represents the full lists for each.
up_cogs= subset(arCOGs, is.element(arCOGs$GeneName, up_name_vec)==TRUE)
up_clust= as.data.frame(summary(up_cogs$funclass_name))
down_cogs= subset(arCOGs, is.element(arCOGs$GeneName, down_name_vec)==TRUE)
down_clust= as.data.frame(summary(down_cogs$funclass_name))

write.table(up_clust, file="arCOG_counts_for_Upregulated_Genes.txt", sep = " ")
write.table(down_clust, file="arCOG_counts_for_Downregulatd_Genes.txt", sep = " ")
