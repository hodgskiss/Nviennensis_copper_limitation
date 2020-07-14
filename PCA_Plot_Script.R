#Set the pathway for the needed files using setwd(). 

#Load needed packages.
library(ggplot2)
library(factoextra)
library(ggrepel)

#Load the rlog transformed data matrix made using DESeq2. 
rlog <- read.csv(file="rlog_mat.txt", header=TRUE, sep="")

#Transform matrix to make each column representative of one gene.
rlog_t <- t(rlog)

#Remove columns (aka genes) with 0 reads in every sample.
rlog_t_edit <- rlog_t[,c(colSums(rlog_t) !=0)]

#Rename rows of matrix to reflect experimental conditions.
row_names <-c("L1","L2","L3","L5","L6","R1","R2","R3","R4","R5")
row.names(rlog_t_edit) <- row_names

#Calculate principal components from rlog dataset. 
pca <- prcomp(rlog_t_edit, scale. = TRUE)
pca_df <- as.data.frame(pca$x)

#Label conditions as limited or replete. 
Condition <- c(rep("Lim", times=5),rep("Rep",times=5))
pca_df$Condition <- Condition

sum <- summary(pca)
var <- round(sum$importance[2,]*100,2)

#Create labels for the x-axis and y-axis.
percentage_label <- paste( colnames(pca_df), "(", paste( as.character(var), "% of explained variance", ")", sep="") )

#Plot the PCA using ggplot. 
plot <- ggplot(pca_df, aes(x=PC1,y=PC2,color=Condition))
plot <- plot + geom_point()
plot <- plot + ggtitle("Copper Conditions PCA Plot") + xlab(percentage_label[1]) + ylab(percentage_label[2])
plot <- plot + theme(plot.title=element_text(hjust=0.5))
plot <- plot + geom_label_repel(aes(label=rownames(pca_df)),point.padding = 0.7)
plot

#Open the PCA plot in a separate window. 
  x11()
  fviz <- fviz_pca_ind(pca, habillage = pca_df$Condition, invisible = "quali", repel = TRUE, show.legend=FALSE, addEllipses=TRUE, ellipse.level=0.95)+ggtitle("PCA of Copper Limited and Copper Replete Conditions") + xlab(percentage_label[1]) + ylab(percentage_label[2])+ theme(plot.title=element_text(hjust=0.5))
  fviz
fviz+geom_text_repel(segment.color=NA)
fviz

###Citations

#ggplot2
#H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

#factoextra
#Alboukadel Kassambara and Fabian Mundt (2017). factoextra: Extract and Visualize the Results of Multivariate Data Analyses. R package version 1.0.5. https://CRAN.R-project.org/package=factoextra

#ggrepel
#Kamil Slowikowski (2019). ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'. R package version 0.8.1. https://CRAN.R-project.org/package=ggrepel



