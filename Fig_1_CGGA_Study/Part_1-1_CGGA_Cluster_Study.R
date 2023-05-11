###CAnalysis for CGGA Gene Cluster and Survival Analysis (Figure 1)
setwd("Your Directory") #Set your file directory

#Add Libraries
library(Rtsne)
library(dplyr)
library(reshape2)  
library(factoextra)
library(ggpubr)
library(Boruta)
library(randomForest)
library(e1071)
library(caret)
library(cluster)

#Tsne analysis on genes for recurrent cases
#Define genes to subset (Supplementary Table 1)
Target.Genes <- c("PROM1", "FUT4", "ITGA6", "CD44", "L1CAM", "SOX2", "NANOG", "OLIG2", "MYC", "BMI1", "MSI1", "NES", "MKI67", "POLD1" , "POLD2", "POLD3", "POLD4", "PCNA", "CCND1", "CCNA2", "ID1", "CD3E", "CD3D", "CD3G", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "CCL5", "PRF1", "CIITA", "JAK2", "GBP5", "IRF2", "GBP4", "IRF1", 
                  "HCST", "IL16", "IL18", "CCL4", "CCL3", "CCL19", "CCL22", "CXCR3", "BCL11B", "IL7R", "KLRG1", "NKG7", "SAMD3", "STYK1", "TBX21", "BANK1", "CD79A", "CD55", "CD19", "CD79B", "CD38", "IGFBP4", "ITM2A", "AMIGO2", "TRAT1", "CD40LG", "ICOS", "NR4A3", "HAVCR2", "KMO", "DNASE1L3", "ANPEP", "CXCL16", "C1QC", "CD5L", "FCGR3A",
                  "ITGB5", "MERTK", "CCL8", "IL2RA", "CTLA4", "FOXP3", "SLC35D1", "GDPD3", "CISH")
length(Target.Genes) #83 Genes

#Read in CGGA Gene and Clinical Data
Genes.CGGA <- read.delim("CGGA Expression Data.csv") #Download from CGGA website (http://www.cgga.org.cn/download.jsp) mRNAseq_693
Filtered.Genes <- Genes.CGGA[Genes.CGGA$Gene_Name %in% Target.Genes,]
rownames(Filtered.Genes) <- Filtered.Genes$Gene #Fix gene label column and remove
Filtered.Genes <- Filtered.Genes[,-1]
Clinical.CGGA <- read.delim("CGGA Clinical Data.csv")  #Download from CGGA website

#Pull rGBM Cases
Clinical.CGGA <- Clinical.CGGA[!is.na(Clinical.CGGA$Histology),]
Recurrent.GBs <- Clinical.CGGA[Clinical.CGGA$Histology == "rGBM",]
recurrent.genes <- Filtered.Genes[,Recurrent.GBs$CGGA_ID]
ncol(recurrent.genes) #109 Cases

#GBM tSNE
set.seed(888)
Recurrent.GBM.Genes <- recurrent.genes
GBM.Recur.Data <- t(Recurrent.GBM.Genes)
GBM.tsne <- Rtsne(GBM.Recur.Data, check_duplicates = FALSE)
GBM.tsne.Points <- as.data.frame(GBM.tsne$Y)

#Identify number of best clusters to identify with kmeans
fviz_nbclust(GBM.tsne.Points, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2) #Picking clusters 3

#kMeans Cluster
GBM.k <- kmeans(GBM.tsne.Points, centers = 3, nstart = 25)
fviz_cluster(GBM.k, GBM.tsne.Points, stand = FALSE, geom = "point", main= "") + theme_classic() + scale_fill_manual(values = c("#20639B", "#1D7874", "#ED553B")) + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B"))+
  xlab("tSNE-1") + ylab("tSNE-2") + theme(legend.position = "none") 

GBM.clin <- data.frame("CGGA_ID" =  colnames(Recurrent.GBM.Genes), "Cluster" = GBM.k$cluster)
GBM.clin <- merge(GBM.clin, Recurrent.GBs, by = "CGGA_ID")

#Suvival Analysis
library(survival)
library(survminer)

fit_gbm <- survfit(Surv(OS, Censor..alive.0..dead.1.) ~Cluster, data = GBM.clin)

#Final Figures
fviz_cluster(GBM.k, GBM.tsne.Points, stand = FALSE, geom = "point", main= "") + theme_classic() + scale_fill_manual(values = c("#20639B", "#1D7874", "#ED553B")) + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B"))+
  xlab("tSNE-1") + ylab("tSNE-2") + theme(legend.position = "none")
ggsurvplot(fit_gbm, data = GBM.clin, pval = TRUE, pval.coord = c(425,0.52), palette = c("#20639B", "#1D7874", "#ED553B"),
                         xlim = c(0,1300),surv.median.line = 'hv', xscale = 30, break.x.by = 150, risk.table = F, legend="none")+xlab("Time (Months)") 

#Immune Deconvolution
rownames(Genes.CGGA) <- Genes.CGGA$Gene_Name
Genes.CGGA <- Genes.CGGA[,-1]
rec_genes <- Genes.CGGA[,Recurrent.GBs$CGGA_ID]
#write.csv(rec_genes, "Recurrent_gene_matrix.csv") #Generate a dataframe of genes to imput into xCell

#Pipelined into TxCell
immune_gbm <- read.delim("xCell_Recurrent_immune.txt") #Included in files
rownames(immune_gbm) <- immune_gbm$X
immune_gbm <- immune_gbm[,-1]

immune_gbm <- data.frame(t(immune_gbm))
immune_gbm$CGGA_ID <- rownames(immune_gbm)
immune_gbm <- merge(immune_gbm, GBM.clin, by = "CGGA_ID")
immune_gbm$Cluster <- as.factor(immune_gbm$Cluster)

#Filter off irrelevant data columns
immune_shape <- subset(immune_gbm, select = -c(CGGA_ID, PRS_type, Histology, Grade, Gender, Age, OS, Censor..alive.0..dead.1., Radio_status..treated.1.un.treated.0., Chemo_status..TMZ.treated.1.un.treated.0., 
                               IDH_mutation_status, X1p19q_codeletion_status, MGMTp_methylation_status, Adipocytes, CD4..Tem, CD4..Tcm, CD8..Tem, CD8..Tcm,
                               Osteoblast, Smooth.muscle, Skeletal.muscle, Sebocytes, Preadipocytes, StromaScore, Pericytes, MPP, Hepatocytes, Chondrocytes,
                               Astrocytes, Keratinocytes, Mesangial.cells, pro.B.cells, mv.Endothelial.cells, ly.Endothelial.cells, Melanocytes, NK.cells, Neurons, 
                               Epithelial.cells, Tgd.cells, Myocytes, Fibroblasts, Erythrocytes, cDC, pDC, GMP, MEP, Megakaryocytes, Platelets))

immune_long <- melt(immune_shape, id = "Cluster")                      
head(immune_long)   

#Fix labels
levels(immune_long$variable) <- c("B Cells", "Basophils", "CD4 T Cells", "Memory CD4 T Cells", "Naive CD4 T Cells", "CD8 T Cells", "Naive CD8 T Cells", 
                                  "Common Lymphoid Progenitor", "Common Myeloid Progenitor", "Class Switched-Memory B Cells", "Dendritic Cells", "Endothelial Cells",
                                  "Eosinophils", "Hematopoietic Stem Cells", "Mesenchymal Stem Cells",
                                  "Macrophages", "M1 Macrophages", "M2 Macrophages", "Mast Cells", "Memory B Cells", "Monocytes", "Natural Killer T Cells", "Neutrophils", 
                                  "Plasma Cells", "Th1 Helper T Cells", "Th2 Helper T Cells", "Tregs", "Activated Dendritic Cells", "Immature Dendritic Cells", "Naive B Cells", "Immune Score", 
                                  "Microenvironment Score")

immune_long$variable <- factor(immune_long$variable, levels = rev(c("Microenvironment Score", "Immune Score", "Endothelial Cells", "Mesenchymal Stem Cells",
                                                                "Hematopoietic Stem Cells","Common Lymphoid Progenitor", "Common Myeloid Progenitor", 
                                                                "CD8 T Cells", "Naive CD8 T Cells", "Natural Killer T Cells", "CD4 T Cells", "Naive CD4 T Cells", "Memory CD4 T Cells",
                                                                "Th1 Helper T Cells", "Th2 Helper T Cells", "Tregs", "B Cells", "Naive B Cells", "Memory B Cells",  "Class Switched-Memory B Cells",
                                                                "Plasma Cells", "Macrophages", "M1 Macrophages", "M2 Macrophages", "Monocytes",
                                                                "Mast Cells", "Eosinophils", "Basophils", "Neutrophils", "Dendritic Cells", "Activated Dendritic Cells", "Immature Dendritic Cells")))



ggplot(immune_long, aes(x = variable, y = value, color = Cluster)) +  # ggplot function
  geom_boxplot(outlier.shape = NA) + coord_flip() + stat_compare_means(method = "anova", label = "p.signif") +
  ylab("Enrichment Score") + xlab("") + theme_classic()   + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) + theme(legend.position = "none")

#Plot significant scores
my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )

ggplot(immune_gbm, aes(x= Cluster, y= MicroenvironmentScore, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("Microenvironment Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggplot(immune_gbm, aes(x= Cluster, y= ImmuneScore, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("Immune Score")  + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggplot(immune_gbm, aes(x= Cluster, y= CLP, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("CLP Enrichment")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggplot(immune_gbm, aes(x= Cluster, y= Macrophages, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("Macrophage Enrichment")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#Supplementary Figure 4: Cluster statistics of CGGA clustering analysis.
#library(caret)
#library(NbClust)
#library(factoextra)
#library(fpc)
#GBM.Recur.Data <- cbind(GBM.Recur.Data,GBM.k$cluster)
#GBM.Recur.Data <- data.frame(GBM.Recur.Data)
#colnames(GBM.Recur.Data)[84] <- "Cluster"
#df <- GBM.Recur.Data
#row_numbers <- c(1:nrow(GBM.Recur.Data))
#making the dummy variables
#df1 <- df[,1:83] #For CGGA
#dummy1 <- dummyVars(~., data = df1)
#dummy <- predict(dummy1, df1)
#note there are 3 clusters, so the random thing is doing the three clusters
#LGG_df <- data.frame()
#for(i in 1:10000)  {
#  my_rows <- sample(row_numbers, size = 50, replace = FALSE)
#  my_d <- c(1, 2, 3)
#  my_cluster_data1 <- cluster.stats(dist(dummy[my_rows,]), clustering = as.integer(df[my_rows, "Cluster"]))
#  d1 <- my_cluster_data1$within.cluster.ss
#  d2 <- my_cluster_data1$average.within
#  d3 <- my_cluster_data1$average.between
#  r_cluster <- sample(my_d, size = 50, replace = TRUE, prob = c(my_cluster_data1$cluster.size[1]/50,
#                                                                my_cluster_data1$cluster.size[2]/50,
#                                                                my_cluster_data1$cluster.size[3]/50) )
#  random_c <- cluster.stats(dist(dummy[my_rows, ]), clustering = r_cluster)
#  r1 <- random_c$within.cluster.ss
#  r2 <- random_c$average.within
#  r3 <- random_c$average.between
#  perm_data <- data.frame("actual_within_cluster_ss" = d1, 
#                          "actual_average_within" =d2,
#                          "actual_average_between" = d3, 
#                          "random_within_cluster_ss" = r1, 
#                          "random_average_within" =r2,
#                          "random_average_between" = r3)
#  LGG_df <- rbind(LGG_df, perm_data)
#  print(i)
#}
#LGG_frame_actual <- data.frame(LGG_df$actual_within_cluster_ss, LGG_df$actual_average_within, LGG_df$actual_average_between, rep("Actual", 10000))
#colnames(LGG_frame_actual) <- c("Within_Cluster_SS", "Average_Within_SS", "Average_Between_SS", "Group")
#LGG_frame_random <- data.frame(LGG_df$random_within_cluster_ss, LGG_df$random_average_within, LGG_df$random_average_between, rep("Random", 10000))
#colnames(LGG_frame_random) <- c("Within_Cluster_SS", "Average_Within_SS", "Average_Between_SS", "Group")
#LGG_frame <- rbind(LGG_frame_actual, LGG_frame_random)
#my_comparisons <- list( c("Actual", "Random") )
# ggplot(LGG_frame, aes(x=Group, y=Average_Within_SS, color=Group))+
#  geom_boxplot()+
#  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") +
#  theme_classic() + theme(legend.position = "none")  + ylab("Within Cluster Distance")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
#  xlab("")  + scale_color_manual(values = c("chartreuse4", "azure4"))+
#  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
#ggplot(LGG_frame, aes(x=Group, y=Average_Between_SS, color=Group))+
#  geom_boxplot()+
#  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") +
#  theme_classic() + theme(legend.position = "none")  + ylab("Between Cluster Distance")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))+
#  xlab("")  + scale_color_manual(values = c("chartreuse4", "azure4"))+
#  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))



