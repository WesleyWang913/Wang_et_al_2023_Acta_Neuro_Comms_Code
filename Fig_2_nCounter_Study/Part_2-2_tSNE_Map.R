#Unsupervised Clustering of gene reads from nCounter
setwd("Your Directory")
#Load Libraries
library(dplyr)
library(Rtsne)
library(fpc)
library(factoextra)
library(dbscan)
library(ggplot2)
library(ggpubr)
library(Boruta)
library(randomForest)
library(e1071)
library(caret)
library(cluster)
library(grid)
library(gridExtra)
library(pROC)
library(survival)

#Tsne analysis on genes for PD and psPD cases
#Question: If we use unsupervised clustering to separate the genes, do they match based on original disease label?


#Read in normalized counts and clinical data
Master.Data <- read.csv("Normalized Counts from GEO.csv") #Normalized Counts found in GEO GSE231994
rownames(Master.Data) <- Master.Data$Name
Master.Data <- Master.Data[,-1]
Master.Data.t <- t(Master.Data)
Samples <- NULL

#Add Disease Labels
Samples$Data <- Master.Data.t
Samples$Disease <- c(rep("PD",27), rep("psPD",21))
Samples$Disease <- as.factor(Samples$Disease)

#Read in differential genes
DGE <- read.csv("Volcano Plot CSV from Rosalind.csv") #Generated in step when from Rosalind DE analysis
DGE <- DGE[DGE$pvalue <= 0.05,] #Filter off non-significant genes
DGE <- DGE[order(abs(DGE$log2FoldChange), decreasing = T),] #Order by top log fold change
gene.150 <- DGE$Symbol[1:150] #Pull top 150 DEG
sf <- Samples$Data[,colnames(Samples$Data) %in% gene.150]

#Supplementary Figure 2: Supplementary Figure 2: Top 150 DEGs from OSU PD vs psPD Cases. 
#library(pheatmap)
#data_subset <- t(sf)
#my_sample_col <- data.frame(Sample = Samples$Disease)
#row.names(my_sample_col) <- colnames(data_subset)
#pheatmap(data_subset, annotation_col = my_sample_col, fontsize = 6)

#tSNE unsupervised clustering
set.seed(999)
PD.tsne <- Rtsne(sf, dims = 2, perplexity = 15)
tsne.Points <- as.data.frame(PD.tsne$Y)

#Find number of best clusters (2)
fviz_nbclust(tsne.Points, kmeans, method = "silhouette") + theme_classic()
PD.k <- kmeans(tsne.Points, centers = 2, nstart = 25)
summary(as.factor(PD.k$cluster))

#Visualize clusters
fviz_cluster(PD.k, tsne.Points, stand = FALSE, geom = "point", main= "") + theme_classic()+
  xlab("tSNE-1") + ylab("tSNE-2") + theme(legend.position = "none") + scale_fill_manual(values = c("#7da27e", "#D2B48C")) + scale_color_manual(values = c("#7da27e", "#D2B48C"))

#Cluster distribution (Figure 2h)
mf <- data.frame("Cluster" = PD.k$cluster, "Status" = Samples$Disease)
mf$Cluster <- as.factor(mf$Cluster )

mf %>%
  count(Cluster, Status) %>%
  group_by(Cluster) %>%
  mutate(n = n/sum(n) * 100) %>%
  ggplot() + aes(Cluster, n, fill = Status, label = paste0(round(n, 0), "%")) + 
  geom_col() +
  geom_text(position=position_stack(0.5), size = 3.6, color = "white")+
  theme_classic()+
  ylab("Percent of Cluster")+
  theme(legend.position = "none")+ scale_fill_manual(values = c("#FF6A4C", "#487BE3"))

all.tab <- table(PD.k$cluster, Samples$Disease)
all.test <- chisq.test(all.tab)
all.test
all.test$observed
all.test$expected

#Supplementary Figure 6: Cluster statistics of nCounter clustering analysis.
#library(caret)
#library(NbClust)
#library(factoextra)
#library(fpc)
#Generate new dataframe for cluster statistics 
#df <- cbind(sf,PD.k$cluster)
#df <- as.data.frame(sf2)
#colnames(df)[147] <- "Cluster"
#row_numbers <- c(1:nrow(df))
#Make the dummy variables
#df1 <- df[,1:146] 
#dummy1 <- dummyVars(~., data = df1)
#dummy <- predict(dummy1, df1)
#Cluster Stats
#LGG_df <- data.frame()
#for(i in 1:10000)  {
#  my_rows <- sample(row_numbers, size = 40, replace = FALSE)
#  my_d <- c(1, 2)
#  my_cluster_data1 <- cluster.stats(dist(dummy[my_rows,]), clustering = as.integer(df[my_rows, "Cluster"]))
#  d1 <- my_cluster_data1$within.cluster.ss
#  d2 <- my_cluster_data1$average.within
#  d3 <- my_cluster_data1$average.between
#  r_cluster <- sample(my_d, size = 40, replace = TRUE, prob = c(my_cluster_data1$cluster.size[1]/40,
#                                                                my_cluster_data1$cluster.size[2]/40) )
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
#ggplot(LGG_frame, aes(x=Group, y=Average_Within_SS, color=Group))+
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

