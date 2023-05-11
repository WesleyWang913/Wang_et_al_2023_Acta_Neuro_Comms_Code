#Analysis of computed features from he Segmentation
setwd("Your Image Directory")

#load libraries (missing libraries are in part 3-1)
library(ggpubr)
library(stringr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(sjstats)
library(emmeans)
library(reshape)
library(rstatix)

##Denoise Images ###########Note: You may skip this step and jump to visualization below############
#Stain Segments
features_list <- list.files(path = "features") #Set path to where your shape features are from the deconvolution
length(features_list) #520 Images

#Pull h features
features_list <- features_list[grepl("_h", features_list)]
length(features_list) #260 Images

#Distance Mapping
filt_list <- features_list
length(filt_list)

###Note: There is some more noise in some images > Cut off based on area > 250 (punctate)
segment.distances <- NULL
segment.mf <- NULL

for (b in filt_list) {
  print(b)
  feats <- read.csv(b) #Read in first image
  feats <- feats[,-1] #remove segment ID
  feats <- feats[feats$R.0.s.area > 250, ]
  feats$Case <- str_extract(b, "[^_]+") #Store Case ID
  feats$Image <- b #Store Image
  segment.mf <- rbind(segment.mf, feats)
  points <- data.frame("X" =feats$R.0.m.cx, "Y"= feats$R.0.m.cy)
  my_dist <- dist(points)
  my_matrix <- as.matrix(my_dist)
  my_matrix1 <- apply(my_matrix, 2, sort)
  distance <- c()
  for(i in 1:ncol(my_matrix1)){
    d50 <- my_matrix1[51,i]
    d100 <- my_matrix1[101,i]
    dist_map <- data.frame("d50" = d50, "d100" = d100)
    distance <- rbind(dist_map, distance)
  }
  dist_mean <- apply(distance, 2, mean)
  dist_median <- apply(distance, 2, median)
  dist_sd <- apply(distance, 2, sd)
    
  dist_mean <- as.data.frame(t(dist_mean))
  dist_median <- as.data.frame(t(dist_median))
  dist_sd <- as.data.frame(t(dist_sd))
    
  names(dist_mean) <- paste("mean", names(dist_mean), sep = "")
  names(dist_median) <- paste("median", names(dist_median), sep = "")
  names(dist_sd) <- paste("sd", names(dist_sd), sep = "")
    
  dist_df <- cbind(dist_mean,dist_median, dist_sd)
  dist_df$Case <- str_extract(b, "[^_]+")
  dist_df$Image <- b
  segment.distances <- rbind(segment.distances, dist_df)
  }

#Start here for reanalysis
segment.distances <- read.csv("he_Segment_Distances.csv")

#Case label based on internal lab labels, please disregard from GEO and Supplement Labels
PD.cases <- c(4, 7, 12, 23)
psPD.cases <- c(36,53, 50, 56)
PD.distances <- segment.distances[segment.distances$Case %in% PD.cases,]
nrow(PD.distances)
PD.distances$Status <- "PD"
psPD.distances <- segment.distances[segment.distances$Case %in% psPD.cases,]
nrow(psPD.distances)
psPD.distances$Status <- "psPD"

PD.distances$Histology <- gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", PD.distances$Image)
PD.distances$Subtype <- paste0(PD.distances$Status, "-", PD.distances$Histology)
PD.distances$Subtype <- factor(PD.distances$Subtype)
levels(PD.distances$Subtype) <- c("PD-Control", "PD-Inflammatory", "PD-Hypercellular")
psPD.distances$Histology <- gsub("^(?:[^_]+_){1}([^_]+).*", "\\1", psPD.distances$Image)
psPD.distances$Subtype <- paste0(psPD.distances$Status, "-", psPD.distances$Histology)
psPD.distances$Subtype <- factor(psPD.distances$Subtype)
levels(psPD.distances$Subtype) <- c("psPD-Control", "psPD-Inflammatory", "psPD-Hypercellular")
df_density <- rbind(PD.distances, psPD.distances)
df_density$Subtype <- factor(df_density$Subtype, levels = c("PD-Control","psPD-Control",
                                                                              "PD-Hypercellular", "psPD-Hypercellular",
                                                                              "PD-Inflammatory","psPD-Inflammatory"))

df_density$Histology <- str_split_i(df_density$Subtype, "-", 2)

#Visualization
my_comparisons <- list( c("PD-Control", "psPD-Control"), c("PD-Hypercellular", "psPD-Hypercellular"), c("PD-Inflammatory", "psPD-Inflammatory") )
ggplot(df_density, aes(x = Subtype, y = meand100, color = Subtype)) +
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("NN-100 Distance") +
  xlab("")  + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

my_comparisons <- list (c("PD", "psPD"))
ggplot(df_density, aes(x = Status, y = meand100, color = Status)) +
  geom_boxplot() + theme_classic() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) + ylab("NN-100 Distance") +
  xlab("")  + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14))  + theme(legend.position = "none") + 
  scale_color_manual(values = c("#FF6A4C", "#487BE3")) +  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#Create mixed affect model (Supplementary Table 3) [Red asterisk means the significance was preserved]
segment.distances$Histology <- str_split_i(segment.distances$Subtype, "-", 2)
distance_df <- subset(segment.distances, select = c(Case, Status, Histology, meand100))
distance_df$Status <- as.factor(distance_df$Status)
distance_df$Case <- as.factor(distance_df$Case)
distance_df$Histology <- as.factor(distance_df$Histology)
summary(distance_df)
dist_eff <- lmer(meand100 ~ Status*Histology + (1|Case), data = distance_df)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = distance_df, adjust = "bonferonni")
























