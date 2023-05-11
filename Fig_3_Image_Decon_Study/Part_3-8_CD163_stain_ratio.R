#Analysis of computed features from CD163 Stain Ratio
library(ggpubr)
setwd("Your Directory")

####Below analyses have been completed and outputted ratio files for visualization is in repository#####
#Read in H Segments
features_list <- list.files(path = "features")
h_list <- features_list[grepl("_h", features_list)] #pull h layer
d_list <- features_list[grepl("_d", features_list)] 
length(h_list)
length(d_list)

#Double check the indices match up
l <- sample(1:260, 1)
d_list[l]
h_list[l]

##Read in the total D and H pixels for each image and then calculate the ratio
PD.cases <- c(4, 7, 12, 23)
psPD.cases <- c(36,53, 50, 56)

ratio_frame <- NULL
d_frame <- NULL
for (n in 1:length(d_list)){
  print(n)
  d_features <- read.csv(paste0("~/Downloads/Genomics Paper/Study 4_Image Analysis/CD163/features/", d_list[n]))
  d_features <- d_features[!d_features$R.0.m.majoraxis == 0,] #There was a random pixel artifact that needs to be filtered off
  d_features$Image <- str_extract(d_list[n], "[^_]+")
  d_features$Case <- sub('.+-(.+)', '\\1', d_features$Image)
  d_features$Image <- d_list[n]
  d_features$Histology <- sapply(strsplit(d_list[n], "_"), function(x) x[2])
  if (d_features$Case[1] %in% PD.cases){
    d_features$Status = "PD"
  } else { 
    d_features$Status = "psPD"
  }
  d_features$Subtype <- paste0(d_features$Status,"-",d_features$Histology)
  d_area <- sum(d_features$B.0.s.area)
  d_objects <- nrow(d_features)
  d_frame <- rbind(d_frame, d_features)
  
  h_features <- read.csv(paste0("~/Downloads/Genomics Paper/Study 4_Image Analysis/CD163/features/", h_list[n]))
  h_area <- sum(h_features$B.0.s.area)
  h_objects <- (nrow(h_features))
  
  ration <- d_area/h_area
  o_ration <- d_objects/h_objects
  
  r_store <- data.frame("Case" = str_extract(d_list[n], "[^_]+"),
                        "Images" = paste(d_list[n], h_list[n]),
                        "DAB" = d_area,
                        "Hem" = h_area,
                        "Ratio" = ration,
                        "D_Objects" = d_objects,
                        "H_Objects" = h_objects,
                        "Object_Ratio" = o_ration,
                        "Histology" = strsplit(d_list[n], "_")[[1]][2],
                        "Status" = d_features$Status[1],
                        "Subtype" = d_features$Subtype[1])
  ratio_frame <- rbind(ratio_frame, r_store)
}

####################Start Here######################
ratio_frame <- read.csv("cd163_ratio_frame.csv")

df_ratio <- ratio_frame
df_ratio$Subtype <- as.factor(df_ratio$Subtype)
levels(df_ratio$Subtype) <- c("PD-Control", "PD-Inflammatory", "PD-Hypercellular", "psPD-Control", "psPD-Inflammatory", "psPD-Hypercellular")
df_ratio$Subtype <- factor(df_ratio$Subtype, levels = c("PD-Control","psPD-Control",
                                                        "PD-Hypercellular", "psPD-Hypercellular",
                                                        "PD-Inflammatory","psPD-Inflammatory"))

df_ratio <- df_ratio[df_ratio$Ratio < 2,] #Outliers were confirmed to be image artifact and need to be removed

my_comparisons <- list( c("PD-Control", "psPD-Control"), c("PD-Hypercellular", "psPD-Hypercellular"), c("PD-Inflammatory", "psPD-Inflammatory") )

ggplot(df_ratio, aes(x = Subtype, y = Ratio, color = Subtype)) +
  geom_boxplot() + stat_compare_means(label = "p.signif", method = "t.test", step.increase = 0, comparisons = my_comparisons) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("CD163 D:H Ratio") +
  xlab("") + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


my_comparisons <- list( c("PD", "psPD") )

ggplot(df_ratio, aes(x = Status, y = Ratio, color = Status)) +
  geom_boxplot() + stat_compare_means(label = "p.signif", method = "t.test", step.increase = 0, comparisons = my_comparisons) +
  scale_color_manual(values = c("#FF6A4C", "#487BE3")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("CD163 D:H Ratio") +
  xlab("")  + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


#Post hoc-testing 
#Create mixed affect model
library(lme4)
library(lmerTest)
library(sjstats)
library(emmeans)

summary(df_ratio) #Make sure groupings are as factors
df_ratio$Histology <- as.factor(df_ratio$Histology)
df_ratio$Status <- as.factor(df_ratio$Status)
df_ratio$Case <- as.factor(df_ratio$Case)

stain_eff <- lmer(Ratio ~ Status*Histology + (1|Case), data = df_ratio)
anova(stain_eff)
effectsize::eta_squared(stain_eff, partial = T)
emmeans(stain_eff, list(pairwise ~ Status | Histology), data = df_ratio, adjust = "bonferonni")













