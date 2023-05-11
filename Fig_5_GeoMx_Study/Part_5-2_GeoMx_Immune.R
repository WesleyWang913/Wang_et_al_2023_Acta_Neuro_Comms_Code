#Immune Analysis of GeoMx Data
setwd("Your Directory")

#Read in Libraries
library(readxl)
library(ggpubr)
library(stringr)

#Load data
immune_decon_final <- read.csv("immune_decon_final.csv") #Immune decon was performed within GeoMx (A copy is found in the repository)
immune_decon_final

#Fix Subtype label
immune_decon_final$Subtype <- clin_data$Subtype
immune_decon_final$Case_ID <- as.factor(immune_decon_final$Case_ID)
immune_decon_final$Status <- as.factor(immune_decon_final$Status)
immune_decon_final$Morphology <- as.factor(immune_decon_final$Morphology )

#Lets make the long figure
#Dissociate the clin data to have only the splitting variable and immune reads
imf <- immune_decon_final[,5:ncol(immune_decon_final)]
imf$Subtype <- factor(imf$Subtype, levels = c("PD-Control", "psPD-Control", "PD-Hypercellular", "psPD-Hypercellular", "PD-Inflammatory", "psPD-Inflammatory"))
imf_long <- melt(imf, id = "Subtype")                      
imf_long

levels(imf_long$variable) <-  c("Macrophages", "Mast Cells", "Naive B Cells", "Memory B Cells", "Plasma Cells", "Naive CD4 T Cells",
                                "Memory CD4 T Cells", "Naive CD8 T Cells", "Memory CD8 T Cells", "NK Cells", "Plasmacytoid Dendritics Cells",
                                "Myeloid Dendritic Cells", "Classical Monocytes" ,"Non-classical Monocytes", "Neutrophils",
                                "Tregs", "Endothelial Cells", "Fibroblasts")

ggplot(imf_long, aes(x = variable, y = value, color = Subtype)) +  # ggplot function
  geom_boxplot(outlier.shape = NA) + coord_flip() + stat_compare_means(method = "anova", label = "p.signif") +
  ylab("Enrichment Score") + xlab("") + theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

my_comparisons <- list( c("PD-Control", "psPD-Control"), c("PD-Hypercellular", "psPD-Hypercellular"), c("PD-Inflammatory", "psPD-Inflammatory") )

#macrophages (more in PD reactive)
ggplot(imf, aes(x= Subtype, y= macrophages, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("Macrophage Enrichment") + 
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  ylim(0,0.8) + xlab("") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

#Fibroblasts (PD all high)
#ggplot(imf, aes(x= Subtype, y= fibroblasts, color = Subtype))+
#  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") +
#  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("Fibroblast Enrichment") + 
#  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) 

#neutrophils (more in the psPD recurrent)
ggplot(imf, aes(x= Subtype, y= neutrophils, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("Neutrophil Enrichment") + 
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4"))+
  ylim(0,0.11)+ xlab("") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

#monocytes.nc.1 ),pre om PD reactive)
ggplot(imf, aes(x= Subtype, y= monocytes.NC.I, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("NC Monocyte Enrichment") + 
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4"))+
  ylim(0,0.09)+ xlab("") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

#naive.cd8 (more in PD control, but then more in psPD reactive)
ggplot(imf, aes(x= Subtype, y= T.CD8.naive, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("Naive CD8 T Cell Enrichment") + 
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4"))+
  ylim(0,0.9)+ xlab("") +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

#Post hoc-testing 
#Create mixed affect model
library(lme4)
library(lmerTest)
library(sjstats)
library(emmeans)

summary(immune_decon_final) #Make sure groupings are as factors

immune_eff <- lmer(neutrophils ~ Status*Morphology + (1|Case_ID), data = immune_decon_final)
anova(immune_eff)
effectsize::eta_squared(immune_eff, partial = T)
emmeans(immune_eff, list(pairwise ~ Status | Morphology), data = immune_decon_final, adjust = "bonferonni")

immune_eff <- lmer(T.CD8.naive ~ Status*Morphology + (1|Case_ID), data = immune_decon_final)
anova(immune_eff)
effectsize::eta_squared(immune_eff, partial = T)
emmeans(immune_eff, list(pairwise ~ Status | Morphology), data = immune_decon_final, adjust = "bonferonni")

immune_eff <- lmer(macrophages ~ Status*Morphology + (1|Case_ID), data = immune_decon_final)
anova(immune_eff)
effectsize::eta_squared(immune_eff, partial = T)
emmeans(immune_eff, list(pairwise ~ Status | Morphology), data = immune_decon_final, adjust = "bonferonni")

immune_eff <- lmer(monocytes.NC.I ~ Status*Morphology + (1|Case_ID), data = immune_decon_final)
anova(immune_eff)
effectsize::eta_squared(immune_eff, partial = T)
emmeans(immune_eff, list(pairwise ~ Status | Morphology), data = immune_decon_final, adjust = "bonferonni")

