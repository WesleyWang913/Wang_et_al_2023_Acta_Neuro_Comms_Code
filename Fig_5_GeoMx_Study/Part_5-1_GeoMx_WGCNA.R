#WGCNA over GeoMx Results
setwd("~/Downloads/Genomics Paper/Study 3_GeoMx")

#Load Libraries
library(WGCNA)
library(tidyr)
library(dplyr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(sjstats)
library(emmeans)

options(stringsAsFactors = FALSE)

#read in data
DSP_reads <- read.csv("GeoMx Expression Data.csv", row.names = 1, header = T) #Data can be found in GEO GSE231994
clin_data <- loadRDS("GeoMx_Clin.RDS") #A RDS copy of cleaned clinical data is in repository 
clin_data #Data is numeric scale, below is a key
#Key:
# 1 = PD; 2 = psPD
# 1= Control; 2 = Reactive; 3 = Recurrent
# 1= PD-control; 2 = PD-Reactive; 3 = PD-Recurrent 4= psPD-control; 5 =ps PD-Reactive; 6 = psPD-Recurrent

traitData <- clin_data
geneData <- DSP_genes

#Check if any samples have missing values
gsg = goodSamplesGenes(geneData, verbose = 3);
gsg$allOK

#Check samples clustering
datTraits <- traitData
datExpr <- geneData
sampleTree2 = hclust(dist(geneData), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Novel Enhancement Status Dendrogram and Heat Map")

#Optimization
#powers <- seq(0, 10, by= 0.2)
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#par(mfrow=c(2,1), mar= c(2,2,2,2))
#cex1=0.6;
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#     main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels=powers,cex=cex1,col="red");
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
#     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#     main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Network parameteres after validation
power=9.5; 
minModSize=30
enforceMMS=FALSE
net <- blockwiseModules(datExpr,power=power,deepSplit=4,minModuleSize=minModSize,
                        mergeCutHeight=0.07, TOMDenom="mean", #detectCutHeight=0.9999,
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,reassignThresh=0.05,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=12000)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
feature = as.data.frame(datTraits$Status);
names(feature) = "Status"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, feature, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(feature), sep="");
names(GSPvalue) = paste("p.GS.", names(feature), sep="")

#Gene Ontology
probes = names(datExpr)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
probe_2 <- AnnotationDbi::select(hs, keys = probes, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
GOenr = GOenrichmentAnalysis(moduleColors, probe_2$ENTREZID, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
#####GO Analysis Stored in Repository as tab.RDS################

tab <- readRDS("tab.RDS")
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols]

#Filter off modules that have non-sig enrichments
sigtab <- screenTab[screenTab$BonferoniP < 0.05,]
#Filter off CC
proctab <- sigtab[!sigtab$termOntology == "CC",]
proctab$module <- as.factor(proctab$module)
target_modules <- levels(proctab$module)
target_modules <- paste0("ME", target_modules)

#WGCNA Subtype comparison
ncol(MEs) #11 Color Modules
fMEs <- MEs[,colnames(MEs) %in% target_modules]
ncol(fMEs) #6 Color modules kept
me_data <- cbind(fMEs, DSP_reads[,5])
colnames(me_data)[7] <- "Subtype"

#Visualization
library(reshape)
me_long <- melt(me_data, id = "Subtype")       
me_long$Subtype <- factor(me_long$Subtype, levels = c("PD-Control", "psPD-Control", "PD-Hypercellular", "psPD-Hypercellular", "PD-Inflammatory", "psPD-Inflammatory"))

ggplot(me_long, aes(x = variable, y = value, color = Subtype)) +  # ggplot function
  geom_boxplot(outlier.shape = NA) + coord_flip() + stat_compare_means(method = "anova", label = "p.signif") +
  ylab("Module Eigengene Score") + xlab("") + theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4"))    

me_data$Subtype <- factor(me_data$Subtype, levels = c("PD-Control", "psPD-Control", "PD-Hypercellular", "psPD-Hypercellular", "PD-Inflammatory", "psPD-Inflammatory"))
my_comparisons <- list( c("PD-Control", "psPD-Control"), c("PD-Hypercellular", "psPD-Hypercellular"), c("PD-Inflammatory", "psPD-Inflammatory") )

library(tidyverse)
me_data$RoI <- rownames(me_data)
me_data$Case <- as.factor(str_split_i(me_data$RoI, "_",1))
me_data$Status <- as.factor(str_split_i(me_data$Subtype, "-",1))
me_data$Histology <- as.factor(str_split_i(me_data$Subtype, "-",2))
me_data


#For each color module, a boxplot comparison was made and post hoc mixed effect testing was done
#Brown (pd hiher in cellular only)
ggplot(me_data, aes(x= Subtype, y= MEbrown, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEbrown Score")+
  ylim(-0.2, 0.2) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEbrown ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")

#Turquoise (pd alsmost always higher excepy in inflammatory zones)
ggplot(me_data, aes(x= Subtype, y= MEturquoise, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEturquoise Score")+
  ylim(-0.2, 0.2)+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEturquoise ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")

#Magenta (PD hi)
ggplot(me_data, aes(x= Subtype, y= MEmagenta, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEmagenta Score")+
  ylim(-0.2, 0.2)+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEmagenta ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")

#Greenyellow (PD always higher)
ggplot(me_data, aes(x= Subtype, y= MEgreenyellow, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEgreenyellow Score")+
  ylim(-0.2, 0.24)+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEgreenyellow ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")

#Red (PD high in inlammatory)
ggplot(me_data, aes(x= Subtype, y= MEred, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEred Score")+
  ylim(-0.15, 0.25)+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEred ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")

#Pink (PD hi in control or inflammed)
ggplot(me_data, aes(x= Subtype, y= MEpink, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("MEpink Score")+
  ylim(-0.1, 0.23)+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

#Post hoc
dist_eff <- lmer(MEpink ~ Status*Histology + (1|Case), data = me_data)
anova(dist_eff)
effectsize::eta_squared(dist_eff, partial = T)
emmeans(dist_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")


