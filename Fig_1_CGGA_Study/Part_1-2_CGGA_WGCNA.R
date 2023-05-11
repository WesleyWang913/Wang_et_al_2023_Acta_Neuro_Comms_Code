#CGGA WGCNA Study (Please run Part 1-1 first to have all libraries)
setwd("Your_Directory")

# Load the WGCNA package
library(WGCNA);
options(stringsAsFactors = FALSE);

#Read in cluster data
GBM_clin <- readRDS("GBM_clin.RDS") #A RDS copy of cleaned clinical data is in repository
traitData <- data.frame("Cluster" = GBM_clin$Cluster)
rownames(traitData) <- GBM_clin$CGGA_ID

#Read in gene
geneData <- read.delim("CGGA Expression Data.csv") #Download from CGGA website (http://www.cgga.org.cn/download.jsp) mRNAseq_693
rownames(geneData) <- geneData$Gene_Name
geneData <- geneData[,colnames(geneData) %in% rownames(traitData)] #109 rGBM Samples
datExpr0 = as.data.frame(t(geneData));

Samples = rownames(datExpr0);
traitRows = match(Samples, rownames(traitData));
datTraits = traitData
collectGarbage();


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Recurrent CGGA GBM Dendrogram and Heat Map")

datExpr <- datExpr0

#Below is the optimization of the network
# Choose a set of soft-thresholding powers
#powers <- seq(4,15,by=0.5)
# Call the network topology analysis function
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#     main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
#     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#     main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Network analysis (note: network has been constructed if you want to skip)
power=10.5; #first tried: 7.6
minModSize=30
enforceMMS=FALSE
net <- blockwiseModules(datExpr,power=power,deepSplit=4,minModuleSize=minModSize,
                        mergeCutHeight=0.07, TOMDenom="mean", #detectCutHeight=0.9999,
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,reassignThresh=0.05,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=12000)

# open a graphics window
sizeGrWindow(12, 9)
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

#######Start here if you want to skip the network contruction###########
load("networkConstruction-auto_9.RData") #Loads MEs, moduleLabels, moduleColors, geneTree

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#Gene Ontology
probes = names(datExpr)

#Changes names of gene symbols (note the next code chunk will load in the results)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
probe_2 <- AnnotationDbi::select(hs, keys = probes, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
GOenr = GOenrichmentAnalysis(moduleColors, probe_2$ENTREZID, organism = "human", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
#write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

#Read in GO results
tab<- read.csv("GOEnrichmentTable.csv")
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols]

#Filter off modules that have non-sig enrichments
sigtab <- screenTab[screenTab$BonferoniP < 0.05,]
#Filter off CC
proctab <- sigtab[!sigtab$termOntology == "CC",]
proctab$module <- as.factor(proctab$module)
target_modules <- levels(proctab$module)
target_modules <- paste0("ME", target_modules)

#Comparing enrichment of a module t o cluster
ncol(MEs) #85 Color Modules
fMEs <- MEs[,colnames(MEs) %in% target_modules]
ncol(fMEs) #25 Color modules kept
me_data <- cbind(fMEs, traitData)

me_long <- melt(me_data, id = "Cluster")       
me_long$Cluster <- as.factor(me_long$Cluster)

#Plot WGCNA Results
ggplot(me_long, aes(x = variable, y = value, color = Cluster)) +  # ggplot function
  geom_boxplot(outlier.shape = NA) + coord_flip() + stat_compare_means(method = "anova", label = "p.signif") +
  ylab("Module Eigengene Score") + xlab("") + theme_classic()   + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) + theme(legend.position = "none")

#Plot specific groups
my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )
me_data$Cluster <- as.factor(me_data$Cluster)

ggplot(me_data, aes(x= Cluster, y= MEnavajowhite1, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEnavajowhite1 Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
 
ggplot(me_data, aes(x= Cluster, y= MElightyellow, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MElightyellow Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

 ggplot(me_data, aes(x= Cluster, y= MEplum2, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEplum2 Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

 ggplot(me_data, aes(x= Cluster, y= MElightgreen, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MElightgreen Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

 ggplot(me_data, aes(x= Cluster, y= MEyellow, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEyellow Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

 ggplot(me_data, aes(x= Cluster, y= MEcyan, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEcyan Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

 ggplot(me_data, aes(x= Cluster, y= MEpink, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEpink Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

  ggplot(me_data, aes(x= Cluster, y= MEhoneydew1, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("MEhoneydew1 Score")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

  
  #Supplementary Figure 5: WGCNA module and immune proliferation gene overlap.
  ####Asessing genes in a module with GO term overlap
  #library(pheatmap)
  #library(org.Hs.eg.db)
  #library(GO.db)
  #library(ggvenn)
  #Pull Lymphocyte Proliferation Genes
  #results <- AnnotationDbi::select(org.Hs.eg.db, keys=c("GO:0046651"), columns = c('SYMBOL'), keytype = "GOALL")
  #gene_symbols <- unique(results$SYMBOL)
  #length(gene_symbols) #296 Genes
  #lymph_genes <- gene_symbols
  #Pull Macrophage Proliferation Genes
  #results <- AnnotationDbi::select(org.Hs.eg.db, keys=c("GO:0061517"), columns = c('SYMBOL'), keytype = "GOALL")
  #gene_symbols <- unique(results$SYMBOL)
  #length(gene_symbols) #12 Genes
  #macro_genes <- gene_symbols
  #Import genes in the WGCNA modules (Pink contains th)
  #library(readxl)
  #pink_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "pink"))
  #nrow(pink_genes) #163 Genes
  #honeydew1_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "honeydew1"))
  #nrow(honeydew1_genes) #110 Genes
  #All these modules were testested but plum2 and light yellow were not interesting so we omitted them
  #plum2_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "plum2"))
  #nrow(plum2_genes) #46 Genes
  #navajowhite1_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "navajowhite1"))
  #nrow(navajowhite1_genes) #752 Genes
  #lightyellow_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "lightyellow"))
  #nrow(lightyellow_genes) #298 Genes
  #lightgreen_genes <- data.frame(read_xlsx("CGGA_WGCNA_Modules.xlsx", sheet = "lightgreen"))
  #nrow(lightgreen_genes) #163 Genes
  #Inspect overlap
  #Cell Cycle Modules
  #x <- list(
  #  'Lymphocytes' = lymph_genes,
  #  'Pink Module' = pink_genes$Gene, 
  #  'Honeydew1 Module' = honeydew1_genes$Gene, 
  #  'Macrophages' = macro_genes)
  #ggvenn(x, 
  #  fill_color = c("#0073C2FF", "pink", "honeydew1", "#CD534CFF"),
  #  stroke_size = 0.5, set_name_size = 4)
  #Immune Modules
  #y <- list(
  #  'Lymphocytes' = lymph_genes,
  #  'Lightgreen Module' = lightyellow_genes$Gene,
  #  'Navajowhite1 Module' = navajowhite1_genes$Gene,
  #  'Macrophages' = macro_genes)
  #ggvenn(
  # y, fill_color = c("#0073C2FF", "lightgreen", "navajowhite1", "#CD534CFF"),
  #  stroke_size = 0.5, set_name_size = 4)
  
 
  
  
  
  



