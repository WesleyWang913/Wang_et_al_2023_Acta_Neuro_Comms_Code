#Assessment of Immune Exhaustion Genes
library(ggpubr)
library(Boruta)
library(randomForest)

#To assess immune exhaustion lets look at PD-1 (CD274)/ CTLA4/ LAG3/ TIM3 (HAVCR2)
#CGGA Cases
setwd("Your Directory")

target <- c("HAVCR2", "LAG3", "CD274") #Designate the Genes

#Read in CGGA Gene and Clinical Data
Genes.CGGA <- read.delim("CGGA Expression Data.txt") #Download from CGGA website (http://www.cgga.org.cn/download.jsp) mRNAseq_693
Filtered.Genes <- Genes.CGGA[Genes.CGGA$Gene_Name %in% target,]
nrow(Filtered.Genes) #All Pulled
rownames(Filtered.Genes) <- Filtered.Genes$Gene
Filtered.Genes <- Filtered.Genes[,-1]
Clinical.CGGA <- read.delim("CGGA_Clinical_Data.txt") #Download from CGGA website

#Pull rGBM Cases
Clinical.CGGA <- Clinical.CGGA[!is.na(Clinical.CGGA$Histology),]
Recurrent.GBs <- Clinical.CGGA[Clinical.CGGA$Histology == "rGBM",]
nrow(Recurrent.GBs) #109 Cases
recurrent.genes <- Filtered.Genes[,Recurrent.GBs$CGGA_ID]
ncol(recurrent.genes) #109 Cases
head(recurrent.genes)
ttgs <- data.frame(t(recurrent.genes))

ttgs$ID <- rownames(ttgs)
colnames(ttgs) <- c("PD-1", "TIM-3", "LAG-3", "ID")

#Pull clinical data of clusters from WGCNA workflow
traitData$ID <- rownames(traitData)

#Merge with cluster identity
tf <- merge(ttgs, traitData, by = "ID")
tf$Cluster <- as.factor(tf$Cluster)
tf <- tf[order(tf$Cluster),]

#Plot CGGA comparisons
my_comparisons <- list( c("1", "2"), c("2", "3"), c("1", "3") )

ggplot(tf, aes(x= Cluster, y= tf$`PD-1`, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("PD-1 Expression")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

-ggplot(tf, aes(x= Cluster, y= tf$`TIM-3`, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("TIM-3 Expression")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggplot(tf, aes(x= Cluster, y= tf$`LAG-3`, color = Cluster))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test") + scale_color_manual(values = c("#20639B", "#1D7874", "#ED553B")) +
  theme_classic() + theme(legend.position = "none")  + ylab("LAG-3 Expression")+ theme(axis.text = element_text(size = 10), axis.title = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#Data imputation of GeoMx Genes
#Identify an appropriate housekeeping gene (low variance across samples that is consistent in both)
#CGGA
cases <- Clinical.CGGA[Clinical.CGGA$Histology == "rGBM",]$CGGA_ID
reads.cgga <- as.matrix(Genes.CGGA[,colnames(Genes.CGGA) %in% cases])
rownames(reads.cgga) <- Genes.CGGA$Gene_Name

#Remove CGGA gene with zero expression in sample
reads.cgga <- reads.cgga[rowSums(reads.cgga == 0) == 0,]

cgga.row.vars <- rowVars(reads.cgga)
quantile(cgga.row.vars, .20)
cgga.low.var.genes <- Genes.CGGA$Gene_Name[cgga.row.vars <= 8.03354]

#GeoMx
setwd("Your Directory")
Genes.GEO <- read.csv("GeoMx Expression Data.csv", row.names = 1, header = T) #Data can be found in GEO GSE231994
head(Genes.GEO)

reads.geo <- as.matrix(Genes.GEO)
rownames(reads.geo) <- Genes.GEO$Sample
reads.geo <- t(reads.geo)

geo.row.vars <- rowVars(reads.geo)
quantile(geo.row.vars, .20)
geo.low.var.genes <-rownames(reads.geo)[geo.row.vars <= 0.09282427]

low.var.genes <- cgga.low.var.genes[cgga.low.var.genes %in% geo.low.var.genes]

#Identify which gene has the lowest variace across both data sets
cgga.var <- data.frame("Gene" = rownames(reads.cgga[rownames(reads.cgga) %in% low.var.genes,]),
"CGGA_Variance" = rowVars(reads.cgga[rownames(reads.cgga) %in% low.var.genes,]))

cgga.var <- cgga.var[order(cgga.var$CGGA_Variance),]
cgga.var$Score_CGGA <- c(1:nrow(cgga.var))

geo.var <- data.frame("Gene" = rownames(reads.geo[rownames(reads.geo) %in% low.var.genes,]),
                       "GEO_Variance" = rowVars(reads.geo[rownames(reads.geo) %in% low.var.genes,]))

geo.var <- geo.var[order(geo.var$GEO_Variance),]
geo.var$Score_Geo<- c(1:nrow(geo.var))

all.var <- merge(cgga.var, geo.var, by = "Gene")

#Sum scores
all.var$Total_Score <- all.var$Score_CGGA + all.var$Score_Geo
all.var[order(all.var$Total_Score),] #Our winner is PRELID2

#Next step, normalize both datasets by their PRELID2 expression 
ref.cgga <- as.vector(reads.cgga[rownames(reads.cgga) == "PRELID2",])
div.cgga <- mapply('/', data.frame(reads.cgga), ref.cgga)
rownames(div.cgga) <- rownames(reads.cgga)
div.cgga[rownames(div.cgga) == "PRELID2",] #We have normalized to PRELID2

ref.geo <- as.vector(reads.geo[rownames(reads.geo) == "PRELID2",])
div.geo <- mapply('/', data.frame(reads.geo), ref.geo)
rownames(div.geo) <- rownames(reads.geo)
div.geo[rownames(div.geo) == "PRELID2",] #We have normalized to PRELID2

#Now we build borutas to predict the expression of our missing genes
library(randomForest)
library(Boruta)

#####CD274 (PD-1)
#Boruta
test <- data.frame(t(div.cgga))
test <- test[,colnames(test) %in% c("CD274", rownames(div.geo))]

pd1.boruta <- Boruta(CD274~., data= test, doTrace= 2, maxRuns=1000)
getConfirmedFormula(pd1.boruta)

#Split Data
pd1.rf <- randomForest(getConfirmedFormula(pd1.boruta), data = test, preProcess = c("center","scale"))
predict.geo <- data.frame(t(div.geo))
pd1.pred <- predict(pd1.rf, predict.geo)

#####HAVCR2 (TIM-3)
#Boruta
test <- data.frame(t(div.cgga))
test <- test[,colnames(test) %in% c("HAVCR2", rownames(div.geo))]

tim3.boruta <- Boruta(HAVCR2~., data= test, doTrace= 2, maxRuns=1000)
getConfirmedFormula(tim3.boruta)

#Split Data
tim3.rf <- randomForest(getConfirmedFormula(tim3.boruta), data = test, preProcess = c("center","scale"))
predict.geo <- data.frame(t(div.geo))
tim3.pred <- predict(tim3.rf, predict.geo)

#####LAG3 
#Boruta
test <- data.frame(t(div.cgga))
test <- test[,colnames(test) %in% c("LAG3", rownames(div.geo))]

LAG3.boruta <- Boruta(LAG3~., data= test, doTrace= 2, maxRuns=1000)
getConfirmedFormula(LAG3.boruta)

#Split Data
LAG3.rf <- randomForest(getConfirmedFormula(LAG3.boruta), data = test, preProcess = c("center","scale"))
predict.geo <- data.frame(t(div.geo))
LAG3.pred <- predict(LAG3.rf, predict.geo)


####Merge Predicted Immune Values
pred.exp <- data.frame(rbind(pd1.pred, tim3.pred, LAG3.pred))
rownames(pred.exp) <- c("PD1", "TIM3", "LAG3")


geo.pred.exp <- data.frame("Subtype" = Genes.GEO$Subtype, t(pred.exp))
me_data <- geo.pred.exp

me_data$Subtype <- factor(me_data$Subtype, levels = c("PD-Control", "psPD-Control", "PD-Hypercellular", "psPD-Hypercellular", "PD-Inflammatory", "psPD-Inflammatory"))
my_comparisons <- list( c("PD-Control", "psPD-Control"), c("PD-Hypercellular", "psPD-Hypercellular"), c("PD-Inflammatory", "psPD-Inflammatory") )

ggplot(me_data, aes(x= Subtype, y= PD1, color = Subtype))+
  geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
  scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("PD1 Imputed Expression")+
  ylim(-0.2, 0.2) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")
ave(file="test5.svg", p, width=2.6041666666666665, height=3.125)
 
ggplot(me_data, aes(x= Subtype, y= LAG3, color = Subtype))+
   geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
   scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
   theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("LAG3 Imputed Expression")+
   ylim(-0.2, 0.2) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

-ggplot(me_data, aes(x= Subtype, y= TIM3, color = Subtype))+
   geom_boxplot() + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", step.increase = 0) +
   scale_color_manual(values = c("#005923","#34A56F", "#1A1ABF", "#7F7FFF", "#c31d3e", "#FF83A4")) +
   theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))  + ylab("TIM3 Imputed Expression")+
   ylim(-0.2, 0.2) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + xlab("")

 #Post hoc-testing 
 #Create mixed affect model
 library(lme4)
 library(lmerTest)
 library(sjstats)
 library(emmeans)
 library(tidyverse)
 
 me_data$Status <-  as.factor(str_split_i(me_data$Subtype, "-", 1))
 me_data$Histology <-  as.factor(str_split_i(me_data$Subtype, "-", 2))
 me_data$Case <- as.factor(str_split_i(rownames(me_data), "_", 1))
 summary(me_data) #Make sure groupings are as factors
 
 #PD1
 immune_eff <- lmer(PD1 ~ Status*Histology + (1|Case), data = me_data)
 anova(immune_eff)
 effectsize::eta_squared(immune_eff, partial = T)
 emmeans(immune_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")
 
 #TIM3
 immune_eff <- lmer(TIM3 ~ Status*Histology + (1|Case), data = me_data)
 anova(immune_eff)
 effectsize::eta_squared(immune_eff, partial = T)
 emmeans(immune_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")
 
 #LAG3
 immune_eff <- lmer(LAG3 ~ Status*Histology + (1|Case), data = me_data)
 anova(immune_eff)
 effectsize::eta_squared(immune_eff, partial = T)
 emmeans(immune_eff, list(pairwise ~ Status | Histology), data = me_data, adjust = "bonferonni")
 
 
 
 
 