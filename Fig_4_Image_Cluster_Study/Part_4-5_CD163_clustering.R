#Evaluation of cluster distribution in relation to PD ans psPD status
setwd("Your Directory")

##############Please note UMAP is computationally heavy and a finalized version of processed data is present in the repository below#################
seg_reads <- readRDS("cd163_Segmentation_masterframe.RDS")
nrow(seg_reads)

#Split off PD Cases
PD.table <- seg_reads[seg_reads$Case < 30,]
nrow(PD.table) 

#Split off psPD cases
psPD.table <- seg_reads[seg_reads$Case > 30,]
nrow(psPD.table)

#Sample data to remove class inbalance 
set.seed(888)
PD.sample <- PD.table[sample(nrow(PD.table), 12000),]
psPD.sample <- psPD.table[sample(nrow(psPD.table), 12000),]
cell.table <- rbind(PD.sample, psPD.sample)
base.table <- cell.table #Store all the selected reads for use in the image pulling

#UMAP (We should only keep morphology features [s and m features])
cell.table <- cell.table[,!grepl("\\.h\\.", colnames(cell.table))] #remove haralick features (pixel texture)
cell.table <- cell.table[,!grepl("\\.b\\.", colnames(cell.table))] #remove intensity features
cell.table <- cell.table[,!grepl("cx", colnames(cell.table))] #remove x positions
cell.table <- cell.table[,!grepl("cy", colnames(cell.table))] #remove y positions
cell.table <- subset(cell.table, select = -c(X, Histology, Subtype, Case, Image)) #Remove miscellaneous

#####Note I run some analysis and found the color layers are just redundant, so I'm just gonna take the red layers
cell.table <- cell.table[,!grepl("G\\.", colnames(cell.table))] 
cell.table <- cell.table[,!grepl("B\\.", colnames(cell.table))] 
colnames(cell.table)
ncol(cell.table)

###############UMAP Analysis#############
library(umap)
library(tidyverse)
library(ggpubr)

#Create new data tables for UMAP analysis
rownames(cell.table) <- NULL 
cell.table <- cell.table %>%
  drop_na() %>%
  mutate(ID=row_number())
cell.meta <- cell.table %>% #store information of cells
  select(Status, ID)

#UMAP Calculation
set.seed(999)
umap_fit <- cell.table %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>%
  umap()

#Store UMAP positions
umap_df <- umap_fit$layout %>%
  as.data.frame() %>%
  mutate(ID=row_number())%>%
  inner_join(cell.meta, by="ID")

colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2") #the code was glitching for some reason

#Plot UMAP (Note there are massive outliers so I focused on the primary clusters area)
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = Status))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2")+
  theme(legend.position="none")+
  xlim(-15,15)+
  ylim(-10,10)

base.table <- base.table[umap_df$UMAP1 > -15 & umap_df$UMAP1 < 15,]
nrow(base.table)
umap_df1 <- umap_df[umap_df$UMAP1 > -15 & umap_df$UMAP1 < 15,]
nrow(umap_df1)

base.table <- base.table[umap_df1$UMAP2 > -10 & umap_df1$UMAP2 < 10,]
nrow(base.table)
umap_df2 <- umap_df1[umap_df1$UMAP2 > -10 & umap_df1$UMAP2 < 10,]
nrow(umap_df2)

umap_df <- umap_df2

umap_df <- umap_df[umap_df$UMAP1 > -15 & umap_df$UMAP1 < 15,]
umap_df <- umap_df[umap_df$UMAP2 > -10 & umap_df$UMAP2 < 10,]

##############Clustering UMAP#################
library(dbscan)
library(fpc)
library(RColorBrewer)

#This plot helps you identify what you eps should be but tbh it only helps get a starting point and then you need to test multiple eps
dbscan::kNNdistplot(umap_df[,1:2], k =  200)
abline(h = 0.75, lty = 2)

# Compute DBSCAN using fpc package
set.seed(123)
db <- fpc::dbscan(umap_df[,1:2], eps = 0.9, MinPts = 200)
umap_df$Cluster <- as.factor(db$cluster)

##########################Start Here For Reanalysis####################
umap_df <- readRDS("CD163_umap_df.RDS")
cell.table <- readRDS("CD163_sampled_cells.RDS")

#Add Cluster Data (Note: make 0 cluster great and rename NA)
umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = Cluster))+
  geom_point(size=0.011)+
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_classic()+
  scale_color_manual(values = c("grey", brewer.pal(n = 8, name = "Dark2")))+
  theme(legend.position = "none") 

#Evaluate the distribution of clusters
filtered_umap_df <- umap_df[!umap_df$Cluster == 0,]
filtered_umap_df$Cluster <- droplevels(filtered_umap_df$Cluster)

#Show distribution using chart
filtered_umap_df %>%
  count(Cluster, Status) %>%
  group_by(Cluster) %>%
  mutate(n = n/sum(n) * 100) %>%
  ggplot() + aes(Cluster, n, fill = Status, label = paste0(round(n, 0), "%")) + 
  geom_col() +
  geom_text(position=position_stack(0.5), size = 3.5, color = "white", angle = 90)+
  theme_classic()+
  ylab("Percent of Cluster")+
  theme(legend.position = "none")+ scale_fill_manual(values = c("#FF6A4C", "#487BE3"))+
  geom_hline(yintercept = 50, linetype = "dashed", size = 0.8)

library(gplots)
dt <- table(filtered_umap_df$Cluster, filtered_umap_df$Status)
dt

chisq <- chisq.test(dt)
chisq
chisq$expected

#Evaluate what the morphology is ########Please note this portion requires source images (by request from authors)############
ttable <- cbind(base.table, umap_df$Cluster) #Please note due to outliers, this needs to be pulled from above processing code
ttable[ttable$`umap_df$Cluster` == 1,] 

#Read in images and spit out a cropped version of the segmented cells
library(EBImage)
setwd("Your Directory")  
test <- ttable[ttable$`umap_df$Cluster` == X,] #Modify X to be the cluster you are evaluating
nrow(test)
for (i in 1:400){
  tryCatch({
    print(i)
    im_name <- gsub("feature.*$", "overlay_h.png", test$Image[i])
    img <- readImage(im_name)
    x <- round(test$R.0.m.cx[i], digits = 0)
    y <- round(test$R.0.m.cy[i], digits = 0)
    im2 <- img[(x-20):(x+20), (y-20):(y+20), 1:3]
    writeImage(im2, files = paste0("our Directory", i, ".png"))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}





