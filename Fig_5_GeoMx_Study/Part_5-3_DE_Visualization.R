#Volcano Plot Visulization
setwd("Your Directory")

#Load Libraries
library(EnhancedVolcano)

#DE datasets are generated within the GeoMx software using data provided in GEO GSE231994
dat_react <- read_xlsx("Reactive DE Analysis.xlsx")
dat_recur <- read_xlsx("Recurrent DE Analysis.xlsx")
dat_control <- read_xlsx("Control DE Analysisl.xlsx")

#Set the active comparison group (comment our the one not used)
res <- dat_recur
res <- dat_react
res <- dat_recur

#Volcano Plot
 EnhancedVolcano(res,
                lab = res$`Target name`,
                x = 'Log2',
                y = 'Pvalue',
                title = NULL,
                FCcutoff = 0.5,
                pCutoff = .05,
                drawConnectors = T,
                widthConnectors = 0.1,
                colConnectors = 'black', 
                labSize = 5, 
                pointSize = 1,
                legendPosition = 'none',
                subtitle = "",
                ylim = c(0, -log10(10e-5)),
                xlim = c(-2,2))#+ scale_x_reverse()

