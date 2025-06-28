options(repos = BiocManager::repositories())

#library(shinysky)

library(shinydisconnect)

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(rintrojs)
library(fresh)
library(shinyWidgets) # for useShinyDashboard()
library(DT)
library(heatmaply)
library(plotly)
library(limma)
library(edgeR)
library(DESeq2)
library(visNetwork)
library(shinyBS)
library(shinycssloaders)
library(EnrichmentBrowser)
library(stringr)
library(markdown)
library(msigdbr)
library(hypeR)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(gprofiler2)
library(RColorBrewer)
library(EnhancedVolcano)
#library(scales)

keggpat <-  readRDS("keggpathwaysandgenes.rds")
# reading files in global
## -----------READING DATASETS -----------

ST.rawdata <- read.csv("data_files/ST_Rawdata.csv", header=TRUE)
ST_allfiltered <- read.csv("data_files/ST_MCF7_ZR751_Filtered.csv", header=TRUE)
ST_onlyMCF.filtered <- read.csv("data_files/ST_only_MCF7_Filtered_sep.csv", header=TRUE)
ST_onlyZR75.filtered <- read.csv("data_files/ST_only_ZR751_Filtered_sep.csv", header=TRUE)
GSE172478.rawdata <- read.csv("data_files/GSE172478_entrez_hgnc_rawdata.csv", header=TRUE)


GSE172478.filtered <- read.csv("data_files/GSE172478_Filtered.csv", header=TRUE)
GSE151321.rawdata <- read.csv("data_files/GSE151321_Rawdata_withentrez.csv", header=TRUE)[,c(1,8,2:7)]
GSE151321.filtered=read.csv("data_files/GSE151321_RNAseq_Filtered.csv", header=TRUE)

GSE113659.rawdata <- read.csv("data_files/GSE113659_Rawdata_withentrez.csv", header=TRUE)
GSE113659.filtered  <- read.csv("data_files/GSE113659_Rawdata_Filtered.csv", header=TRUE)
###Arrays
GSE160497_Array=read.csv("data_files/GSE160497_Array_filtered.csv")
GSE28009_Array=read.csv("data_files/GSE28009_Array_filtered.csv")
GSE84992_Array=read.csv("data_files/GSE84992_Array_filtered.csv")
GSE84992_Array <- GSE84992_Array[complete.cases(GSE84992_Array), ] # remove na values

GSE70822_Array=read.csv("data_files/GSE70822_Array_filtered.csv")
GSE70822_Array <- GSE70822_Array[complete.cases(GSE70822_Array), ] # remove na values

gene.list <- read.csv("data_files/genename_ensembleID_63682genes.csv",header = T)
rownames(gene.list) <- gene.list$genename
#GSE151321 <- read.csv("data_files/GSE151321_kidney_Filtered.csv", header=TRUE)
