library(limma)
library(edgeR)
# treatment_group <- "Drsp"
# control_group <- "Cont"
limmaall <- function( experiment_title, control_group, treatment_group) {

  #list("ST_MCF7_ZR751","GSE172478_RPE","GSE151321_HK","GSE113659_PAE","GSE160497_MCF7","GSE28009_EA","GSE84992_SM","GSE70822_PM")
  if(experiment_title=="ST_MCF7_ZR751"){
    data <- read.csv("data_files/ST_MCF7_ZR751_Filtered.csv", header=TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("MCF7_EV" ,    "MCF7_EVA"  ,  "MCF7_EVAS"  , "MCF7_OV"   ,  "MCF7_OVA" , "MCF7_OVA" , "MCF7_OVAS" ,"MCF7_OVAS" , "ZR751_EV"  ,  "ZR751_EVA" ,  "ZR751_OV"  ,  "ZR751_OVA"  ))
    #design model matrix
    model_matrix <- model.matrix(~0+group)
    
    colnames(model_matrix) <-sort(c("MCF7_EV", "MCF7_OV", "MCF7_EVA","MCF7_EVAS",
                                    "MCF7_OVA", 
                                    "MCF7_OVAS",
                                    "ZR751_EV", "ZR751_EVA",
                                    "ZR751_OV"  ,  "ZR751_OVA"))
    
    
    #dge list formation
    dge <- DGEList(data)
    dge <- calcNormFactors(dge)
    #voom transformation
    v <- voom(dge, model_matrix)
    fit <- lmFit(v, model_matrix)
    
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    
    contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = colnames(coef(fit)))
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_ST " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
  }
  if(experiment_title=="GSE113659_PAE"){
    #drsp_RawCountData_Filtered <- read.csv("files/Seniye_MCF7_Filtered.csv", header=TRUE)
    
    data <- read.csv("data_files/GSE113659_Rawdata_Filtered.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("NEDD9siRNA_ALDO","NEDD9siRNA_ALDO", "NEDD9siRNA_ALDO", 
                      "NEDD9siRNA_Control", "NEDD9siRNA_Control", "NEDD9siRNA_Control",
                      "Control",  "Control" , "Control" , "Control" , "Control" , "Control",
                      "Aldo", "Aldo", "Aldo", "Aldo", "Aldo", "Aldo"
    ))
    
    #design model matrix
    model_matrix<- model.matrix(~0+group)

    
    colnames(model_matrix) <- sort(c("NEDD9siRNA_ALDO", "NEDD9siRNA_Control",
                                     "Control",  "Aldo"))
    
    #dge list formation
    dge <- DGEList(data)
    dge <- calcNormFactors(dge)
    #voom transformation
    v <- voom(dge, model_matrix)
    fit <- lmFit(v, model_matrix)
    
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    
    contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = colnames(coef(fit)))
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE113659 " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
  }

  if(experiment_title=="GSE151321_HK"){
    
    data <- read.csv("data_files/GSE151321_RNAseq_Filtered.csv", header=TRUE)
    
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("Control", "Aldo","Fine.Spiro","Fine.Spiro","Aldo_Fine","Aldo_Spiro")) #spiro ya fine.spiro gibi davran ki replikalı olsun
    #design model matrix
    model_matrix <- model.matrix(~0+group)
    
    colnames(model_matrix) <- sort(c( "Control", "Aldo","Fine.Spiro","Aldo_Fine","Aldo_Spiro" ))
    
    
    #dge list formation
    dge <- DGEList(data)
    dge <- calcNormFactors(dge)
    #voom transformation
    v <- voom(dge, model_matrix)
    fit <- lmFit(v, model_matrix)
    
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    
    contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = colnames(coef(fit)))
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE151321" ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
  }
  if(experiment_title=="GSE172478_RPE"){
    data <- read.csv("data_files/GSE172478_Filtered.csv", header=TRUE)
    rownames(data) <- data[,1]
    data <- data[,-1]
    group <- factor(c("Etoh","Etoh","Etoh","Methanol","Methanol","Methanol","Methanol",
                      "Aldo","Aldo","Aldo","Cortisol","Cortisol",
                      "Cortisol",  "Cortisol",  "Cortisol_RU486",  "Cortisol_RU486","Cortisol_RU486","Cortisol_RU486"
    )
    )
    #design model matrix
    model_matrix <- model.matrix(~0+group)
    
    colnames(model_matrix) <- sort(c("Etoh","Methanol",
                                     "Aldo",
                                     "Cortisol",  "Cortisol_RU486"
    )
    )
    
    
    #dge list formation
    dge <- DGEList(data)
    dge <- calcNormFactors(dge)
    #voom transformation
    v <- voom(dge, model_matrix)
    fit <- lmFit(v, model_matrix)
    
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    
    contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = colnames(coef(fit)))
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE172478 " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
  }
  if(experiment_title=="GSE160497_MCF7"){
    data=read.csv("data_files/GSE160497_Array_filtered.csv")[-c(1:11),]
    rownames(data)=data$Symbol
    data=data[,-1]
    sampleInfo=readRDS("data_files/sampleInfo.RDS")
    design <- model.matrix(~0+sampleInfo$group)
    
    ## the column names are a bit ugly, so we will rename
    colnames(design) <- sort(c("Dox","minusDox","Dox_Aldo","minusDox_Aldo","Dox_E2","Dox_E2_Aldo","Dox_RA","Dox_RA_Aldo"))
   
    fit <- lmFit(data, design)
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    contrasts <- makeContrasts(contrasts = contrast_name, levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE160497_Array " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
    
    
  }
  if(experiment_title=="GSE28009_EA"){
    data=read.csv("data_files/GSE28009_Array_filtered.csv")
    rownames(data)=data$Symbol
    data=data[,-1]
    sampleInfo=readRDS("data_files/sampleInfoGSE28009.RDS")
    design <- model.matrix(~0+sampleInfo$group)
    
    ## the column names are a bit ugly, so we will rename
    colnames(design) <- sort(c("Control","Aldo"))
    
    fit <- lmFit(data, design)
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    contrasts <- makeContrasts(contrasts = contrast_name, levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE28009_Array " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
    
    
  }
  if(experiment_title=="GSE84992_SM"){
    data=read.csv("data_files/GSE84992_Array_filtered.csv")
    rownames(data)=data$Symbol
    data=data[,-1]
    sampleInfo=readRDS("data_files/sampleInfoGSE84992.RDS")
    design <- model.matrix(~0+sampleInfo$group)
    
    
    colnames(design) <- sort(c("Aldo_48h",
                               "Spironolactone_48h",
                               "Control_48h",
                               "Eplerenone_48h",
                               "Mifepristone_48h",  
                               "Prednisolone_48h",
                               "Aldo_24h",
                               "Aldo_Spiro_24h", 
                               "Aldo_Eplerenone_24h",
                               "Aldo_Mifepristone_24h"
    ))
    
    fit <- lmFit(data, design)
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    contrasts <- makeContrasts(contrasts = contrast_name, levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    result_df=na.omit(result_df)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE84992_Array " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
    
    
  }
  if(experiment_title=="GSE70822_PM"){
    data=read.csv("data_files/GSE70822_Array_filtered.csv")
    rownames(data)=data$Symbol
    data=data[,-1]
    sampleInfo=readRDS("data_files/sampleInfoGSE70822.RDS")
    design <- model.matrix(~0+sampleInfo$group)
    
    
    colnames(design) <- sort(c("Aldo","Spironolactone","Control"  ))
    
    fit <- lmFit(data, design)
    contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
    contrasts <- makeContrasts(contrasts = contrast_name, levels=design)
    
    fit2 <- contrasts.fit(fit, contrasts)
    fit2 <- eBayes(fit2)
    
    
    result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
    contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE70822_Array " ,sep = "")
    colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
    result_df <- round(result_df,10)
    
    
  }
  return(result_df)
}


# if(experiment_title=="GSE151321"){
# 
#   data <- read.csv("data_files/GSE151321_RNAseq_Filtered.csv", header=TRUE)
#   
#   rownames(data) <- data[,1]
#   data <- data[,-1]
#   group <- factor(c("Control", "Aldo","Fine.Spiro","Fine.Spiro","Aldo_Fine.Aldo_Spiro","Aldo_Fine.Aldo_Spiro"))
#   #design model matrix
#   model_matrix <- model.matrix(~0+group)
#   
#   colnames(model_matrix) <- sort(c( "Control", "Aldo","Fine.Spiro","Aldo_Fine.Aldo_Spiro" ))
#   
#   
#   #dge list formation
#   dge <- DGEList(data)
#   dge <- calcNormFactors(dge)
#   #voom transformation
#   v <- voom(dge, model_matrix)
#   fit <- lmFit(v, model_matrix)
# 
#   contrast_name <- paste(treatment_group,  "-", control_group, sep = "")
#   
#   contrast_matrix <- makeContrasts(contrasts = contrast_name, levels = colnames(coef(fit)))
#   
#   fit2 <- contrasts.fit(fit, contrast_matrix)
#   fit2 <- eBayes(fit2)
#   
#   result_df <- topTable(fit2, sort.by ="logFC", n = Inf)
#   contrast_namechange <- paste(treatment_group,  "_vs_", control_group,"_GSE151321 " ,sep = "")
#   colnames(result_df) <- c( paste( contrast_namechange , c("logFC", "AveExpr", "t value", "p value", "p adj value", "B value")))
#   result_df <- round(result_df,10)
# }