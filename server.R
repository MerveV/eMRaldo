
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
library(shinyalert)



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



# setwd("~/Desktop/eMRald")

server <- function(input,output,session){
  ldat <- list("ST_MCF7_ZR751","GSE172478_RPE","GSE151321_HK","GSE113659_PAE","GSE160497_MCF7","GSE28009_EA","GSE84992_SM","GSE70822_PM")
  onlyrnaseq=ldat[1:4]
  

  # reactiveValues
  reac <- reactiveValues()
 # source("data_read.R")
  
  
  ##### Analysis of RNAseq 
  
  observeEvent(input$analyze,
               {
                
                 
                 reac$exp1 <- input$experiment_title_1
                 reac$exp2 <- input$experiment_title_2
                 reac$exp3 <- input$experiment_title_3
                 reac$exp4 <- input$experiment_title_4
                 reac$metod <- input$metapvalue
                 reac$control_group <- input$control_group
                 reac$control_group_2<- input$control_group_2
                 reac$control_group_3<- input$control_group_3
                 reac$control_group_4<- input$control_group_4
                 reac$treatment_group<- input$treatment_group
                 reac$treatment_group_2<- input$treatment_group_2
                 reac$treatment_group_3<- input$treatment_group_3
                 reac$treatment_group_4<- input$treatment_group_4
                 reac$uploadoption <- input$uploadoption
                 
               })
  #####
  ST.Group.all <- list( "MCF7_EV" ,    "MCF7_EVA"  ,  "MCF7_EVAS"  , "MCF7_OV"   ,  "MCF7_OVA" , "MCF7_OVAS" , "ZR751_EV"  ,  "ZR751_EVA" ,  "ZR751_OV"  ,  "ZR751_OVA"  ) ####------------------Seniye

  GSE172478_groups=list("Etoh","Methanol","Aldo","Cortisol","Cortisol_RU486")
  array_group=list("Dox","minusDox","Dox_Aldo","minusDox_Aldo","Dox_E2","Dox_E2_Aldo","Dox_RA","Dox_RA_Aldo")
  GSE28009_Array_groups=list("Control","Aldo")
  GSE84992_Array_groups=list("Control_48h","Aldo_48h","Spironolactone_48h","Eplerenone_48h","Mifepristone_48h",  "Prednisolone_48h","Aldo_24h","Aldo_Spiro_24h", "Aldo_Eplerenone_24h","Aldo_Mifepristone_24h")
  
  GSE70822_Array_groups=list("Control","Aldo","Spironolactone")
  GSE151321_groups=list( "Control", "Aldo","Fine.Spiro","Aldo_Fine.Aldo_Spiro")
  GSE151321_groups2=list( "Control", "Aldo","Aldo_Fine","Aldo_Spiro")
  GSE113659_groups=list("Control","Aldo",  "NEDD9siRNA_Control","NEDD9siRNA_ALDO")
  #GSE151321_kidneyl=list("Control","Aldo","Finerenone","Spironolactone","Aldo_Spiro","Aldo_Fine")
  
  
  rawdatafunction <- function(exp){
    raw_data <- switch(exp,
                       "ST_MCF7_ZR751"=ST.rawdata,
                     
                       "GSE172478_RPE"=GSE172478.rawdata,
                       "GSE160497_MCF7"=GSE160497_Array[-c(1:11),], # remove unrelevant rows
                       "GSE28009_EA"=GSE28009_Array,
                       "GSE84992_SM"=GSE84992_Array,
                       "GSE70822_PM"=GSE70822_Array,
                      # "GSE151321"=GSE151321.rawdata,
                       "GSE151321_HK"=GSE151321.rawdata,
                       "GSE113659_PAE"=GSE113659.rawdata
    )
  }
  ########## Select Raw count data
  switc <- function(exp){
    group_names <- switch(exp,
                          "ST_MCF7_ZR751"=ST.Group.all,
                      
                          "GSE172478_RPE"=GSE172478_groups,
                          "GSE160497_MCF7"=array_group,
                          "GSE28009_EA"=GSE28009_Array_groups,
                          "GSE84992_SM"=GSE84992_Array_groups,
                          "GSE70822_PM"=GSE70822_Array_groups,
                         # "GSE151321"=GSE151321_groups,
                          "GSE151321_HK"=GSE151321_groups2,
                          "GSE113659_PAE"=GSE113659_groups
    )
  }
  
  # control and treatment dropdown menus for the first comparison
  output$select_control_output <- renderUI({
    
    group_names <- switc(input$experiment_title_1)
    
    selectInput(inputId = "control_group",
                label = "Control:",
                choices = group_names,
                selected = group_names[1])
    
    
  })
  output$select_treatment_output <- renderUI({
    
    group_names <- switc(input$experiment_title_1)
    
    selectInput(inputId = "treatment_group",
                label = "Treatment:",
                choices = group_names,
                selected = group_names[2])
  })
  
  # control and treatment dropdown menus for the second comparison
  output$select_control_output_2 <- renderUI({
    group_names <- switc(input$experiment_title_2)
    selectInput(inputId = "control_group_2",
                label = "Control:",
                choices = group_names,
                selected = group_names[1])
  })
  
  output$select_treatment_output_2 <- renderUI({
    group_names <- switc(input$experiment_title_2)
    
    selectInput(inputId = "treatment_group_2",
                label = "Treatment:",
                choices = group_names,
                selected = group_names[2])
  })
  
  # control and treatment dropdown menus for the third comparison
  output$select_control_output_3 <- renderUI({
    group_names <- switc(input$experiment_title_3)
    
    selectInput(inputId = "control_group_3",
                label = "Control:",
                choices = group_names,
                selected = group_names[1])
  })
  
  output$select_treatment_output_3 <- renderUI({
    group_names <- switc(input$experiment_title_3)
    
    selectInput(inputId = "treatment_group_3",
                label = "Treatment:",
                choices = group_names,
                selected = group_names[2])
  })
  
  # control and treatment dropdown menus for the 4th comparison
  output$select_control_output_4 <- renderUI({
    group_names <- switc(input$experiment_title_4)
    
    selectInput(inputId = "control_group_4",
                label = "Control:",
                choices = group_names,
                selected = group_names[1])
  })
  
  output$select_treatment_output_4 <- renderUI({
    group_names <- switc(input$experiment_title_4)
    selectInput(inputId = "treatment_group_4",
                label = "Treatment:",
                choices = group_names,
                selected = group_names[2])
  })
  # exp selection is shorten with numbers

  conditions <- reactive({
    req(input$analyze)
    if (reac$exp2== "Select"&&reac$exp3 == "Select"&& reac$exp4 != "Select"){ 
      conditions <- "14"
    }
    else if (reac$exp2 == "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
      conditions <- "13"
    }    # 1 ve 2
    else if(reac$exp3 == "Select"&& reac$exp4 == "Select"&& reac$exp2 != "Select"){
      conditions <- "12"
    }
    else if(reac$exp2 != "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
      conditions <- "123"
    } # 1,3,4
    else if(reac$exp3 != "Select"&& reac$exp2 == "Select"&& reac$exp4 != "Select"){
      conditions <- "134"
    }
    else if(reac$exp2 != "Select"&& reac$exp3 == "Select"&& reac$exp4 != "Select"){
      conditions <- "124"
    }
    else if(reac$exp2 != "Select"&& reac$exp3 != "Select"&& reac$exp4 != "Select"){
      conditions <- "1234"
    }
    else if(reac$exp2 == "Select"&& reac$exp3 == "Select"&& reac$exp4 == "Select"){
      conditions <- "1"
    }
  })

  # titles are formed here
  titles <- reactive({
    req(input$analyze)
    if (reac$exp2 == "Select"&& reac$exp3 == "Select"&& reac$exp4 != "Select"){ 
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group ," and ",
                   reac$exp4," ",reac$treatment_group_4,  "_vs_", reac$control_group_4 ,  sep = "")
    }
    else if (reac$exp2 == "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group  ," and ",reac$exp3," ",reac$treatment_group_3,  "_vs_", reac$control_group_3,  sep = "")
      
    }    # 1 ve 2
    else if(reac$exp3 == "Select"&& reac$exp4 == "Select"&& reac$exp2 != "Select"){
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group  ," and ",reac$exp2," ",reac$treatment_group_2,  "_vs_", reac$control_group_2,  sep = "")
    }
    else if(reac$exp2 != "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
      titles= paste(reac$exp1," ",reac$treatment_group,  "_vs_", 
                    reac$control_group  ," and ",reac$exp2," ",reac$treatment_group_2,  "_vs_", reac$control_group_2," and "
                    ,reac$exp3," ",reac$treatment_group_3,  "_vs_", reac$control_group_3,  sep = "")
    } # 1,3,4
    else if(reac$exp3 != "Select"&& reac$exp2 == "Select"&& reac$exp4 != "Select"){
      titles= paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group  ," and ",reac$exp3," ",
                    reac$treatment_group_3,  "_vs_", reac$control_group_3," and ",reac$exp4," ",reac$treatment_group_4,  "_vs_", 
                    reac$control_group_4,  sep = "")
    }
    else if(reac$exp2 != "Select"&& reac$exp3 == "Select"&& reac$exp4 != "Select"){
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group  ," and ",reac$exp2," ",reac$treatment_group_2,  "_vs_", reac$control_group_2," and ",reac$exp4," ",reac$treatment_group_4,  "_vs_", reac$control_group_4,  sep = "")
    }
    else if(reac$exp2 != "Select"&& reac$exp3 != "Select"&& reac$exp4 != "Select"){
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group  ," and ",
                   reac$exp2," ",reac$treatment_group_2,  "_vs_", reac$control_group_2," and ",
                   reac$exp3," ",reac$treatment_group_3,  "_vs_", reac$control_group_3," and ",
                   reac$exp4," ",reac$treatment_group_4,  "_vs_", reac$control_group_4,
                   sep = "")
    }
    else if(reac$exp2 == "Select"&& reac$exp3 == "Select"&& reac$exp4 == "Select"){
      titles=paste(reac$exp1," ",reac$treatment_group,  "_vs_", reac$control_group , sep = "")
    }
    
  }  ) 
  source("limmaall.R")
  
  func_meta <- function(vec,method){
    
    if (method=="fisher"){res <- as.numeric(poolr::fisher(vec)[["p"]])}
    else if (method=="stouffer"){res <- as.numeric(poolr::stouffer(vec)[["p"]])}
    else if (method=="invchisq"){res <- as.numeric(poolr::invchisq(vec)[["p"]])}
    else if (method=="binomtest"){res <- as.numeric(poolr::binomtest(vec)[["p"]])}
    else if (method=="bonferroni"){res <- as.numeric(poolr::bonferroni(vec)[["p"]])}
    else{res <- as.numeric(poolr::tippett(vec)[["p"]])}
    
    return(res)
  }
  
  metaanalysis=function(data,methodmeta,columnames){
    # data=reac$limma
    # methodmeta = reac$metod 
    result <- data.frame(data)
    result3 <- result %>% dplyr::select(matches(c("p.value","p value"),ignore.case = TRUE)) # select columns containing p-values
    logfcres <- result %>% dplyr::select(matches("logFC",ignore.case = TRUE)) # select columns containing logFC
 
    empt <- c()
    upc <- c()
    downc <- c()
    metod <- methodmeta
    for (i in 1:nrow(result3)) {
      k <- as.vector(unlist(unname(result3[i,]))) # make all values in same row as vector
      empt <- c(empt,func_meta(k,metod)) # calculate meta-analysis
      l <- as.vector(unlist(unname(logfcres[i,]))) # logFC values of same row
      if(any(k<=0.05)==TRUE){
        sig <- l[which(k<=0.05)]
        up <- length(which(as.vector(sig) > 0)) # number of upregulated genes
        down <- length(which(as.vector(sig) < 0))# number of downregulated genes
      }else {
        up <- 0
        down <- 0}
      upc <- c(upc,up)
      downc <- c(downc,down)
    }
    matchExpression <- paste(columnames, collapse = "|")
    result2 <- result #%>% dplyr::select(matches(matchExpression))
    
    result <- cbind(result[,1],round(result2[,-1],3),round(empt,3), upc,downc)
    colnames(result)[1] <- "hgnc_symbol"
    print(head(result))
    colnames(result)[ncol(result)-2] <- "Meta–analysis of p–values"
    colnames(result)[ncol(result)-1] <- c("Sig upregulated")
    colnames(result)[ncol(result)] <- c("Sig downregulated")
    return(result)
  }
  
  observeEvent(input$analyze, {
    if (conditions()=="1") shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"))
    if (conditions()=="12")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"))
    if (conditions()=="14")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    
    
    if (conditions()=="13")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_3 !=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="123")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="124")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="134")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="1234")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_3!=input$treatment_group_4, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
      withProgress(message = 'Analyzing the selected data with Limma-Voom and merging them...', style = "notification", value = 0.1, {
        Sys.sleep(0.25)
        
        if (conditions()=="1"){
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom...', detail = "First data")
          l<- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.6)
          reac$limma  <- cbind(rownames(  l ), l  )
          rownames(reac$limma ) <- 1:dim(reac$limma  )[1]
          colnames(reac$limma )[1] <- "hgnc_symbol"
        }
        else if (conditions()=="12"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma2 <- limmaall( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
          reac$limma <-base::merge(limma1,limma2, by="row.names")
          colnames(reac$limma)[1] <- "hgnc_symbol"
        }
        else if (conditions()=="14"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma4 <- limmaall( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
          reac$limma <-base::merge(limma1,limma4, by="row.names")
          colnames(reac$limma)[1] <- "hgnc_symbol"
        }
        else if (conditions()=="13"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma3 <- limmaall( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
          reac$limma <-base::merge(limma1,limma3, by="row.names")
          colnames(reac$limma)[1] <- "hgnc_symbol"
        }
        else if (conditions()=="123"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma2 <- limmaall( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
          incProgress(0.7, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "Third data")
          limma3 <- limmaall( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
          l12  <- merge(limma1,limma2, by="row.names")
          colnames(  l12)[1] <- "hgnc_symbol"
          limma3 <-cbind(rownames(limma3),limma3) 
          colnames(limma3)[1] <- "hgnc_symbol"
          reac$limma <- merge(l12 ,limma3, by="hgnc_symbol")
        }
        else if (conditions()=="124"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma2 <- limmaall( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
          incProgress(0.7, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "Third data")
          limma4 <- limmaall( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
          l12  <- merge(limma1,limma2, by="row.names")
          colnames(  l12)[1] <- "hgnc_symbol"
          limma4 <-cbind(rownames(limma4),limma4) 
          colnames(limma4)[1] <- "hgnc_symbol"
          reac$limma <- merge(l12 ,limma4, by="hgnc_symbol")
        }
        else if (conditions()=="134"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma2 <- limmaall( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
          incProgress(0.7, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "Third data")
          limma4 <- limmaall( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
          l12  <- merge(limma1,limma2, by="row.names")
          colnames(  l12)[1] <- "hgnc_symbol"
          limma4 <-cbind(rownames(limma4),limma4) 
          colnames(limma4)[1] <- "hgnc_symbol"
          reac$limma <- merge(l12 ,limma4, by="hgnc_symbol")
        }
        else if (conditions()=="1234"){ 
          incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
          limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
          incProgress(0.4, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
          limma2 <- limmaall( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
          incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "Third data")
          limma3 <- limmaall( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
          incProgress(0.7, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "Fourth data")
          limma4 <- limmaall( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
          l12<- merge(limma1,limma2, by="row.names")
          colnames(l12 )[1] <- "hgnc_symbol"
          l3 <-cbind(rownames(limma3),limma3) 
          incProgress(0.8)
          colnames(l3)[1] <- "hgnc_symbol"
          l123  <- merge(  l12 ,l3, by="hgnc_symbol")
          l4<-cbind(rownames(limma4),limma4) 
          colnames(l4)[1] <- "hgnc_symbol"
          reac$limma  <- merge(  l123 ,l4, by="hgnc_symbol")
        }
        
        shiny::setProgress(1)
      })
    
   
    
    
  })

  ##########                               ####################                               ###############
  # result part 
  output$results_ui_output <- renderUI({
    if (conditions()=="1") shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"))
    if (conditions()=="12")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"))
    if (conditions()=="14")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    
    
    if (conditions()=="13")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_3 !=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="123")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="124")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="134")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="1234")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_3!=input$treatment_group_4, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
   req( input$analyze)
 
      tabsetPanel(type = "tabs", id = "tabsetss",
                  tabPanel(title="Analysis",
           
                         # setSliderColor(c(rep("#00838f",600)), c(1:600)),
                           
                          panel(class="myclass1",fluidRow( div(style="display: inline-block;vertical-align:top; width: 20px;",HTML("<br>")),
                                         div(style="display: inline-block;vertical-align:top; width: 400px;", selectInput('in4', 'Display Pathway-Specific Genes', c("Display all genes", names(keggpat)[-c(2,53)]), selectize=TRUE)),
                                         
                                          
                                         div(style="display: inline-block;vertical-align:top;width: 300px;",fileInput("genefile",  tippy::tippy(text = "Upload Genes",tooltip = "Gene names should be in the first column",arrow = TRUE, placement = "auto") ,multiple = FALSE,accept = c("text/csv","text/comma-separated-values,text/plain",".csv"))),
                                         div(style="display: inline-block;vertical-align:top;width: 150px;",  awesomeRadio("sep", "Separator",choices = c(Comma = ",",Semicolon = ";",Tab = "\t"),selected = ",", status = "warning")),
                                         div(style="display: inline-block;vertical-align:top;width: 100px;",  prettyCheckbox("header", "Header", FALSE,status = "warning",shape="round") ),
                                         div(style="display: inline-block;vertical-align:middle;width: 200px;",column(12, fluidRow(HTML("<br>")),fluidRow(materialSwitch("usegenes","Use Uploaded Genes",value = FALSE, status = "warning"))))
                          ))   ,
                          
                         panel( fluidRow(div(style = 'overflow-x: scroll',DT::DTOutput("tablefirst")%>% shinycssloaders::withSpinner(type = getOption("spinner.type", 6),
                                                                                                                                      color = getOption("spinner.color", default = "#00838f"),
                                                                                                                                      size = getOption("spinner.size", default = 1) ))), 
                                downloadButton(outputId = "d1_4", label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;") ,
                                
                                fluidRow(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="#17414E"><b>Select the rows to further visualize</b></font>'),
                                                  bsButton("bsb355", label="", icon=icon("question"), style="info", size="small")),
                                         shinyBS::bsPopover("bsb355", 'When you select the genes, new table containing selected genes is created.', trigger = "focus")),
                                fluidRow(DTOutput("subsettable"),
                                         downloadButton(outputId = "subset_d1_4", label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                         )
                           ),# Analysis tab ends
                  tabPanel(title="Scatter Plot",
                    br(),
                   
                       panel(class="myclass1", fluidRow(column(3,selectInput("select_xaxis",label = "Select x-axis",choices = updatednames())),
                                       column(3,selectInput("select_yaxis","Select y-axis",choices = updatednames())),
                                       column(3,sliderInput("dotsize","Point size",min=1,max=20,value = 2))),
                              
                              fluidRow(column(4,materialSwitch( inputId = "scatter_material",label = "Use selected genes",value = FALSE, status = "warning")),
                                       column(2,actionBttn("scatter","Draw" ,icon =icon("bolt"),color = "warning")))),
                       
                  fluidRow( plotlyOutput("scatterplot",height = '800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))),
                  fluidRow(br())
                    
                  ),
                  tabPanel(title="Heatmap",br(),
                           panel(class="myclass1",
                                               fluidRow(column(width = 2,htmltools::h4('Row dendrogram')),column(width = 2,offset = 2,htmltools::h4('Column dendrogram'))),
                                               fluidRow(shiny::column(width=2,shiny::selectizeInput("distFun_row", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                                        shiny::column(width=2,shiny::selectizeInput("hclustFun_row", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                                        shiny::column(width=2,shiny::selectizeInput("distFun_col", "Distance method", c(Euclidean="euclidean",Maximum='maximum',Manhattan='manhattan',Canberra='canberra',Binary='binary',Minkowski='minkowski'),selected = 'euclidean')),
                                                        shiny::column(width=2,shiny::selectizeInput("hclustFun_col", "Clustering linkage", c(Complete= "complete",Single= "single",Average= "average",Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",Ward.D= "ward.D",Ward.D2= "ward.D2"),selected = 'complete')),
                                                        shiny::column(width=4,shiny::sliderInput("r", "Number of Clusters", min = 1, max = 15, value = 2)))),
                                               
                                               fluidRow( htmltools::tags$a(id = 'downloadData', class = paste("btn btn-default shiny-download-link",'mybutton'), href = "", target = "_blank", download = NA, shiny::icon("clone"), 'Download Heatmap as HTML'),
                                                          plotlyOutput("heatmapplot",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1)))),# Heatmap tab ends
                  ###############################                 ###############################"                ###############################
                  tabPanel(title="Functional Profile",
                           
                           fluidRow( column(4,panel(heading="Interactive Gost plot", status = "primary",
                                                    actionBttn("gostbut","Run",style = "jelly",icon =icon("bolt"), color = "warning"),plotlyOutput("gost1",height = "600px"),
                                                    style = "height: 700px;")),
                                     column(8, panel( heading="DisGeNET ",style = "height: 700px;",status = "primary",
                                          div(style="display: inline-block;vertical-align:top; width: 200px;",   
                                              column(3, dropdown(inputId="disgenetdropdown" ,label ="Parameters", status = "success", width = "800px",circle = FALSE,margin = "20px",
                                                             icon = icon("gear"),# style = "stretch", status = "default", width = "800px",
                                                             # animate = animateOptions(enter = animations$fading_entrances$fadeInDown,exit = animations$fading_exits$fadeOutUp),
                                                             tooltip = tooltipOptions(title = "Click to select!"), 
                                                             
                                                             fluidRow(
                                                                      column(4,sliderInput("pvaluecutof", label = h4("p-value cutoff"), min = 0, max = 1, value = 0.05)),
                                                                      column(4,sliderInput("qvaluecutof", label = h4("q-value cutoff"), min = 0, max = 1, value = 0.2)),
                                                                      column(4,awesomeRadio("circular", label = h4("Use circular layout"), choices = c("TRUE"=TRUE,"FALSE"=FALSE),selected="TRUE", status = "warning", inline = T))
                                                                      ),
                                                             
                                                             fluidRow(column(6,sliderInput("showCategory", label = h4("A number of terms to display"), min = 1, max = 20, value = 3) ),
                                                                      column(6,sliderInput("minGSSize1", label = h4("Min. size of genes for testing"), min = 2, max = 50, value = 2) )
                                                                      )
                                          ))),      
                                          div(style="display: inline-block;vertical-align:left; width: 600px;", 
                                              column(3,actionBttn("actionfordisgenet","Draw",icon =icon("bolt"),style = "jelly", color = "warning")),
                                                        column(6,uiOutput("LogFC1UI"))),
                                       
                                          br(),fluidRow( tabsetPanel(
                                                                tabPanel(title = "Gene-Concept Network",
                                                                         plotOutput("network1",height = "520px") ),
                                                                tabPanel(title="Heatplot",
                                                                         
                                                                         plotOutput("heatplot1",height = "520px") )
                                                                
                                          ))
                                     ))
                           ),
                           tabsetPanel(
                                  tabPanel(title="GO Enrichment Analysis ",
                                           br(), panel(class="myclass1",
                                           fluidRow(#column(width=5,awesomeRadio("selectont1", label = h4("Select subontology"),  choices = c("ALL"="ALL","Molecular Function" = "MF", "Biological Process" = "BP", "Cellular Component"="CC"),inline = TRUE,selected = "ALL",status = "danger")),  
                                             column(width=4,sliderInput("selectlevel1", label = h4("Min. size of genes for testing"), min = 1, max = 10, value = 2)),
                                             column(width=2,actionBttn("got","Create Table",icon =icon("bolt"),style = "jelly", color = "warning")))),
                                           br(),br(),
                                           fluidRow(  dataTableOutput("GOtable1")),br(),br(),br(),br() ),
                                  
                                  tabPanel(title="Visualization of Enrichment",
                                           panel(class="myclass1", 
                                                 fluidRow(column(3,sliderInput("slid", label = ("Number of displayed terms"), min = 2, max = 60, value = 30)),
                                                    column(2,uiOutput("goLogFC1UI")), 
                                                    column(2,actionBttn("dotplotbttn","Run",icon =icon("play"),style = "jelly", color = "warning")),
                                                    column(4,awesomeRadio("layot","Select layout",c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl','lgl'),  inline = TRUE,selected = "circle", status = "warning"))
                                                   )),  
                                           tabBox(
                                             title = "",
                                             id = "tabset1", height = "250px",
                                             tabPanel("Dot plot", plotOutput("dotplt",height = "800px")%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))),
                                             tabPanel("Enrichment Map", plotOutput("emapplot",height = "800px",width = "1400px")%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))),
                                             tabPanel("Gene-Concept Network",plotOutput("cnplt",width = "1500px",height = "1000px")%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1)) )
                                           )
                                          )
                           )
                           
                  ) , tabPanel(title="Transcription Factors",  
                              panel(class="myclass1",
                               fluidRow(column(width=2,awesomeRadio("selectlist", "Select geneset list",choices = c("All"="all", "NR targets"="nrtargets"),inline = TRUE,  status = "warning",checkbox = TRUE)),
                                        conditionalPanel(condition = "input.selectlist == 'nrtargets'",column(width=2,awesomeRadio("selectbackground", "Select Background",choices = c("All genes"="allgenes", "NR target genes"="nrtargetsgenes"),inline = F, status = "warning", checkbox = TRUE))),
                                        column(width=2,actionBttn("DoRothEAbutton","RUN",icon =icon("play"),style = "jelly", color = "warning" )),
                                               column(width=5, HTML(paste0(" <p>The number of all transcription factors is 1333.  <p>"," <p> The number of only nuclear receptors (NR) is 47. <p>"))
                                               )
                                        )),
                               
                               conditionalPanel(condition = "input.DoRothEAbutton != 0", 
                                                tabsetPanel(
                                                       tabPanel(title="Table ",
                                                                fluidRow( DT::dataTableOutput("DoRothEAtable")),
                                                                 fluidRow(   downloadButton(outputId = "DoRothEAdownload", label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;")),
                                                                br(),br(),
                                                                tabsetPanel(
                                                                  tabPanel(title="Heatmap ",
                                                                fluidRow(   column(width = 3,plotOutput("dotplotfordorothea"),br(),br(),br(),
                                                                                   wellPanel(style="font-size: 16px; line-height: 24px;",HTML(paste0(" <p> In the annotation part of the Heatmap, values show the confidence scores classifying regulons based on their quality (ranging from A (highest confidence) to E (lowest confidence)). <p>",
                                                                                                         "<p> A: Interactions that receive support from all lines of evidence are regarded as highly reliable. <p>",
                                                                                                         "<p> B, C, D: For curated and/or ChIP-seq interactions with varying degrees of additional evidence.<p>",
                                                                                                         "<p> E: Interactions that are solely supported by computational predictions. <p>",
                                                                                                         "<p> To better estimate TF activities, it is recommended to select regulons from the confidence levels A, B and C.<p>")))) ,
                                                                            column(width = 9,fluidRow(column(2,fluidRow(h4(" Use dendrogram")),
                                                                                                             fluidRow(   materialSwitch("usedendrog","For samples",value = TRUE, status = "warning")),
                                                                                                             fluidRow(   materialSwitch("usedendrogrow","For genes",value = TRUE, status = "warning"))
                                                                                                             
                                                                                                             ),
                                                                                                      column(6, dropdown(inputId="dropdowndorothea" ,label ="More Parameters for Heatmap", status = "success", width = "800px",circle = FALSE,margin = "20px",
                                                                                                                         icon = icon("folder-open",lib = "glyphicon"),# style = "stretch", status = "default", width = "800px",
                                                                                                                         
                                                                                                                         tooltip = tooltipOptions(title = "Click to select!"), 
                                                                                                                         fluidRow(awesomeRadio("selecteddist", "Distance method",choices = c("dist (default)"="dist","pearson", "spearman" , "kendall"),inline = TRUE,  status = "warning", checkbox = TRUE)),
                                                                                                                         
                                                                                                                         
                                                                                                                         fluidRow(column(3,sliderInput("col_text_Ang", label = h4("Column text angle"), min = 0, max = 90, value = 45)),
                                                                                                                                  column(3,sliderInput("fontsize_rowc", label = h4("Row font size"), min = 1, max = 20, value = 10)),
                                                                                                                                  column(3,sliderInput("fontsize_colc", label = h4("Column font size"), min = 1, max = 20, value = 10))
                                                                                                                                  
                                                                                                                                  
                                                                                                                         )
                                                                                                                         
                                                                                                      )
                                                                                                             )),
                                                                                   
                                                                                  fluidRow( plotlyOutput("heatmapplotdorothea",height='800px')),
                                                                                   fluidRow(column(4,offset=1,verbatimTextOutput("textTFsheatmap")))
                                                                                   ) ),br(),br(),br(),br() ) , 
                                                                tabPanel(title="MsigDB ",br(),
                                                                         DT::DTOutput("MSigtableTF"),br(),
                                                                         fluidRow(plotlyOutput("msigdbTFheatmap",height='600px')), fluidRow(column(4,offset=1,verbatimTextOutput("textmisgdbTF"))),
                                                                         
                                                                         br()) )),
                                                       
                                                       tabPanel(title="Network",
                                                                sidebarLayout(sidebarPanel(width = 3,
                                                                                        
                                                                                           fluidRow(align = 'center',column(10,awesomeRadio("nodesize","Node size based on", choices = c("FDR"="fdr", "Number of overlapped genes"="numbgen"),selected = "numbgen", status = "warning",inline = TRUE))),
                                                                                           fluidRow(align = 'center',column(10, 
                                                                                                                            awesomeRadio("nodeclr","Node color based on", choices = c("Unique Data"="udta", "Difference"="diffr"),selected = "udta", status = "warning",inline = TRUE)
                                                                                                                            
                                                                                                                            
                                                                                                                            )),
                                                                                           fluidRow(align = 'center',column(10, conditionalPanel(condition = "input.nodeclr=='udta'",
                                                                                                                                                 uiOutput("nodecolorui")),
                                                                                                                            conditionalPanel(condition = "input.nodeclr=='diffr'",
                                                                                                                                             uiOutput("differencecolor"))
                                                                                                                          
                                                                                                                            
                                                                                                                            
                                                                                           )),
                                                                                           fluidRow(align = 'center',column(10,awesomeRadio("fdr_tf","TF selection based on",choices = c("FDR","TF"),selected = "TF", status = "warning",inline = TRUE))),
                                                                                            fluidRow(align = 'center',column(10, conditionalPanel(condition = "input.fdr_tf=='FDR'",
                                                                                                                             sliderInput("fdrcutof", label = ("Select FDR cutoff"), min = 0, max = 1, value = 0.05,step = 0.01)     ),
                                                                                                            conditionalPanel(condition = "input.fdr_tf=='TF'",sliderInput("tfnumber", label = ("Number of TFs"), min = 1, max = 30, value = 5,step = 1) ))),
                                                                                           fluidRow(column(10,awesomeRadio("cutnetwork","Cut network based on ",choices = c("Jaccard index"="jacc","Number of shared genes"="shrdgene","Percentage"="percnt","Edge's logFC"="edgelogfc","Significancy"="sign"),selected = "jacc", status = "warning",inline = FALSE))),
                                                                                           fluidRow(align = 'center',column(10, conditionalPanel(condition = "input.cutnetwork=='jacc'",
                                                                                                                                                 sliderInput("jaccardindex","Cut network based on Jaccard index", min = 0, max = 1, value = 0,step = 0.1)     ),
                                                                                                                            conditionalPanel(condition = "input.cutnetwork=='shrdgene'",
                                                                                                                                             sliderInput("targetnumber", label = ("Number of shared genes"), min = 0, max = 30, value = 0,step = 1) ),
                                                                                                                            conditionalPanel(condition = "input.cutnetwork=='percnt'",
                                                                                                                                             sliderInput("percntgene", label = ("Percentage of shared genes"), min = 0, max = 100, value = 10,step = 1) ),
                                                                                                                            conditionalPanel(condition = "input.cutnetwork=='edgelogfc'",
                                                                                                                                             sliderInput("edgelogfcslider", label = ("Absolute mean logFC of shared genes"), min = 0, max = 5, value = 2,step = 0.1) ),
                                                                                                                            conditionalPanel(condition = "input.cutnetwork=='sign'",
                                                                                                                                             awesomeRadio("signifcnt", label = ("Wilcoxon test"),choices = c("Significant","Not Significant"),selected = "Significant", status = "warning",inline = TRUE ) )
                                                                                                                            
                                                                                                                            
                                                                                                                            )),
                                                                                        
                                                                                           
                                                                                            fluidRow(align = 'center',column(8,offset=1,actionBttn("networkDorothe","Create network",icon =icon("play"),style = "jelly", color = "warning" ))),
                                                                                           br(),
                                                                                          
                                                                                           fluidRow(actionButton("visexp111",  "Explanation" ,icon = icon("info"), class = "btn-success"),
                                                                                                    
                                                                                                    bsModal("visexp112221", "", "visexp111", size = "large",includeMarkdown("mdfiles/visnetworkexplain.md"))  )
                                                                                           
                                                                                           
                                                                                            ),
                                                                              mainPanel(fluidRow(visNetworkOutput("networkDorotheplot",height = "800px")),fluidRow(br())
                                                                                        ))
                                                          
                                                       ) # tabpanel network
                                                )# tabsetPanel
                                                
                               )  # conditionalPanel
                               
                               
                  )
                  ,tabPanel(id="xxx",title="MSigDB",
                          
                              tabsetPanel(id="tabboxtabs2",
                                     tabPanel("Table",
                                              fluidRow(
                                                shiny::sidebarPanel(width = 2, fluidRow( align = 'center',actionBttn("collect","Create Table",icon =icon("play"),style = "jelly", color = "warning")),
                                                                    
                                                                    
                                                                    awesomeCheckboxGroup(inputId = "msigselect",label = h4('Select Gene Sets'), status = "warning",
                                                                                         choices = list("H"="H","C1"="C1", "C2_CGP "="CGP", "C2_CP"="CP","C2_CP:BIOCARTA"="BIOCARTA","C2_CP:KEGG"="KEGG",
                                                                                                        "C2_CP:PID"="PID","C2_CP:REACTOME"="REACTOME","C2_CP:WIKIPATHWAYS"="WIKIPATHWAYS","C3_MIR:MIR_Legacy"="MIR_Legacy",
                                                                                                        "C3_MIR:MIRDB"="MIRDB","C3_TFT:GTRD"="GTRD","C3_TFT:TFT_Legacy"="TFT_Legacy",
                                                                                                        "C4_CGN"="CGN","C4_CM"="CM","C5_GO:BP"="BP","C5_GO:CC"="CC","C5_GO:MF"="MF","C5_HPO"="HPO","C6"="C6","C7_IMMUNESIGDB"="IMMUNESIGDB","C7_VAX"="VAX","C8"="C8"),
                                                                                         selected = "H"),
                                                                    actionButton("collec",  tags$b(em("MSigDB Collections")) ,icon = icon("info"), class = "btn-success"),
                                                                    
                                                                    bsModal("modalExample12", "", "collec", size = "large",includeMarkdown("mdfiles/MSigDBCollections.md"))  ),
                                                mainPanel(width = 10,   fluidRow(DT::DTOutput("MSigtable"),br(),downloadButton(outputId = "MSigtabledown", label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;")),
                                                            )),
                                              fluidRow(
                                             column(2,offset=1, fluidRow(h4(" Use dendrogram")),
                                                                fluidRow(   materialSwitch("usedendrog_msig","For samples",value = TRUE, status = "warning")),
                                                                fluidRow(   materialSwitch("usedendrogrow_msig","For genes",value = FALSE, status = "warning"))
                                                                 ),
                                                column(6, dropdown(inputId="dropdownmsig" ,label ="More Parameters for Heatmap", status = "success", width = "800px",circle = FALSE,margin = "20px",
                                                                   icon = icon("folder-open",lib = "glyphicon"),# style = "stretch", status = "default", width = "800px",
                                                                   
                                                                   tooltip = tooltipOptions(title = "Click to select!"), 
                                                                   fluidRow(awesomeRadio("selecteddist2", "Distance method",choices = c("dist (default)"="dist","pearson", "spearman" , "kendall"),inline = TRUE,  status = "warning", checkbox = TRUE)),
                                                                   
                                                                   fluidRow(column(3,sliderInput("col_text_Ang_msig", label = h4("Column text angle"), min = 0, max = 90, value = 45)),
                                                                            column(3,sliderInput("fontsize_rowc_msig", label = h4("Row font size"), min = 1, max = 20, value = 10)),
                                                                            column(3,sliderInput("fontsize_colc_msig", label = h4("Column font size"), min = 1, max = 20, value = 10)) )
                                                )   )),
                                                
                                          
                                              fluidRow(plotlyOutput("msigdbheatmap",height='800px')), fluidRow(column(4,offset=1,verbatimTextOutput("textmisgdb")))
                                              ),
                                             
                              tabPanel("Network",
                                       br(), tags$style(HTML(".tooltip > .tooltip-inner {width: 400px;color: black;background-color: #FCAA1C;}")),  
                                       fluidRow( column(3,  panel(class="myclass3", 
                                                                  fluidRow(align = 'center',sliderInput("fdrcutofmsig"," q-value",min = 0,max = 1,value = 0.05)),
                                                                  fluidRow( align = 'center',sliderInput("termnumber","Term number",min = 2,max = 20,value = 10),bsTooltip("termnumber", "This will be active if the number of nodes are greater than 20 after qvalue selection.",
                                                                                                                                                                           "bottom", options = list(container = "body"))),
                                                                  fluidRow( align = 'center',sliderInput("fontsizenet"," Font size",min = 2,max = 30,value = 12)),
                                                                  fluidRow(align = 'center',selectInput("select_geneset", label = "Select Collections", choices = NULL, multiple = TRUE) ),
                                                                  fluidRow(align = 'center',actionBttn("networkplot", label=tags$div( icon("play"),HTML(" <b>Create Network</b> "),style = "font-size:20px;"),style = "jelly", color = "warning") ),
                                                                 br(),
                                                                  fluidRow( style="font-size: 16px; line-height: 16px;",HTML(paste0( "<p> q-value is selected to draw network for terms having q-value less than selected. <p>",
                                                                    " <p> The width of edges are represented by the number of shared genes between two terms. Thicker edges means more shared genes. <p> ",
                                                                                                 "<p> Node sizes are represented by p-value (bigger nodes means smaller p-values). </p>",
                                                                                                 "<p>If the number of nodes exceeds 20 after applying the q-value cutoff, it limits the visualization to only include the first selected-number of terms. </p>"
                                                                                                 
                                                                            ) ))),
                                                        downloadButton( outputId = "downloadnetwork",label = "Download Netwotk as HTML",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;")
                                                        ),
                                                 column(9,fluidRow(visNetworkOutput("msigdbNetwork", width = "100%", height = "700px")%>% shinycssloaders::withSpinner(type = getOption("spinner.type", 6),
                                                                                                                                                                       color = getOption("spinner.color", default = "#00838f"),
                                                                                                                                                                       size = getOption("spinner.size", default = 1) )
                                                 ) )
                                         )
                                       )
                              
                              )
                  ),
                  tabPanel(title="Compare Clusters",
                           panel(class="myclass1",
                           fluidRow(column(3,sliderInput("font.slider", label = "Font size", min = 2, max = 12, value = 6)),
                                    column(3,sliderInput("logFCthreshold", label = "Absolute logFC threshold", min = 0.0000000001, max = 4, value = 1) ),
                                    column(3, actionBttn("ccaction","Plot Enrichment",icon =icon("bolt"),style = "jelly", color = "warning")))),
                           plotOutput("kegg1",height = "800px")%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))
                           
                  )
      )
  }) 
  
  ####################### PART 1. TABLE ----------------------------------

  # this applies meta analysis for p-values
  result_table <- reactive({
    
 columnames <- input$colname
    
    if (conditions()=="1"){
      result <- data.frame(reac$limma)
      
      matchExpression <- paste(columnames, collapse = "|")
      result2 <- result %>% dplyr::select(matches(matchExpression))
      
      result2 <- data.frame(round(result2,3))
      
      result <- cbind(result[,1],result2)
   
      
      colnames(result)[1] <- "hgnc_symbol"  
      result_table <-  result
      
    }else{ # Meta–analysis will be performed when only more than one data is selected
     
      withProgress(message = 'Performing Meta-analysis of p values...', style = "notification", value = 0.1, {
        Sys.sleep(0.05)
        req(reac$limma)
        # method for Meta–analysis of p–values
        result <- data.frame(reac$limma)
      
        result3 <- result %>% dplyr::select(matches("p.value")) # select columns containing p-values
        logfcres <- result %>% dplyr::select(matches("logFC")) # select columns containing logFC
        empt <- c()
        upc <- c()
        downc <- c()
        metod <- reac$metod 
        print(head(result3))
        for (i in 1:nrow(result3)) {
          k <- as.vector(unlist(unname(result3[i,]))) # make all values in same row as vector
         
          empt <- c(empt,func_meta(k,metod)) # calculate meta-analysis
          l <- as.vector(unlist(unname(logfcres[i,]))) # logFC values of same row
          if(any(k<=0.05)==TRUE){
            sig <- l[which(k<=0.05)]
            up <- length(which(as.vector(sig) > 0)) # number of upregulated genes
            down <- length(which(as.vector(sig) < 0))# number of downregulated genes
          }else {
            up <- 0
            down <- 0}
          upc <- c(upc,up)
          downc <- c(downc,down)
        }
        matchExpression <- paste(columnames, collapse = "|")
        result2 <- result %>% dplyr::select(matches(matchExpression))
        
        result <- cbind(result[,1],round(result2,3),round(empt,3), upc,downc)
        colnames(result)[1] <- "hgnc_symbol"
        colnames(result)[ncol(result)-2] <- "Meta–analysis of p–values"
        colnames(result)[ncol(result)-1] <- c("Sig upregulated")
        colnames(result)[ncol(result)] <- c("Sig downregulated")
        result_table <-  result
        
        incProgress(0.5)
        shiny::setProgress(1)
      })
    }
    result_table
  })
  ### KEGG patway genes. or user upload genes
  d <- reactive({
    selectedpat <- input$in4
    print(is.null(input$in4))
    
    dt <- data.frame(result_table())
    if(input$in4!="Display all genes"){
      mylist_sub <- keggpat[grep(selectedpat, names(keggpat))]
      airSE <- idMap(mylist_sub, org = "hsa", from = "ENTREZID", to = "SYMBOL")
      patgenes <- airSE[[1]]
      dt2 <- dt[dt[,1]%in%patgenes,   ]
    }
    else{dt2=dt}

    
    if(input$usegenes==TRUE){
    
      if(!is.null(input$genefile)){
      
        tryCatch(
          {
            #patgenes2 <- read.csv(input$genefile$datapath,header = input$header)
            
            patgenes2 <- data.frame(read_delim(input$genefile$datapath,
                       delim = input$sep, escape_double = FALSE, col_names = input$header, 
                       trim_ws = TRUE))
            
         print(head(patgenes2))
            dt2 <- dt2[dt2[,1]%in%(patgenes2[,1]),   ]
            print(dt2)
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          }
        )
      }else {return(dt2)}
      
      
      
    }else {return(dt2)}
    
    
    
  })
  # column display
  result_tbl <- reactive({
    columnames <- input$colname 
    
    if (conditions()!="1"){
     
      result <- result_table()
      
      if (input$in4=="Display all genes" ){
        if(input$usegenes==FALSE){
          result2 <- result
        }
        else {
          result2 <- d()
        }
        
      }
      else{
        result2 <- d()
      }
    }
    
    else {    # Meta–analysis will be performed when only more than one data is selected
      #result <- data.frame(reac$limma)
      result <- data.frame( result_table())
      matchExpression <- paste(columnames, collapse = "|")
      result2 <- result %>% dplyr::select(matches(matchExpression))
      result2 <- data.frame(round(result2,3))
      result <- cbind(result[,1],result2)
      colnames(result)[1] <- "hgnc_symbol"
      
      if (input$in4=="Display all genes" ){
        if(input$usegenes==FALSE){
          result2 <- result
        }
        else {
          result2 <- d()
        }
        
      }
      else{
        result2 <- d()
      }
    }
    
  })

  output$tablefirst <-DT::renderDT({
   
    result2 <- data.frame(result_tbl())
   
   
    DT::datatable(result2, escape=F,plugins = 'accent-neutralise',
                  rownames=T,filter = 'top',caption = titles(),class = "compact hover row-border",
                  extensions = c('Scroller','Select', 'Buttons'),
                  callback=JS('$("button.buttons-copy").css("background","#cedbd3");
                    $("button.buttons-excel").css("background","#cedbd3");
                     $("button.selectAll").css("background","#cedbd3");
                      $("button.deselectAll").css("background","#cedbd3");
                    return table;'),
                  options = list(initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#cedbd3', 'color': '164036'});","}"),
                                 select = list(style = "multi", items = "row"), autoWidth = TRUE,
                                 columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                 language = list(info = 'Showing _START_ to _END_ of _TOTAL_ variables'),
                                 deferRender = TRUE,

                                 scrollX = TRUE,scrollY = 500,scroller = TRUE,dom = "Blfrtip",
                                 buttons = list('copy', 'excel',list(extend='selectAll',className='selectAll',
                                                                            text="Select All Rows",
                                                                            action=DT::JS("function(e, table, node, config) {",
                                                                                          "  table.rows({ search: 'applied'}).deselect();",
                                                                                          "  table.rows({ search: 'applied'}).select();",
                                                                                          "}")
                                 ), list(extend='selectNone',className='deselectAll',text="Deselect All",
                                         action=DT::JS("function(e, table, node, config) {",
                                                       "  table.rows({ search: 'applied'}).select();",
                                                       "  table.rows({ search: 'applied'}).deselect();",
                                                       "}")))),selection="none"
    )%>%
      formatStyle(columns =colnames(result2)[which(grepl("logFC",colnames(result2)))], target = "cell", backgroundColor = "#cedbd3")


    
  },server=FALSE) 
  
  
  output$d1_4 <- shiny::downloadHandler(
    filename = function(){"Result of selected data.csv"}, 
    content = function(fname){
      write.csv(format(reac$limma,decimal.mark=",",big.mark="."), fname)
    }
    
  )
  
  reactivedata <- reactiveValues(subsetdt=NULL)
  
  data.subset=shiny::eventReactive(input$tablefirst_rows_selected,{
    
    selected <- input$tablefirst_rows_selected
    
    
    if (input$in4=="Display all genes" ){
      if(input$usegenes==FALSE){
        result <- data.frame( result_table())
      }
      else {
        result <- d()
      }
      
    }
    else{
      result <- d()
    }
    if(length(result)){
      newresult=result[selected , ]
    }
    reactivedata$subsetdt=newresult
    newresult
    
  })  # reactive subset
  output$subsettable <- DT::renderDT({
    print( reactivedata$subsetdt)
    newresult <- data.subset()
    matchExpression <- paste(input$colname, collapse = "|")
    newresult2 <- newresult %>% dplyr::select(matches(matchExpression))
    newresult <- cbind(newresult[,1],newresult2)
    colnames(newresult)[1] <- "hgnc_symbol"
    DT::datatable(
      newresult,
      escape = FALSE,
      callback=JS('$("button.buttons-copy").css("background","#cedbd3"); 
                    $("button.buttons-excel").css("background","#cedbd3"); 
                    return table;'),
      extensions = c('Select', 'Buttons', 'ColReorder'),filter = "top",
      class = "compact hover row-border",
      options = list(
        columnDefs = list(list(className = 'dt-center', targets = "_all")),
        colReorder = TRUE,
        lengthMenu = c(10, 30, 50,100,200,500),
        pageLength = 10,
        scrollX = TRUE,
        dom = "Blfrtip",
        buttons = list('copy',  'excel' ), initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color':'#cedbd3', 'color': '164036'});",
          "}")
      )
    ) %>% 
      formatStyle(columns =colnames(newresult)[which(grepl("logFC",colnames(newresult)))], target = "cell", backgroundColor = "#cedbd3") 
    
    
  })
  output$subset_d1_4 <- downloadHandler(
    filename = function(){"Result of selected genes.csv"}, 
    content = function(fname){
      write.csv(format(data.subset(),decimal.mark=",",big.mark="."), fname)
    }
  )
  
  ####################### PART 2. Scatter plot ----------------------------------
  x_axis <- reactive({input$select_xaxis})
  y_axis <- reactive({input$select_yaxis})
  updatednames <- reactive({
    dt <- data.frame(result_table())
    
    colnames(dt) <- str_replace_all(colnames(dt), "[..]", " ")
    
    updatednames <- c(colnames(dt)[grep("logFC",colnames(dt))],
                      colnames(dt)[grep("p value",colnames(dt) )])
  })
  scatter_plot <- shiny::eventReactive(input$scatter,{
    # if user wants to use all genes, it will use result_table(), but if user wants to use subsetted genes, it uses data.subset()
    if(input$scatter_material==FALSE){
      dt <- data.frame(result_table())
    }else {  dt <- data.subset() }
  
    # column names needs to be adjusted so that selected x-y axis names can be searched in them.
    # example: column name at first:"Drsp Cont_our_MCF7..p.value"
    colnames(dt) <- str_replace_all(colnames(dt), "[..]", " ")# remove any .. 
    
    # do same for selected x and y axis names, remove space
    xx <- str_replace_all(x_axis(), "[..]", " ")
    yy <- str_replace_all(y_axis(), "[..]", " ")
    
    # search selected axis names in columns of data
    x=data.frame(dt)[,grep(xx,colnames(dt))]
    
    y=data.frame(dt)[,grep(yy,colnames(dt))]
    sze=as.numeric(input$dotsize)
    # plot interactive plot 
    scatter_plot <-  plot_ly(data = data.frame(cbind(x,y)), x = ~x, y = ~y,text = ~dt[,1] , textposition = 'middle right',
                             marker = list(size = sze,color = 'rgba(255, 182, 193, .9)',
                                           line = list(color = 'rgba(152, 0, 0, .8)',width = 2)))%>%
      layout(plot_bgcolor = "#EEF2F7", xaxis = list(title = x_axis()), 
             yaxis = list(title = y_axis()))
    
  })
  output$scatterplot <- renderPlotly({
    shiny::validate(shiny::need(input$scatter>0, "Click Draw to see the plot."))
    scatter_plot()
  })
  
  ####################### PART 3. HEATMAP  ----------------------------------
  distfun_row = function(x) stats::dist(x, method = input$distFun_row)
  distfun_col =  function(x) stats::dist(x, method = input$distFun_col)
  
  hclustfun_row = function(x) stats::hclust(x, method = input$hclustFun_row)
  hclustfun_col = function(x) stats::hclust(x, method = input$hclustFun_col)
  
  output$heatmapplot <- plotly::renderPlotly({
    shiny::validate(
      need(!is.null(reactivedata$subsetdt), "Please select genes to create heatmap!"))
    NULL
  })
  ##### heatmap
  shiny::observeEvent(data.subset(),{
    req(!is.null(reactivedata$subsetdt))
    interactiveHeatmap<- shiny::reactive({
      shiny::validate(need(base::nrow(data.subset())<500,"If the number of selected genes is more than 500, the heatmap works very slowly. Please select less number of genes (<500)"),
                      need(length(grep( "logFC" , colnames( data.subset() ) ))>1,"Heatmap cannot be drawn when there is only one comparison selected."    )
                     )
   
      
      df <- data.subset()
      onlyFC <- df[ , grepl( "logFC" , colnames( df ) ) ]
      rownames(onlyFC) <- df[,1]
      experiment_names=colnames(onlyFC)
      colnames(onlyFC)=paste(rep("Data ",length(colnames(onlyFC))),1:length(colnames(onlyFC))) # make col names like data1 data2...
      tit=paste(paste(colnames(onlyFC),experiment_names,sep =": " ),collapse = "\n")
      
      p <- heatmaply:: heatmaply(onlyFC, 
                                 distfun_row =  distfun_row,
                                 hclustfun_row = hclustfun_row,
                                 distfun_col = distfun_col,
                                 hclustfun_col = hclustfun_col,key.title=tit,label_names = c("Gene", "Data", "LogFC"),
                                 scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                                   low = "blue", 
                                   high = "red", 
                                   midpoint = 0
                                 ),
                                 k_row = input$r)
      p
    })
    
    
    output$heatmapplot <- plotly::renderPlotly({
     
      interactiveHeatmap()
      
    })
    
    output$downloadData <- shiny::downloadHandler(
      filename = function() {
        paste0("heatmaply-", strftime(Sys.time(),'%Y%m%d_%H%M%S'), ".html")
      },
      content = function(file) {
        h <- interactiveHeatmap()
        h$width='100%'
        h$height='800px'
        libdir <- paste0(tools::file_path_sans_ext(basename(file)),"_files")
        htmltools::save_html(htmltools::browsable(htmltools::tagList( h)),file=file,libdir = libdir)
        rmarkdown::pandoc_self_contained_html(file, file)
        unlink(libdir, recursive = TRUE)
      }
    )
  })
  ####################### PART 4. GOST plot------------
  gostplotre <- eventReactive(input$gostbut,{
    withProgress(message = 'Generating Gost plot...', style = "notification", value = 0.1, {
      Sys.sleep(0.3)
      isolate({
        tryCatch({
          newresult <- isolate(data.subset())
          print(dim(newresult))
          rownames(newresult) <- newresult[,1]
          
          merged <- merge(newresult,gene.list, by="row.names")
          print(dim(merged))
          incProgress(0.5)
          p = gost(merged$ensembleID, organism = "hsapiens")
          
        },
        error = function(e) {
          # showNotification(id="errorNotify1", "Error is detected, rerun the analysis", type = "error", duration = NULL)
          return(NULL)
        }
        
        )
      })
      return(p)
    })
  })
  output$gost1 <- renderPlotly({
    shiny::validate(shiny::need(!is.null(data.subset()),"Please select genes first!"),
                    shiny::need(input$gostbut>0,"Click 'Run'!")
                    )
    p <- gostplotre()
    if(is.null(p$result$term_name))
    {    shiny::validate("Please select more/different genes!")
      return(NULL)
    }
    else {
      gostplot(p, interactive = T)
    }
    
    
  })
  
  levelselect <- reactiveValues()
  observeEvent(input$got,
               {
                 levelselect$lev <- input$selectlevel1
                 levelselect$ont <- input$selectont1
                 levelselect$data <- data.subset()
               }
  )
  
  ####################### PART 5. GO   ----------
  ### data update for the selection of LogFC---
  output$LogFC1UI <- renderUI({
    req(input$analyze)
    dt <- data.frame(result_table())
    sub <- dt[grep("logFC",colnames(dt))]
    tagList(awesomeRadio("LogFC1", label = h4("Use logFC from data"),choices =  c(1:ncol(sub)),selected = 1, inline = TRUE, status = "warning" ,checkbox = TRUE))
  })
  
  output$goLogFC1UI <- renderUI({
    dt <- data.frame(result_table())
    sub <- dt[grep("logFC",colnames(dt))]
    tagList(awesomeRadio("selectlogFCGO", label = "Use logFC from data",choices =  c(1:ncol(sub)),inline = TRUE, status = "warning", checkbox = TRUE))
  })
  # ---- GO
  goenrich <- eventReactive(input$got,{
    req(input$tabsetss == "Functional Profile")
    ggo <- NULL
    
    withProgress(message = 'Finding GO terms...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      if(is.null(input$tablefirst_rows_selected)){
        showNotification("Please select genes from Analysis tab!")
      }else{
        print("dsfdfdf")
        newresult <- data.subset()
        rownames(newresult) <- newresult[,1]
        merged <- merge(newresult,gene.list, by="row.names")
        cha <- unlist(as.character(merged$entrezgene_id))
        ggo <- clusterProfiler::enrichGO(gene= as.vector(cha),
                                         OrgDb         = "org.Hs.eg.db",
                                         ont           = "ALL",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.05,
                                         minGSSize =levelselect$lev ,
                                         readable      = TRUE)
        
      }
      incProgress(0.5)
      shiny::setProgress(1)
    })
    ggo
  })
  output$GOtable1 <- DT::renderDataTable({
  req(input$got)
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
  
    shiny::validate(shiny::need(input$got>0,"Click 'Create GO table'!"),
                    shiny::need(try(is.null(goenrich()) == FALSE), "No gene can be mapped. Please select different or more genes."))
    ggo <- goenrich()
    print("GO table")
    ggo <- data.frame(ggo)
    ggo[,6] <-data.frame(round(ggo[,6],3))
    ggo[,7] <-data.frame(round(ggo[,7],3))
    ggo[,8] <- data.frame(round(as.numeric(ggo[,8]),5))
    print(input$tabsetss)
    #ggo <- ggo[ggo[,3]!=0,]
    DT::datatable(
      as.data.frame(ggo),
      escape = FALSE,rownames = F,
      extensions = c('Select', 'Buttons', 'ColReorder'),filter = "top",
      callback=JS('$("button.buttons-copy").css("background","#cedbd3");
                    $("button.buttons-excel").css("background","#cedbd3");
                    $("button.buttons-csv").css("background","#cedbd3");
                     $("button.buttons-print").css("background","#cedbd3");
                    return table;'),
      options = list(
        colReorder = TRUE,
        lengthMenu = c(10, 30, 50,100),
        pageLength = 10,
        scrollX = TRUE,
        #autoWidth = TRUE,
        dom = 'Blfrtip',
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color':  '#cedbd3', 'color': '164036'});",
          "}")
      ))
  })
  # reactiveValues
  num <- reactiveValues()
  observeEvent(input$dotplotbttn,
               {
              
                 num$n <- as.numeric(input$slid)
                 num$l <- input$layot
                 
               })
  
  
  
  
  observe ({
    if (req(input$tabsetss) == "Functional Profile"){
      
      observeEvent(input$dotplotbttn,{
    
        
        output$emapplot <- renderPlot({
          shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
          num <-num$n
          edo <- pairwise_termsim(goenrich())
          enrichplot::emapplot(edo, cex_category=1.5, showCategory=num)+ theme(legend.text=element_text(size=16, face = "bold"),legend.title=element_text(size=16, face = "bold"))  
        })
        output$dotplt <- renderPlot({
          shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
          num <-num$n
          ego <- goenrich()
          enrichplot::dotplot(ego, showCategory=num)#+ theme(legend.text=element_text(size=16, face = "bold"),legend.title=element_text(size=16, face = "bold"))  
        })
        output$cnplt <- renderPlot({
          shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
          k <- num$l
          num <-num$n
          n <- data.subset()
          sub <- n[grep("logFC",colnames(n))]
          s <- sub[,as.numeric(input$selectlogFCGO)]
          names(s) <- n[,1]
          ego <- goenrich()
          enrichplot::cnetplot(ego,foldChange=s,showCategory=num,layout = k) + 
            scale_color_gradient2(name='LogFC', low='darkgreen', high='red',mid = "yellow")+ 
            theme(legend.text=element_text(size=16, face = "bold"),legend.title=element_text(size=16, face = "bold"))  
          
        })
        
        
      })
    }
  })
  ####################### PART 6. Disgenet   ---------- 
  
 
  disgenetnetworks <- eventReactive(input$actionfordisgenet,{
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"),
                     shiny::need(input$actionfordisgenet>0,"Click 'Draw Network'!")
                    )
    withProgress(message = 'Performing enrichment analysis based on the DisGeNET ...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      isolate({
        tryCatch({
          newresult <- data.subset()
          rownames(newresult) <- newresult[,1] # rownames are now hgnc symbols
          merged <- merge(newresult,gene.list, by="row.names")
          
          ent <- as.character(merged$entrezgene_id)
          res <- DOSE::enrichDGN(ent,minGSSize = input$minGSSize1, qvalueCutoff = input$qvaluecutof,pvalueCutoff = input$pvaluecutof,readable = TRUE)
          
          
          if(is.null(res@result[["ID"]]))
          {
            showNotification(id="warnNotify", "No gene can be mapped ...", type = "warning", duration = NULL)
            showNotification(id="warnNotify2", "Change the parameters inside dropdown button and try again.", type = "warning", duration = NULL)
            return(NULL)
          }
          
        },
        error = function(e) {
          
          #showNotification(id="errorNotify1", "Error is detected, rerun the analysis", type = "error", duration = NULL)
          
          return(NULL)
        }
        
        )
      })
      return(res)
    })
  })
  
  output$network1 <- renderPlot({
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"),
                     shiny::need(input$actionfordisgenet>0,"Click 'Draw Network'!"),
                     need(!is.null(disgenetnetworks()),"No gene can be mapped. Please select more/different genes!"))
    
    res <- disgenetnetworks()
    
    newresult <- data.subset()
    rownames(newresult) <- newresult[,1] # rownames are now hgnc symbols
    merged <- merge(newresult,gene.list, by="row.names")
    x <- as.numeric(input$LogFC1)
    merged2=data.frame(merged)[,grep("logFC",colnames(merged))]
 
    if (conditions()=="1"){
      geneList=merged2
    }
    else{
      geneList <- unlist(merged2[,x])
    }

    ent <- as.character(merged$entrezgene_id)
    names(geneList) <- ent
    enrichplot::cnetplot(res, categorySize="pvalue", foldChange=geneList,showCategory=input$showCategory,
                         circular=as.logical(input$circular), colorEdge=T,
                          #color.params = list(edge = TRUE),
                       node_label="all")
    
    
    
    
    
  })
  output$heatplot1 <- renderPlot({
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!") ,  
                     shiny::need(input$actionfordisgenet>0,"Click 'Draw Network'!"),
                     need(!is.null(disgenetnetworks()),"No gene can be mapped. Please select more/different genes!"))
    res <- disgenetnetworks()
    
    newresult <- data.subset()
    rownames(newresult) <- newresult[,1] # rownames are now hgnc symbols
    merged <- merge(newresult,gene.list, by="row.names")
    ent <- as.character(merged$entrezgene_id)
    x <- as.numeric(input$LogFC1)
    merged2=data.frame(merged)[,grep("logFC",colnames(merged))]

    if (conditions()=="1"){
      geneList=merged2
    }
    else{
      geneList <- unlist(merged2[,x])
    }
    
    names(geneList) <- ent
    # p <-  enrichplot::cnetplot(res, categorySize=input$categorySize, foldChange=geneList,showCategory=input$showCategory,circular=as.logical(input$circular),  colorEdge = TRUE,node_label="all")
    enrichplot::heatplot(res,showCategory=input$showCategory, foldChange=geneList)
    
    
    
    
  })
 
  ####################### PART 7. DoRothea  ----------------------------------
  netwrk <- reactiveValues()
  observeEvent(input$networkDorothe,
               {
                 netwrk$fdrcutof <- input$fdrcutof
                 netwrk$tfnumber=input$tfnumber
                 netwrk$nodecolor=input$nodecolor
                 
               }
  )
  

output$nodecolorui <- renderUI({
  req(input$analyze)
  dt <- data.frame(result_table())
  sub <- dt[grep("logFC",colnames(dt))]
  dats=paste(rep("Data ",ncol(sub)),1:ncol(sub))
  tagList(awesomeRadio("nodecolor", label = "Node color",choices =  c(dats),selected = dats[1], inline = TRUE, status = "warning" ,checkbox = TRUE))
})


output$differencecolor <- renderUI({
  req(input$analyze)
  if (conditions()=="1"){
    NULL
  }else{
    
 
  dt <- data.frame(result_table())
  sub <- dt[grep("logFC",colnames(dt))]
  dats=paste(rep("Data ",ncol(sub)),1:ncol(sub))
  tagList(column(6, selectInput(inputId = "firstdataset",label = "First Data",choices = dats,selected = dats[1])),
          column(6, selectInput(inputId = "seconddataset",label = "Second Data", choices = dats,selected = dats[2])))
  }
})

  listfinal <- readRDS("genesets.RDS")
  listfinalALL <- readRDS("genesetsALL.RDS")
  net <- readRDS("net.RDS")
  alldataTF <- readRDS("alldataTF.RDS")
 
  
  observeEvent(input$DoRothEAbutton,{
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!") )
    withProgress(message = 'Calculating overlapped genes...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      
      analysistable <- result_table()
      #  net <- dorothea::dorothea_hs
      analysistable <- analysistable[analysistable[,1]%in% unique(net$target),]
      bggens <- dim(analysistable)[1]
      ####-----------
      data.subset <- data.subset()
      base::rownames(data.subset) <- data.subset[,1]
      if (!is.null(input$tablefirst_rows_selected)){
        if(input$selectlist=="all"){ # if all genes not only NF targets are used
          hyp_obj <- hypeR( base::rownames(data.subset), listfinalALL,background=20295)
        }
        else {
          if(input$selectbackground=="allgenes"){
            hyp_obj <- hypeR( base::rownames(data.subset), listfinal,background=20295)
          }else { hyp_obj <- hypeR( base::rownames(data.subset), listfinal,background=8338)}
          
        }
        
        dt <- data.frame(hyp_obj[["data"]])
        dt <- dt[dt$overlap>0,]
        
        colnames(dt) [1]<-"TF"
    
        output$DoRothEAtable <- DT::renderDataTable({
          shiny::validate(shiny::need(!is.null(dt),"No overlaps found for selected genes"),
                          shiny::need(!is.null(data.subset()),"Please select genes first!"))
      
          DT::datatable( dt,rownames= FALSE,extensions = c( 'Buttons', 'ColReorder'),filter = "top",
                         callback=JS('$("button.buttons-copy").css("background","#cedbd3"); 
                    $("button.buttons-excel").css("background","#cedbd3"); 
                     $("button.buttons-pdf").css("background","#cedbd3"); 
                      $("button.buttons-colvis").css("background","#cedbd3"); 
                    return table;'),
                         options = list(initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#cedbd3', 'color': '164036'});","}"),
                           colReorder = TRUE,
                           lengthMenu = c(10, 30, 50,100,200,500),
                           pageLength = 10,
                           scrollX = TRUE,
                           dom = "Blfrtip",
                           buttons = list('colvis','copy', 'excel', 'pdf')
                         )
          )
        })
        
      }
      else{
        showNotification("Please select genes from Analysis tab!")
      }
      incProgress(0.5)
      shiny:: setProgress(1)
    })
    ##### msigdb 
    observeEvent(input$DoRothEAtable_rows_selected,{
    selected <-input$DoRothEAtable_rows_selected # selected rows from table
    hit.genes <- as.vector(dt[selected,8]) # select the genes of rows clicked
    newVec <- unlist(sapply(hit.genes, strsplit, "\\s+", USE.NAMES = FALSE)) 
    xx <- gsub(",", " ", newVec) ; xxx <-unique( unlist(strsplit(xx," "))) # separate the genes as char vector
    
    msigdt <- msigdbr(species = "Homo sapiens", category = "H" )
    nn = msigdt  %>% dplyr::distinct(gs_name, gene_symbol,gs_description) %>% as.data.frame()
    
    hh <- enricher(gene = unlist(xxx ), TERM2GENE = nn)
    shiny::validate(shiny::need(!is.null(hh),"No overlaps found for selected genes! Please select more genes or different collections.")   )
    h <- data.frame(hh@result)[,-c(2,5,6,7)]
    h[,5] <- data.frame(round(h[,5],3))
    h[,4] <- data.frame(round(h[,4],3))
    h[,6] <- data.frame(round(h[,6],3))
    nn2<-unique(nn[,c(1,3)]) # unique description of gene sets
    
    h=h %>% dplyr::left_join(nn2, by = c("ID" = "gs_name")) 
  
    output$MSigtableTF <- DT::renderDT({
      shiny::validate(need(length(input$DoRothEAtable_rows_selected)>0,  "Select a row first!"),
                      
                      shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!")
       )
     
      
      tf.genes <- as.vector(dt[selected,1])  # TFs
      tf.genes <- unlist(sapply(tf.genes, strsplit, "\\s+", USE.NAMES = FALSE)) 
      print(tf.genes)
      tf.genes= paste(tf.genes,collapse=', ')
      tit=paste("Selected TF: ",tf.genes )
      DT::datatable( h,rownames= FALSE,extensions = c( 'Buttons', 'ColReorder'),filter = list(position = 'top'),class='stripe cell-border hover', caption = tit,
                     callback=JS('$("button.buttons-copy").css("background","#cedbd3"); 
                    $("button.buttons-excel").css("background","#cedbd3"); 
                     $("button.buttons-pdf").css("background","#cedbd3"); 
                      $("button.buttons-colvis").css("background","#cedbd3"); 
                    return table;'),
                     options = list(initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#cedbd3', 'color': '164036'});","}"),
                                    colReorder = TRUE,
                                    lengthMenu = c(5,10, 30, 50,100,200,500),
                                    pageLength = 10,
                                    scrollX = TRUE,
                                    dom = "Blfrtip",
                                    buttons = list('colvis','copy', 'excel', 'pdf')
                                    
                                    
                     )
      )
    })
    output$msigdbTFheatmap <- plotly::renderPlotly({
      shiny::validate(need(length(input$DoRothEAtable_rows_selected)>0,  ""),
        shiny::need(!is.null(input$MSigtableTF_rows_selected),"Please select terms first!"),
        shiny::need(!is.null(hh),"No overlaps found for selected genes! Please select more genes or different collections."),
        need(conditions()!= "1", "To create heatmap, please select more data comparison."))
      req(input$MSigtableTF_rows_selected )
   
      selectedrow<-input$MSigtableTF_rows_selected # selected rows from table
      hit.genes <- as.vector(h[selectedrow,7]) # select the genes of rows clicked
      newVec <- unlist(sapply(hit.genes, strsplit, "/", USE.NAMES = FALSE)) 
  
      
      logFCfor.selectedgenes <- data.subset[data.subset[,1] %in% newVec, ] # find logFC values of selected genes
   
      
      data.heatmap <- logFCfor.selectedgenes[ , grepl( "logFC" , colnames( logFCfor.selectedgenes ) ) ]
      base::rownames(data.heatmap) <- logFCfor.selectedgenes[,1]
      # for annotation part in heatmap
      selected.geneset <- h[selectedrow,c(1,7)]
      
      for (i in 1:nrow(data.heatmap)){
        gene <- base::rownames(data.heatmap)[i] 
        for (j in 1:nrow(selected.geneset)){
          
          resultofcheck=match(1, str_detect( unlist(sapply(as.vector(selected.geneset[j,2]), strsplit, "/", USE.NAMES = FALSE)), gene))
          if (!is.na(resultofcheck)){
            data.heatmap[i,selected.geneset[j,1]]=selected.geneset[j,1]
          }
        }
      }
      
      data.heatmap <- data.frame(data.heatmap)
      
      col.logfc <- grep("logFC", colnames(data.heatmap)) # select the columns containing logFC
      
      selecteddist2 <- input$selecteddist2
      dataX= data.heatmap[,col.logfc]
      colnames(dataX)=paste(rep("Data ",length(colnames(dataX))),1:length(colnames(dataX))) # make col names like data1 data2...
      tit=paste(paste(colnames(dataX),colnames(data.heatmap[,col.logfc]),sep =": " ),collapse = "\n")
      
  
      
      if(base::nrow(data.heatmap)>1 && base::nrow(data.heatmap)<500){
        
        if ((ncol(data.heatmap)-length(col.logfc)) >1){
          annotpart <-  data.heatmap[,-col.logfc]
          
          abbreviatednames= abbreviate(gsub("_"," ",colnames(annotpart)), minlength = 1, use.classes = F,
                                       dot = F, strict = FALSE,
                                       method =  "left.kept", named = TRUE)
          colnames(annotpart)=abbreviatednames # make abreviation for annotation part in colnames
          # annotlabel=paste(paste(abbreviatednames,colnames(data.heatmap[,-col.logfc]),sep =": " ),collapse = "\n")
          heatmaply:: heatmaply(as.matrix(dataX[,col.logfc]), 
                                Colv=TRUE,Rowv=TRUE, 
                                row_side_colors =  annotpart,
                                distfun= stats::dist,
                                
                                column_text_angle = 45,
                                fontsize_row = 10,
                                fontsize_col =10,
                                
                                key.title=tit,
                                label_names = c("Gene", "Data", "LogFC"),
                                #  key=TRUE,
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
                                # colors = colr #bluered(length(as.matrix(dataX))), plot_method="plotly",
          )
          #%>%layout(showlegend = TRUE)
        }
        else {  
          heatmaply:: heatmaply(as.matrix(dataX[,col.logfc]), Colv=TRUE,Rowv=TRUE, distfun= stats::dist,key.title=tit,
                                column_text_angle = 45,
                                fontsize_row = 10,
                                fontsize_col = 10,
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
          )
          
        }
        
        
        
      }else{
        showNotification(paste ("Select different terms from the above table which contains different genes!", "Heatmap cannot be drawn if there are more than 500 genes!"))
        shiny::validate("Heatmap cannot be drawn if there are more than 500 genes!")
        return(NULL)}
      
    })
    output$textmisgdbTF <- renderText({
      req(input$MSigtableTF_rows_selected )
 
      df <- data.subset()
      base::rownames(df) <- df[,1]
      selected <-input$MSigtableTF_rows_selected # selected rows from table
      hit.genes <- as.vector(h[selected,7]) # select the genes of rows clicked
      newVec <- unlist(sapply(hit.genes, strsplit, "/", USE.NAMES = FALSE)) 
      logFCfor.selectedgenes <- df[df[,1] %in% newVec, ] # find logFC values of selected genes
      
      
      data.heatmap <- logFCfor.selectedgenes[ , grepl( "logFC" , colnames( logFCfor.selectedgenes ) ) ]
      base::rownames(data.heatmap) <- logFCfor.selectedgenes[,1]
      # for annotation part in heatmap
      selected.geneset <- h[selected,c(1,7)]
      
      for (i in 1:nrow(data.heatmap)){
        gene <- base::rownames(data.heatmap)[i] 
        for (j in 1:nrow(selected.geneset)){
          
          resultofcheck=match(1, str_detect( unlist(sapply(as.vector(selected.geneset[j,2]), strsplit, "/", USE.NAMES = FALSE)), gene))
          if (!is.na(resultofcheck)){
            data.heatmap[i,selected.geneset[j,1]]=selected.geneset[j,1]
          }
        }
      }
      
      data.heatmap <- data.frame(data.heatmap)
      
      col.logfc <- grep("logFC", colnames(data.heatmap)) 
      if(base::nrow(data.heatmap)>1 && base::nrow(data.heatmap)<500){
        
        if ((ncol(data.heatmap)-length(col.logfc)) >1){
          annotpart <-  data.heatmap[,-col.logfc]
          abbreviatednames= abbreviate(gsub("_"," ",colnames(annotpart)), minlength = 1, use.classes = F,
                                       dot = F, strict = FALSE,
                                       method =  "left.kept", named = TRUE)
          HTML(paste(paste(abbreviatednames,colnames(annotpart) ,sep =": " ),collapse  = "\n"))
        }}
      
      
    })
    })
    ########### Heatmap dorothea
    
    output$heatmapplotdorothea <- plotly::renderPlotly({
      shiny::validate(need(length(input$DoRothEAtable_rows_selected)>0,  "Select a row to plot heatmap!"),
                      need(conditions()!= "1", "To create heatmap, please select more data comparison."),
                      shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
      
      selected <-input$DoRothEAtable_rows_selected # selected rows from table
      hit.genes <- as.vector(dt[selected,8]) # select the genes of rows clicked
      newVec <- unlist(sapply(hit.genes, strsplit, "\\s+", USE.NAMES = FALSE)) 
      xx <- gsub(",", " ", newVec) ; xxx <- unlist(strsplit(xx," ")) # separate the genes as char vector
      
      NRtargets <- data.subset[data.subset[,1] %in% xxx, ] # find logFC values of selected genes.  filters the rows of a data frame data.subset where the value in the first column (i.e., data.subset[,1]) is present in a vector xxx. 
      
      dorothehatmap <- NRtargets[ , grepl( "logFC" , colnames( NRtargets ) ) ]
      base::rownames(dorothehatmap) <- NRtargets[,1]
      # for annotation part in heatmap
      net <- data.frame(net)
      net2 <- net[net[,3] %in% xxx,]
      selected.TFs <- dt[selected,1]
      net2 <- net2[net2[,1] %in% selected.TFs, ]
      
      for(j in selected){
        for (i in 1: nrow(dorothehatmap))  {
          net3 <- net2[net2$tf==dt[j,1],]
          target.gene <- rownames(dorothehatmap)[i]
          target.gene.tf <- net3[net3$target==target.gene,]
          
          if(!is.na(target.gene.tf[1,1])){
            
            dorothehatmap[i,dt[j,1]]=(target.gene.tf$confidence)
          }
        }
      }
      
      dorothehatmap <- data.frame(dorothehatmap)
     
      col.logfc <- grep("logFC", colnames(dorothehatmap)) # select the columns containing logFC
      
      
      dataX= dorothehatmap[,col.logfc] # only logFC columns
      colnames(dataX)=paste(rep("Data ",length(colnames(dataX))),1:length(colnames(dataX)))
      
      selecteddist <- input$selecteddist
      if (input$usedendrogrow==TRUE){ dendo=TRUE } else {dendo=FALSE}
      if (input$usedendrog==TRUE){ dendoColv=TRUE } else {dendoColv=FALSE}
      tit=paste(paste(colnames(dataX),colnames(dorothehatmap[,col.logfc]),sep =": " ),collapse = "\n")
  
      #colorRampPalette(c("blue", "white","red"),space = "Lab",bias=bia)(length(a)-1)
      #colr=c(colorRampPalette(c("blue","white"),space = "Lab")(length(a1)), colorRampPalette(c("white","red"),space = "Lab")(length(a2)))
      if(base::nrow(dataX)>1 && base::nrow(dataX)<500){
      if ((ncol(dorothehatmap)-length(col.logfc)) >1){
          annotpart <-  dorothehatmap[,-col.logfc]
          
          heatmaply:: heatmaply(as.matrix(dataX), Colv=dendoColv,Rowv=dendo,
                              
                                column_text_angle = input$col_text_Ang,
                                fontsize_row = as.numeric(input$fontsize_rowc),
                                fontsize_col = as.numeric(input$fontsize_colc),
                                side_color_colorbar_len = 0.4,
                                colorbar_xpos = 1.0,
                                colorbar_ypos = 0.65,
                                row_side_colors =  annotpart,row_side_palette = c(" " = "#E0E0E0", "A" = "#4BDE4B",  "B" = "#FF8000", "C" = "#33DFC3", "D" = "#B74BDE","E" = "#4B5ADE"),
                                distfun=selecteddist,  key.title="logFC",label_names = c("Gene", "Data", "LogFC"),
                                key=TRUE,
                                #colors = colr, plot_method="plotly",
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
                                )#%>%plotly::layout(showlegend = TRUE)
        }
        else {  
          heatmaply:: heatmaply(as.matrix(dataX), Colv=dendoColv,Rowv=dendo, 
                                column_text_angle = input$col_text_Ang,
                                fontsize_row = as.numeric(input$fontsize_rowc),
                                fontsize_col = as.numeric(input$fontsize_colc),
                                distfun=selecteddist,  key.title="logFC",label_names = c("Gene", "Data", "LogFC"),
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
                                
                                )}
      }else{
        showNotification(paste ("Select different terms from the above table which contains different genes!", "Heatmap cannot be drawn if there are more than 500 genes!"))
        shiny::validate("Heatmap cannot be drawn if there are only one gene or more than 500 genes!")
        return(NULL)}
      
    })
    output$textTFsheatmap <- renderText({
      req(input$DoRothEAtable_rows_selected )
      selected <-input$DoRothEAtable_rows_selected # selected rows from table
      hit.genes <- as.vector(dt[selected,8]) # select the genes of rows clicked
      newVec <- unlist(sapply(hit.genes, strsplit, "\\s+", USE.NAMES = FALSE)) 
      xx <- gsub(",", " ", newVec) ; xxx <- unlist(strsplit(xx," ")) # separate the genes as char vector
      
      NRtargets <- data.subset[data.subset[,1] %in% xxx, ] # find logFC values of selected genes.  filters the rows of a data frame data.subset where the value in the first column (i.e., data.subset[,1]) is present in a vector xxx. 
      
      dorothehatmap <- NRtargets[ , grepl( "logFC" , colnames( NRtargets ) ) ]
      base::rownames(dorothehatmap) <- NRtargets[,1]
      # for annotation part in heatmap
      net <- data.frame(net)
      net2 <- net[net[,3] %in% xxx,]
      selected.TFs <- dt[selected,1]
      net2 <- net2[net2[,1] %in% selected.TFs, ]
      
      for(j in selected){
        for (i in 1: nrow(dorothehatmap))  {
          net3 <- net2[net2$tf==dt[j,1],]
          target.gene <- rownames(dorothehatmap)[i]
          target.gene.tf <- net3[net3$target==target.gene,]
          
          if(!is.na(target.gene.tf[1,1])){
            
            dorothehatmap[i,dt[j,1]]=(target.gene.tf$confidence)
          }
        }
      }
      
      dorothehatmap <- data.frame(dorothehatmap)
      
      col.logfc <- grep("logFC", colnames(dorothehatmap)) # select the columns containing logFC
      
      
      dataX= dorothehatmap[,col.logfc] # only logFC columns
      colnames(dataX)=paste(rep("Data ",length(colnames(dataX))),1:length(colnames(dataX)))
      
      if(base::nrow(dataX)>1 && base::nrow(dataX)<500){
        
     
    
          HTML(paste(paste( colnames(dataX),colnames(dorothehatmap[,col.logfc]),sep =": " ),collapse  = "\n"))
        }
      
      
    })
    output$DoRothEAdownload <- downloadHandler(
      filename = function(){"DoRothEA.csv"}, 
      content = function(fname){
        write.csv(format(data.frame(dt),decimal.mark=",",big.mark="."), fname)
      }
    )
    output$dotplotfordorothea <- renderPlot({
      shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
      hyp_dots(hyp_obj)
    })
    
    output$networkDorotheplot <- renderVisNetwork({
      shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
     # if (input$nodeclr=="diffr"){validate(need(input$firstdataset != input$seconddataset,"Select different datasets!"))}
       
      req(input$networkDorothe)
      ## this part is to obtain data containing only logfc columns---------------------------------------
      data.subset_onlylogFC <- as.data.frame(data.subset[ , grepl( "logFC" , colnames( data.subset ) ) ])
 
      base::rownames(data.subset_onlylogFC) <- data.subset[,1]

      if (conditions()=="1"){
        dataX= data.subset_onlylogFC
        colnames(dataX)="Data 1"
      }else{
      col.logfc <- grep("logFC", colnames(data.subset_onlylogFC)) # select the columns containing logFC
      
      dataX= data.subset_onlylogFC[,col.logfc] # only logFC columns for all selected genes
      colnames(dataX)=paste(rep("Data ",length(colnames(dataX))),1:length(colnames(dataX))) }
      ###### --------------------------------------------------------------------------------------------
      
      dt <- data.frame(dt) # this is TF table 
      # number of nodes or TFs is determined with FDR or selected TF number
   
      if (input$fdr_tf=="FDR"){
        dt2 <- dt[dt$fdr<=netwrk$fdrcutof,]
      }else {
        if(dim(dt)[1] > as.numeric(netwrk$tfnumber)){
          dt2 <- dt[1:netwrk$tfnumber,]
        }else {  dt2 <- dt }
      }
      ###### --------------------------------------------------------------------------------------------  
      shiny::validate(need(dim(dt2)[1]>0, "Change the FDR cutoff or number of TFs!"))
      #write.csv(dt2,"dt2.csv")

      #  dt:  TF       pval    fdr     signature geneset overlap   background                hits
    #       ZNF541  7.2e-06   0.0096        30      85       4      20295            FBP2,IDNK,RGN,TKTL1
      ###### --------------------------------------------------------------------------------------------
      # this creates default node table
      if(dim(dt2)[1]>0){
        if ( input$nodesize=="fdr"){
          nodes <- data.frame(id = rownames(dt2),  label = rownames(dt2), color = "lightblue",title=dt2[,8],value=((dt2$fdr)*(-1)))        
        }else{
          nodes <- data.frame(id = rownames(dt2),  label = rownames(dt2), color = "lightblue",title=dt2[,8],value=((dt2$overlap)))        
        }
     
      if (input$nodeclr=="diffr"){
     
        if(input$firstdataset != input$seconddataset){
         newdata=data.frame(dataX[,input$firstdataset]- dataX[,input$seconddataset])
          rownames(newdata)=rownames(dataX)
          colnames(newdata)="Data diff"
          dataX=data.frame(newdata)
      
        }else {
           shiny::validate("Select different datasets!")
         # showNotification("Select different datasets")
        }
        
      }  
        
       z.score_val=c() # contains all z-score for selected dataset
       z.score_val.all=c() # contains all z-score for all dataset
      for (i in 1:nrow(dt2)){
        newVec <- unlist(sapply(as.vector( dt2[i,8]), strsplit, "\\s+", USE.NAMES = FALSE)) 
        xx <- gsub(",", " ", newVec) ; xxx <- unlist(strsplit(xx," ")) # separate the genes as char vector
        
        NRtargets <- data.frame(dataX[rownames(dataX)%in% xxx, ]) # find logFC values of selected genes.  filters the rows of a data frame df where the value in the first column (i.e., df[,1]) is present in a vector xxx. 
        colnames(NRtargets)=colnames(dataX)
        
        
        if (input$nodeclr=="diffr"){ 
          zscor=round(base::mean(NRtargets[,1])/stats::sd(NRtargets[,1]),6)
          
          z.score_val=append(z.score_val,zscor)
          }else{
          zscor=round(base::mean(NRtargets[,netwrk$nodecolor])/stats::sd(NRtargets[,netwrk$nodecolor]),6)
          
          z.score_val=append(z.score_val,zscor)
        }
         
  
        nodes$title[i]=paste(paste("z-score:",round(zscor,3),sep = " "),nodes$title[i],collapse = ",") # node's title
        nodes[i,"borderWidth"]=abs(zscor)*10                # node's borderwidth
        nodes[i,"zscor"]=(round(zscor,6))                            # zscore which will be used for coloring
        # for (j in 1:ncol(NRtargets)){
        #   zscorfor.all=round(base::mean(NRtargets[,j])/stats::sd(NRtargets[,j]),6)
        #   z.score_val.all=append(z.score_val.all,zscorfor.all)
        # }
       
        if (input$nodeclr=="diffr"){ 
          nodes[i,"wilcox_nodes_meanlogFC"]=if(nrow(NRtargets)<2) 1 else    wilcox.test(as.vector(NRtargets[,1]), mu = 0)[["p.value"]]
       
        }else{
          nodes[i,"wilcox_nodes_meanlogFC"]=if(nrow(NRtargets)<2) 1 else    wilcox.test(as.vector(NRtargets[,netwrk$nodecolor]), mu = 0)[["p.value"]]
     
        }
        
        
        nodes[i,"shape"]=if( nodes[i,"wilcox_nodes_meanlogFC"]<=0.05) "diamond" else "dot"
      }
       
       
       
       #### for colors of nodes----------
       pos_val=z.score_val[z.score_val>0]
       pos_val=pos_val[!duplicated(pos_val)]
       neg_val=z.score_val[z.score_val<=0]
       neg_val=neg_val[!duplicated(neg_val)]
       ramp1 <- colorRamp(c( "#FBEDED","#DE514C"))
       ramp2 <- colorRamp(c("#4DBFF7","#EDF8FE"))
       poscolor=rgb( ramp1(seq(0, 1, length.out = length(pos_val))), maxColorValue = 255)
       negcolor=rgb( ramp2(seq(0, 1, length.out = length(neg_val))), maxColorValue = 255)
       colortable=rbind( cbind(as.numeric(sort(pos_val)),poscolor), cbind(as.numeric(sort(neg_val)),negcolor))
       colnames(colortable)=c("val","colors")
    
       for (k in 1:nrow(nodes)){
         nodes[k,"color"]= colortable[ colortable[,1] == nodes[k,"zscor"] ,2]
       }
     
      
        
        jaccard <- function(a, b) {
          intersection = length(intersect(a, b))
          union = length(a) + length(b) - intersection
          return (intersection/union)
        }
        x=1 # indices for new column in newdt
        newdt <- data.frame()
        for (i in 1:(dim(dt2)[1]) ){
          a <- unlist(strsplit(dt2[i,8],",")) # genes are in this column
          for(k in 1:(dim(dt2)[1]) ){
            if (i !=k){
              b <-  unlist(strsplit(dt2[k,8],",")) # genes for second TF
              ovrlapgene=intersect(a,b)
              numberof.common.genes <- length(ovrlapgene) # check overlapped genes between two TFs
              jc=jaccard(a, b)
              if (numberof.common.genes>0) # if there is overlapped genes
              {
                if (dim(newdt)[1]<1){
                  newdt[x,1] <- dt2[i,1] # first TF
                  newdt[x,2] <- dt2[k,1] # other TF that has common genes with first one
                  newdt[x,3] <-numberof.common.genes # number of overlapped genes
                  newdt[x,4] <- jc # jaccard
                  newdt[x,5] <- paste(ovrlapgene,collapse = ",")  # shared genes
                  newdt[x,6] <- 100*(numberof.common.genes/length(union(a, b))) 
                 # newdt[x,7] <- paste("jc:" ,round(jc,2), "# gene:",numberof.common.genes,"\n", "% gene:", round(100*(numberof.common.genes/length(union(a, b))),1)  , sep = " " )
                  x=x+1
                }
                else{
                  selected=c(dt2[i,1],dt2[k,1])
                  nna <- rownames(newdt[(newdt[,1]==selected[2] & newdt[,2]==selected[1]),])
                  
                  if( isEmpty(nna)){
                    newdt[x,1] <- dt2[i,1]
                    newdt[x,2] <- dt2[k,1]
                    newdt[x,3] <-numberof.common.genes
                    newdt[x,4] <- jc # jaccard
                    newdt[x,5] <- paste(ovrlapgene,collapse = ",")  # shared genes
                    newdt[x,6] <- 100*(numberof.common.genes/length(union(a, b))) 
                    #newdt[x,7] <-  paste("jc:" ,round(jc,2), "# gene:",numberof.common.genes,"\n", "% gene:",  round(100*(numberof.common.genes/length(union(a, b))),1)  , sep = " " )
                    x=x+1
                  }
                }
              }
            }
          }
        }
        colnames(newdt) <- c("from","to","width","jaccard","title","percent")
        mean_val.all=c() # vector of mean of logfc for all data
        meanlogFCvector=c()# mean of logfc for selected data
        for (i in 1:nrow(newdt)){
          strng=newdt$title[i]
          newVec <- unlist(sapply(as.vector( strng), strsplit, "\\s+", USE.NAMES = FALSE)) 
          xx <- gsub(",", " ", newVec) ; xxx <- unlist(strsplit(xx," ")) # separate the genes as char vector

          NRtargets <- data.frame(dataX[rownames(dataX)%in% xxx, ] )# find logFC values of selected genes.  filters the rows of a data frame df where the value in the first column (i.e., df[,1]) is present in a vector xxx. 
          colnames(NRtargets)=colnames(dataX)
          if (input$nodeclr=="diffr") {colnum=1} else {colnum=netwrk$nodecolor}
          
         
            meanlogFC=round(base::mean(NRtargets[,colnum]),6)
         meanlogFCvector=append(meanlogFCvector,meanlogFC)
          newdt[i,"meanlog"]=meanlogFC
          newdt[i,"title"]=paste("Mean LogFC:",meanlogFC,newdt[i,"title"], sep = " ")
          
          # for (j in 1:ncol(NRtargets)){
          #   mean.all=round(base::mean(NRtargets[,j]),6)
          #   mean_val.all=append(mean_val.all,mean.all)
          # }
          newdt[i,"wilcox_meanlogFC"]= if(base::nrow(NRtargets)<2) 1 else wilcox.test(as.vector(NRtargets[,colnum]), mu = 0)[["p.value"]]
          newdt[i,"dashes"]= if( newdt[i,"wilcox_meanlogFC"]<=0.05) FALSE else TRUE
          newdt[i,"smooth"] =  FALSE
        }
        # for colors of edges -------
        pos_valedge=meanlogFCvector[meanlogFCvector>0]
        pos_valedge=pos_valedge[!duplicated(pos_valedge)]
        neg_valedge=meanlogFCvector[meanlogFCvector<=0]
        neg_valedge=neg_valedge[!duplicated(neg_valedge)]
        ramp1 <- colorRamp(c( "#FBDCDC","red"))
        ramp2 <- colorRamp(c("blue","#DCE6FE"))
        poscolor=rgb( ramp1(seq(0, 1, length.out = length(pos_valedge))), maxColorValue = 255)
        negcolor=rgb( ramp2(seq(0, 1, length.out = length(neg_valedge))), maxColorValue = 255)
        colortable=rbind( cbind(as.numeric(sort(pos_valedge)),poscolor), cbind(as.numeric(sort(neg_valedge)),negcolor))
        colnames(colortable)=c("val","colors")
        
        for (k in 1:nrow(newdt)){
          newdt[k,"color"]= colortable[ colortable[,1] == newdt[k,"meanlog"] ,2]
        }
        
        # for edges's arrows
        for(i in 1:nrow(newdt)){
          firstcol= newdt[i,"from"]
          secondcol= newdt[i,"to"]
          
          arrowline1=alldataTF%>% dplyr::filter(source==firstcol)%>% dplyr::filter(target==secondcol)
          arrowline2=alldataTF%>% dplyr::filter(source==secondcol)%>% dplyr::filter(target==firstcol)
          if(nrow(arrowline1)>0){
          # newdt[i,"arrows"]="to"  # if(arrowline1$mor <0) "R" else "A"
          #  newdt[i,"label"]=if(arrowline1$mor <0) "R" else "A"
            newdt[i,"from"] =secondcol
            newdt[i,"to"]=firstcol
            newdt[i,"arrows"]=if(arrowline1$mor <0) "to" else "from"
            newdt[i,"label"]=if(arrowline1$mor <0) "R" else "A"
          }else if (nrow(arrowline2)>0){
           newdt[i,"arrows"]=if(arrowline2$mor <0) "to" else "from"
            newdt[i,"label"]=if(arrowline2$mor <0) "R" else "A"
           
          }
          else {
            newdt[i,"arrows"]=""
            #newdt[i,"dashes"]= FALSE
           
          }
          
        }

        # Edges will be cut based on below
        if(input$cutnetwork=="jacc"){
          edges <- newdt[newdt$jaccard>= as.numeric(input$jaccardindex),]}
       else if(input$cutnetwork=="shrdgene" ){  edges <- newdt[newdt$width >= as.numeric(input$targetnumber),]}
       else if(input$cutnetwork=="edgelogfc" ){ edges <- newdt[abs(newdt$meanlog) >= as.numeric( input$edgelogfcslider),]}
       else if(input$cutnetwork=="sign"){ if(input$signifcnt=="Significant"){ edges <- newdt[newdt$wilcox_meanlogFC <= as.numeric(0.05),]} else {edges <- newdt[newdt$wilcox_meanlogFC > as.numeric(0.05),] }}
        else {edges <- newdt[newdt$percent >= as.numeric(input$percntgene),] }
    
        #if(input$showlabel==TRUE){edges }else {edges=edges[,-7]} # this will be removed later

        edges$width =scales::rescale( edges$width,to=c(2,8))
 print(edges)
        visNetwork(nodes, edges, width = "100%")%>%  
          visOptions(highlightNearest = TRUE, nodesIdSelection = T) %>% 
          visEdges(arrows = list(to = list(enabled = TRUE, type = "bar"))  , physics=FALSE, length=120)%>%
          visPhysics(stabilization = FALSE) %>%
          visLegend(enabled = TRUE,
                    #addNodes =list(
          #   list(label = "mean(logFC)<0", shape = "icon", 
          #        icon = list(code = "f068",size = 50, color = "blue")),
          #   list(label = "mean(logFC)>0", shape = "icon", 
          #        icon = list(code = "f068", size = 50, color = "red"))), 
          #   #addEdges = data.frame(label = c("Significant","Not Significant"),dashes=c(FALSE,TRUE), arrows =  "vee", color=c("red","red"),font.align = "top"), useGroups = FALSE)  %>%
          addEdges =   data.frame(enabled = TRUE,font.align = "top",label = c("mean(logFC)<0","mean(logFC)>0","Significant","Not Significant","Repressor","Activator"),dashes=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
                                  arrows = list(to = list(enabled = TRUE, scaleFactor = c(0,0,0,0,1,1),width=c(1,2,3,4,50,7) ,type = c("vee","vee","vee","vee","bar","arrow"))), color=c("blue","red","black","black","black","black")),useGroups = FALSE) %>%
          
         # addEdges =   data.frame(enabled = TRUE,font.align = "top",label = c("Significant","Not Significant","Repressor","Activator"),dashes=c(FALSE,TRUE,FALSE,FALSE),arrows = list(to = list(enabled = TRUE, scaleFactor = c(0,0,1,1), type = c("vee","vee","bar","arrow"))), color=c("black","black","black","black")),useGroups = FALSE) %>%
          addFontAwesome(version = "4.7.0") %>%
        visLayout(randomSeed = 1212)
      }
      else{
        showNotification("Change the FDR cutoff or number of TFs!")
      }
    })
  
  })
  ####################### PART 8. MSigDB  ----------------------------------

  data_reactive <- reactiveValues(
    dt =NULL,
    nn.gscat=NULL,
    minval=0
  )
  
    observeEvent(input$collect,{
      shiny::validate(
        shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!")
       )
      
    withProgress(message = 'Calculating overlapped genes...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      
      if (!is.null(input$tablefirst_rows_selected)){
        collections <- input$msigselect # selected collections
    
        a <- as.data.frame(msigdbr_collections()) # all collections
        xk <- c()
       
        for (i in 1:length(collections)) {
          if(collections[i]=="C1"|collections[i]=="C6"|collections[i]=="H"|collections[i]=="C8"){
            x <-  a[grep(collections[i],a[,1]),1:2]
            xk <- rbind(xk,x)
            #  msigdt <- msigdbr(species = "Homo sapiens", collection = xk[i,1] )
            #    dt <- rbind(dt,msigdt)
          }
          else {
            x <-  a[grep(collections[i],a[,2]),1:2]
            xk <- rbind(xk,x)
            # msigdt <- msigdbr(species = "Homo sapiens", collection = xk[i,1], subcollection =xk[i,2] )
            # dt <- rbind(dt,msigdt)
          }    }
        dt <- c();  k <- xk
      
         for (i in 1:dim(k)[1]) {
           #  msigdt <- msigdbr(species = "Homo sapiens", collection = k[i,1], subcollection =k[i,2] )
          msigdt <- msigdbr(species = "Homo sapiens", category = k[i,1], subcategory =k[i,2] )
          dt <- rbind(dt,msigdt)
        }
      
        nn.gscat = dt  %>% dplyr::distinct(gs_cat,gs_name,gene_symbol) %>% as.data.frame()
     
        data_reactive$nn.gscat <- nn.gscat
        nn = dt  %>% dplyr::distinct(gs_name, gene_symbol,gs_description) %>% as.data.frame()
    
        hh <- enricher(gene = unlist(data.subset()[1]), TERM2GENE = nn)
        if(!is.null(hh)){  
          data_reactive$dt <-  data.frame(hh@result)[,-c(2,5,6,7)]
          data_reactive$dt[,6] <- data.frame(round(data_reactive$dt[,6],3))
          data_reactive$minval <-  min(data_reactive$dt[,6])
         
        }
        
    
        output$MSigtable <- DT::renderDT({
          shiny::validate(
                        shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"),
                        shiny::need(!is.null(hh),"No overlaps found for selected genes! Please select more genes or different collections."))
          h <- data.frame(hh@result)[,-c(2,5,6,7)]
          h[,5] <- data.frame(round(h[,5],3))
          h[,4] <- data.frame(round(h[,4],3))
          h[,6] <- data.frame(round(h[,6],3))
          nn2<-unique(nn[,c(1,3)]) # unique description of gene sets
         
          h2=h %>% dplyr::left_join(nn2, by = c("ID" = "gs_name")) 
          
       
          
    
          DT::datatable( h,rownames= FALSE,extensions = c( 'Buttons', 'ColReorder'),filter = list(position = 'top'),class='stripe cell-border hover', 
                         callback=JS('$("button.buttons-copy").css("background","#cedbd3"); 
                    $("button.buttons-excel").css("background","#cedbd3"); 
                     $("button.buttons-pdf").css("background","#cedbd3"); 
                      $("button.buttons-colvis").css("background","#cedbd3"); 
                    return table;'),
                         options = list(initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#cedbd3', 'color': '164036'});","}"),
                           colReorder = TRUE,
                           lengthMenu = c(10, 30, 50,100,200,500),
                           pageLength = 10,
                           scrollX = TRUE,
                           dom = "Blfrtip",
                           buttons = list('colvis','copy', 'excel', 'pdf')
                          
                           
                           )
          )
        })
      }
      else{
        showNotification("Please select genes from Analysis tab!")
      }
      incProgress(0.5)
      shiny:: setProgress(1)
    })
    
    output$MSigtabledown <- downloadHandler(
      filename = function(){"MSigtable.csv"}, 
      content = function(fname){
        write.csv(format(data.frame(hh@result),decimal.mark=",",big.mark="."), fname)
      }
    )
    output$msigdbheatmap <- plotly::renderPlotly({
      shiny::validate(
                      shiny::need(!is.null(data.subset()),"Please select genes first!"),
                      shiny::need(!is.null(hh),"No overlaps found for selected genes! Please select more genes or different collections."),
                      need(conditions()!= "1", "To create heatmap, please select more data comparison."))
      req(input$MSigtable_rows_selected )
      h <- data.frame(hh@result)[,-c(2,5,6,7)]
     
      h[,5] <- data.frame(round(h[,5],3))
      h[,4] <- data.frame(round(h[,4],3))
      h[,6] <- data.frame(round(h[,6],3))
      df <- data.subset()
      base::rownames(df) <- df[,1]
      selected <-input$MSigtable_rows_selected # selected rows from table
      hit.genes <- as.vector(h[selected,7]) # select the genes of rows clicked
      newVec <- unlist(sapply(hit.genes, strsplit, "/", USE.NAMES = FALSE)) 
      logFCfor.selectedgenes <- df[df[,1] %in% newVec, ] # find logFC values of selected genes

      data.heatmap <- logFCfor.selectedgenes[ , grepl( "logFC" , colnames( logFCfor.selectedgenes ) ) ]
      base::rownames(data.heatmap) <- logFCfor.selectedgenes[,1]
      # for annotation part in heatmap
      selected.geneset <- h[selected,c(1,7)]
      
      for (i in 1:nrow(data.heatmap)){
        gene <- base::rownames(data.heatmap)[i] 
        for (j in 1:nrow(selected.geneset)){
          
          resultofcheck=match(1, str_detect( unlist(sapply(as.vector(selected.geneset[j,2]), strsplit, "/", USE.NAMES = FALSE)), gene))
          if (!is.na(resultofcheck)){
            data.heatmap[i,selected.geneset[j,1]]=selected.geneset[j,1]
          }
        }
      }
      
      data.heatmap <- data.frame(data.heatmap)
   
      col.logfc <- grep("logFC", colnames(data.heatmap)) # select the columns containing logFC
     
      selecteddist2 <- input$selecteddist2
      dataX= data.heatmap[,col.logfc]
      colnames(dataX)=paste(rep("Data ",length(colnames(dataX))),1:length(colnames(dataX))) # make col names like data1 data2...
      tit=paste(paste(colnames(dataX),colnames(data.heatmap[,col.logfc]),sep =": " ),collapse = "\n")
      
      if (input$usedendrogrow_msig==TRUE){ dendo=TRUE } else {dendo=FALSE}
      if (input$usedendrog_msig==TRUE){ dendoColv=TRUE } else {dendoColv=FALSE}
      # # to make midpoint 0 in color heatmap, bias is calculated
      # a=as.matrix(dataX)
      # b=length(a[a>=0]) /length(a)
      # bia=log(b)/log(0.5)
      # 
      # a1=a[a<0]
      # a2=a[a>0]
      # 
      # #write.csv(dataX,"dataX.csv")
      # colr=c(colorRampPalette(c("blue","white"),alpha=T)(length(a1)), colorRampPalette(c("white","red"),alpha=T)(length(a2)))
      # 
      # 
   
       if(base::nrow(data.heatmap)>1 && base::nrow(data.heatmap)<500){
        
        if ((ncol(data.heatmap)-length(col.logfc)) >1){
          annotpart <-  data.heatmap[,-col.logfc]
        
          abbreviatednames= abbreviate(gsub("_"," ",colnames(annotpart)), minlength = 1, use.classes = F,
                          dot = F, strict = FALSE,
                          method =  "left.kept", named = TRUE)
          colnames(annotpart)=abbreviatednames # make abreviation for annotation part in colnames
         # annotlabel=paste(paste(abbreviatednames,colnames(data.heatmap[,-col.logfc]),sep =": " ),collapse = "\n")
          heatmaply:: heatmaply(as.matrix(dataX[,col.logfc]), 
                                Colv=dendoColv,Rowv=dendo, 
                                row_side_colors =  annotpart,
                                distfun=selecteddist2,
                              
                                column_text_angle = input$col_text_Ang_msig,
                                fontsize_row = as.numeric(input$fontsize_rowc_msig),
                                fontsize_col = as.numeric(input$fontsize_colc_msig),
                               
                                key.title=tit,
                                label_names = c("Gene", "Data", "LogFC"),
                              #  key=TRUE,
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
                               # colors = colr #bluered(length(as.matrix(dataX))), plot_method="plotly",
                                )
          #%>%layout(showlegend = TRUE)
        }
        else {  
          heatmaply:: heatmaply(as.matrix(dataX[,col.logfc]), Colv=dendoColv,Rowv=dendo, distfun=selecteddist2,key.title=tit,
                                column_text_angle = input$col_text_Ang_msig,
                                fontsize_row = as.numeric(input$fontsize_rowc_msig),
                                fontsize_col = as.numeric(input$fontsize_colc_msig),
                                scale_fill_gradient_fun =ggplot2::scale_fill_gradient2(
                                  low = "blue", high = "red", 
                                  midpoint = 0
                                ) 
                                )
          
          }
        
        
   
      }else{
        showNotification(paste ("Select different terms from the above table which contains different genes!", "Heatmap cannot be drawn if there are more than 500 genes!"))
       shiny::validate("Heatmap cannot be drawn if there are more than 500 genes!")
         return(NULL)}
      
    })
    output$textmisgdb <- renderText({
      req(input$MSigtable_rows_selected )
      h <- data.frame(hh@result)[,-c(2,5,6,7)]
      h[,5] <- data.frame(round(h[,5],3))
      h[,4] <- data.frame(round(h[,4],3))
      h[,6] <- data.frame(round(h[,6],3))
      df <- data.subset()
      base::rownames(df) <- df[,1]
      selected <-input$MSigtable_rows_selected # selected rows from table
      hit.genes <- as.vector(h[selected,7]) # select the genes of rows clicked
  
      newVec <- unlist(sapply(hit.genes, strsplit, "/", USE.NAMES = FALSE)) 
      logFCfor.selectedgenes <- df[df[,1] %in% newVec, ] # find logFC values of selected genes
      
      
      data.heatmap <- logFCfor.selectedgenes[ , grepl( "logFC" , colnames( logFCfor.selectedgenes ) ) ]
      base::rownames(data.heatmap) <- logFCfor.selectedgenes[,1]
      # for annotation part in heatmap
      selected.geneset <- h[selected,c(1,7)]
      
      for (i in 1:nrow(data.heatmap)){
        gene <- base::rownames(data.heatmap)[i] 
        for (j in 1:nrow(selected.geneset)){
          
          resultofcheck=match(1, str_detect( unlist(sapply(as.vector(selected.geneset[j,2]), strsplit, "/", USE.NAMES = FALSE)), gene))
          if (!is.na(resultofcheck)){
            data.heatmap[i,selected.geneset[j,1]]=selected.geneset[j,1]
          }
        }
      }
      
      data.heatmap <- data.frame(data.heatmap)
      
      col.logfc <- grep("logFC", colnames(data.heatmap)) 
      if(base::nrow(data.heatmap)>1 && base::nrow(data.heatmap)<500){
        
        if ((ncol(data.heatmap)-length(col.logfc)) >1){
          annotpart <-  data.heatmap[,-col.logfc]
          abbreviatednames= abbreviate(gsub("_"," ",colnames(annotpart)), minlength = 1, use.classes = F,
                                       dot = F, strict = FALSE,
                                       method =  "left.kept", named = TRUE)
          HTML(paste(paste(abbreviatednames,colnames(annotpart) ,sep =": " ),collapse  = "\n"))
        }}
      
      
    })
    
  })
  observe({
    updateSelectInput(session, "select_geneset",
                      label = "Select Collections",
                    choices = unique(data_reactive$nn.gscat$gs_cat ),selected = unique(data_reactive$nn.gscat$gs_cat )[1]
    )
  })
  
  
  
  
  #observe({updateSliderInput(session,"fdrcutofmsig"," q-value",min = as.double(data_reactive$minval),max = 1,value = 1)  })
  observeEvent(input$networkplot,{ 
   
    list_enrichment <- data.frame(data_reactive$dt) # this is table created here
    nn.gscat <- data.frame( data_reactive$nn.gscat) # genes and terms
    nn.gscat2 <- nn.gscat[nn.gscat$gs_cat %in% input$select_geneset,2]

    cutof <- as.double(input$fdrcutofmsig)

    list_enrichment2 <-  list_enrichment[list_enrichment[,6] <= cutof,]
   
    list_enrichment2 <-  list_enrichment2[list_enrichment2[,1] %in% nn.gscat2,]
   print(list_enrichment2[1,4])
   
    output$msigdbNetwork <- renderVisNetwork({
      shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!"))
      termnumber=as.integer(input$termnumber)
   
    if(dim(list_enrichment2)[1]>0){
      if (nrow(list_enrichment2)>21){
        list_enrichment2=list_enrichment2[1:termnumber,]
      }
      print(dim(list_enrichment2))
      nodes <- data.frame(id = rownames(list_enrichment2),  label = rownames(list_enrichment2), color = "#0DCFDC",title= paste( paste0("p-val: ",round( list_enrichment2[,4],3)),paste0(unlist(str_split(list_enrichment2[,7],"/")),sep="",collapse = ", "),sep=". ")
 ,
                          value=( -log10(list_enrichment2[,4]) ))           
  
     x=1 # indices for new column in newdt
      newdt <- data.frame()
      for (i in 1:(dim(list_enrichment2)[1]) ){
        a <- unlist(strsplit(list_enrichment2[i,7],"/")) # genes are in this column
    
        for(k in 1:(dim(list_enrichment2)[1]) ){
          if (i !=k){
            b <-  unlist(strsplit(list_enrichment2[k,7],"/")) # genes for second TF
            numberof.common.genes <- length(intersect(a,b)) # check overlapped genes between two TFs
            sharedgenes=intersect(a,b)
         
            if (numberof.common.genes>0) # if there is overlapped genes
            {
              if (dim(newdt)[1]<1){
                newdt[x,1] <- list_enrichment2[i,1] # first TF
                newdt[x,2] <- list_enrichment2[k,1] # other TF that has common genes with first one
                newdt[x,3] <-numberof.common.genes # number of overlapped genes
                newdt[x,4] <-paste(sharedgenes, sep=" ",collapse = ", ")
                x=x+1
               
              }
              else{
                selected=c(list_enrichment2[i,1],list_enrichment2[k,1])
                nna <- base::rownames(newdt[(newdt[,1]==selected[2] & newdt[,2]==selected[1]),])
               
                if( isEmpty(nna)){
                  newdt[x,1] <- list_enrichment2[i,1]
                  newdt[x,2] <- list_enrichment2[k,1]
                  newdt[x,3] <- numberof.common.genes
                  newdt[x,4] <- paste(sharedgenes,sep=" ",collapse = ", ")
                  x=x+1
                }
              }
            }
          }
        }
      }

      colnames(newdt) <- c("from","to","width","title")
      newdt$width =scales::rescale( newdt$width,to=c(2,8))
      
      edges <- newdt
      fontsize=as.numeric(input$fontsizenet)
    
      visNetwork(nodes, edges, width = "100%")%>% visNodes(size = 'size',font=list(size=fontsize),
                                                           shadow = T,
                                                           color = 'color') %>%  visEdges( physics=FALSE,length = 350) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = T) %>% visPhysics(stabilization = FALSE) %>%
        visLayout(randomSeed = 193)
    }
    else{
      showNotification("Change the q-value cutoff!")
    }
    })
    })

  
  
  
  
  
  ####################### PART 9. Compare Cluster----------------------------------
 
   # reactiveValues
  ccreac <- reactiveValues()
  observeEvent(input$ccaction,
               {
                 ccreac$logFCthreshold <- as.numeric(input$logFCthreshold)
                 
               })
  output$kegg1 <- renderPlot({
    shiny::validate( shiny::need(!is.null(input$tablefirst_rows_selected),"Please select genes first!") )
    print(input$tabsetss)
    req(input$tabsetss == "Compare Clusters")
    req(input$ccaction)
    newresult <- data.subset()
    rownames(newresult) <- newresult[,1] # rownames are now hgnc symbols
 
    #  merged <- merge(newresult,gen, by="row.names")
    #print(colnames(merged))
    #newresult <- merged[,c(ncol(merged),3:ncol(merged))]
    selected.cutoff <-as.numeric(ccreac$logFCthreshold)
    tims=  which(grepl("logFC",colnames(newresult)))
    eg = bitr((newresult[,1]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    
    gcSample <- list()
    listindex=1
    print(tims)
    isolate({tryCatch({
      for (i in tims){
        geneList = newresult[,i]
        names(geneList) = as.character(eg$ENTREZID)
        geneList = sort(geneList, decreasing = TRUE)
        mydf <- data.frame(Entrez=names(geneList), FC=geneList)
        mydf<-na.omit(mydf)
        mydf <- mydf[abs(mydf$FC) >= selected.cutoff,]
        mydf$group <- "upregulated"
        mydf$group[mydf$FC <= (selected.cutoff*(-1)) ] <- "downregulated"
        
        gcSample[[listindex]]=(mydf[mydf$group=="upregulated",1])
        gcSample[[listindex+1]]= (mydf[mydf$group=="downregulated",1])
        listindex=listindex+2
        
        
      }
      namesof.list <- c()
      for (i in tims){
        namesof.list <- append(namesof.list, c(paste0(colnames(newresult)[i]," upregulated" ),paste0(colnames(newresult)[i]," downregulated" ," ")))
      }
      
      names(gcSample) <- namesof.list
      print(gcSample)
      ck <- compareCluster(geneCluster = gcSample, fun = "enrichGO",OrgDb='org.Hs.eg.db', pvalueCutoff=1)
      ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
      
      dotplot(ck,color="pvalue",showCategory = 30, font.size = input$font.slider, title=paste("Comparing GO enrichment results of upregulated and  downregulated genes in all comparison"))
      
    },error=function(e){
      shiny::validate("Please select different logFC threshold or more genes!!")
      return(NULL)
    })   })
    
    
   
    
  })
  
  
  ############################## DATA SUMMARY PART ------
  
  output$data_summary_ui_output <- renderUI({
    input$analyze
    if (conditions()=="1") shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"))
    if (conditions()=="12")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"))
    if (conditions()=="14")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    
    
    if (conditions()=="13")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                           need(input$control_group_3 !=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="123")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"))
    
    if (conditions()=="124")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="134")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_3!=input$treatment_group_3, "Control and Treatment samples cannot be same!"),
                                            need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))
    if (conditions()=="1234")shiny::validate(need(input$control_group!=input$treatment_group, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_2!=input$treatment_group_2, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_3!=input$treatment_group_4, "Control and Treatment samples cannot be same!"),
                                             need(input$control_group_4!=input$treatment_group_4, "Control and Treatment samples cannot be same!"))

 

    # if only First data is selected
    if (conditions()=="1"){
      
      tabsetPanel(type = "tabs",
                  
                  tabPanel(title = input$experiment_title_1,
                         
                       fluidRow(   column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                                tabPanel("Volcano Plot",
                                                                         column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                         column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                         column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                         column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                         plotOutput("volcano1",height='800px')   ),
                                                                tabPanel("Plots", fluidRow(class = "myRow2",
                                                                  column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                         fluidRow(br()
                                                                                  #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                  )),
                                                                  column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                         fluidRow(br()
                                                                                  #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                  ))
                                                                ),fluidRow(br())),
                                                                tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                         downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  )
      )}
    # only 1 and 4
    else if (conditions()=="14"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                       
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                        )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()#, downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                        ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  ),
                  
                  tabPanel(title = input$experiment_title_4,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff4", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff4", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize4", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize4", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano4",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))),
                                                          
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_4",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                        )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_4",height='600px')),
                                                                                               fluidRow(br()#, downloadButton( outputId = "downloadpca4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                        ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_4"),
                                                                      downloadButton(outputId = "raw_download_button_4", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  )
                  
                  
      )
    }
    # 1 ve 3
    else if (conditions()=="13"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                         
                                                              tabPanel("Volcano Plot",
                                                                       column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                       column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                       column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                       column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                       plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))     ),
                                                              tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                         column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                                fluidRow(br()
                                                                                                         #  downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                         )),
                                                                                         column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                                fluidRow(br()#, downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                         ))
                                                              ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  )
                  ,
                  tabPanel(title = input$experiment_title_3,
                          
                           fluidRow( column(width =12,tabsetPanel(id = "tabbox",type="tabs",
                                                          
                                                              tabPanel("Volcano Plot",
                                                                       column(width=3, sliderInput("pCutoff3", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                       column(width=3,sliderInput("FCcutoff3", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                       column(width=3,sliderInput("pointsize3", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                       column(width=3,sliderInput("labsize3", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                       plotOutput("volcano3",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                              tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                         column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_3",height='600px')),
                                                                                                fluidRow(br()
                                                                                                         # downloadButton( outputId = "downloadbar3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                         )),
                                                                                         column(6, fluidRow(plotOutput(outputId = "pca_plot_output_3",height='600px')),
                                                                                                fluidRow(br()#, downloadButton( outputId = "downloadpca3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                         ))
                                                              ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_3"),
                                                                      downloadButton(outputId = "raw_download_button_3", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  )
                  
      )
    }
    # 1 ve 2
    else if (conditions()=="12"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                         
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))      ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  ),
                  tabPanel(title = input$experiment_title_2,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff2", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff2", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize2", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize2", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano2",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_2"),
                                                                      downloadButton(outputId = "raw_download_button_2", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  )
                  
      )
    }
    # 1,2,3
    else if (conditions()=="123"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                      
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))      ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  ),
                  tabPanel(title = input$experiment_title_2,
                          
                         
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff2", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff2", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize2", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize2", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano2",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_2"),
                                                                      downloadButton(outputId = "raw_download_button_2", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  ), tabPanel(title = input$experiment_title_3,
                            
                              fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                                tabPanel("Volcano Plot",
                                                                         column(width=3, sliderInput("pCutoff3", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                         column(width=3,sliderInput("FCcutoff3", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                         column(width=3,sliderInput("pointsize3", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                         column(width=3,sliderInput("labsize3", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                         plotOutput("volcano3",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                                tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                           column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_3",height='600px')),
                                                                                                  fluidRow(br()
                                                                                                           # downloadButton( outputId = "downloadbar3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                  )),
                                                                                           column(6, fluidRow(plotOutput(outputId = "pca_plot_output_3",height='600px')),
                                                                                                  fluidRow(br()#, downloadButton( outputId = "downloadpca3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                  ))
                                                                ),fluidRow(br())),
                                                                tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_3"),
                                                                         downloadButton(outputId = "raw_download_button_3", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                              ) ))
                  )
                  
      )
    }
    # 1,3,4
    else if (conditions()=="134"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                        
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))      ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                  ),
                  tabPanel(title = input$experiment_title_3,
                       
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff3", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff3", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize3", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize3", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano3",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_3",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        # downloadButton( outputId = "downloadbar3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_3",height='600px')),
                                                                                               fluidRow(br()#, downloadButton( outputId = "downloadpca3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_3"),
                                                                      downloadButton(outputId = "raw_download_button_3", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  ),tabPanel(title = input$experiment_title_4,
                       
                             fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                               tabPanel("Volcano Plot",
                                                                        column(width=3, sliderInput("pCutoff4", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                        column(width=3,sliderInput("FCcutoff4", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                        column(width=3,sliderInput("pointsize4", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                        column(width=3,sliderInput("labsize4", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                        plotOutput("volcano4",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                               tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                          column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadbar4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 )),
                                                                                          column(6, fluidRow(plotOutput(outputId = "pca_plot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadpca4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 ))
                                                               ),fluidRow(br())),
                                                               tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_4"),
                                                                        downloadButton(outputId = "raw_download_button_4", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  )
                  
      )
    }
    # 1,2,4
    else if (conditions()=="124"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))      ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                                                             
                           )))
                  ),
                  tabPanel(title = input$experiment_title_2,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff2", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff2", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize2", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize2", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano2",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_2"),
                                                                      downloadButton(outputId = "raw_download_button_2", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                           ) ))
                           
                           
                  ),tabPanel(title = input$experiment_title_4,
                         
                             fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                               tabPanel("Volcano Plot",
                                                                        column(width=3, sliderInput("pCutoff4", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                        column(width=3,sliderInput("FCcutoff4", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                        column(width=3,sliderInput("pointsize4", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                        column(width=3,sliderInput("labsize4", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                        plotOutput("volcano4",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                               tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                          column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadbar4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 )),
                                                                                          column(6, fluidRow(plotOutput(outputId = "pca_plot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadpca4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 ))
                                                               ),fluidRow(br())),
                                                               tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_4"),
                                                                        downloadButton(outputId = "raw_download_button_4", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                             ) ))
                  )
      )
    }
    # 1,2,3,4
    else if (conditions()=="1234"){
      tabsetPanel(type = "tabs",
                  tabPanel(title = input$experiment_title_1,
                        
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano1",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))     ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_1",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca1",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_1"),
                                                                      downloadButton(outputId = "raw_download_button_1", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                                                             
                           ) ))
                  ),
                  tabPanel(title = input$experiment_title_2,
                         
                           fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                             tabPanel("Volcano Plot",
                                                                      column(width=3, sliderInput("pCutoff2", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                      column(width=3,sliderInput("FCcutoff2", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                      column(width=3,sliderInput("pointsize2", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                      column(width=3,sliderInput("labsize2", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                      plotOutput("volcano2",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                             tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                        column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadbar2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               )),
                                                                                        column(6, fluidRow(plotOutput(outputId = "pca_plot_output_2",height='600px')),
                                                                                               fluidRow(br()
                                                                                                        #downloadButton( outputId = "downloadpca2",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                               ))
                                                             ),fluidRow(br())),
                                                             tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_2"),
                                                                      downloadButton(outputId = "raw_download_button_2", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                                                             
                           ) ))
                  ), tabPanel(title = input$experiment_title_3,
                          
                              fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                                tabPanel("Volcano Plot",
                                                                         column(width=3, sliderInput("pCutoff3", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                         column(width=3,sliderInput("FCcutoff3", label = h4("FCcutoff"), min = 0,  max =5, value =1, step = 0.1)),
                                                                         column(width=3,sliderInput("pointsize3", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                         column(width=3,sliderInput("labsize3", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                         plotOutput("volcano3",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                                tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                           column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_3",height='600px')),
                                                                                                  fluidRow(br()
                                                                                                           # downloadButton( outputId = "downloadbar3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                  )),
                                                                                           column(6, fluidRow(plotOutput(outputId = "pca_plot_output_3",height='600px')),
                                                                                                  fluidRow(br()#, downloadButton( outputId = "downloadpca3",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                  ))
                                                                ),fluidRow(br())),
                                                                tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_3"),
                                                                         downloadButton(outputId = "raw_download_button_3", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))) ))
                  ),tabPanel(title = input$experiment_title_4,
                        
                             fluidRow(column(width = 12,tabsetPanel(id = "tabbox",type="tabs",
                                                               tabPanel("Volcano Plot",
                                                                        column(width=3, sliderInput("pCutoff4", label = h4("pCutoff"), min = 0.0001,  max =0.1, value = 0.05)),
                                                                        column(width=3,sliderInput("FCcutoff4", label = h4("FCcutoff"), min = 0,   max =5, value =1, step = 0.1)),
                                                                        column(width=3,sliderInput("pointsize4", label =h4("PointSize"), min = 1,  max =10, value = 3)),
                                                                        column(width=3,sliderInput("labsize4", label = h4("LabSize"), min = 1,  max =10, value = 3)),
                                                                        plotOutput("volcano4",height='800px')%>% withSpinner(type = getOption("spinner.type", 6), color = getOption("spinner.color", default = "#00838f"), size = getOption("spinner.size", default = 1))   ),
                                                               tabPanel("Plots", fluidRow(class = "myRow2",
                                                                                          column(6, fluidRow(plotOutput(outputId = "library_size_barplot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadbar4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 )),
                                                                                          column(6, fluidRow(plotOutput(outputId = "pca_plot_output_4",height='600px')),
                                                                                                 fluidRow(br()
                                                                                                          #downloadButton( outputId = "downloadpca4",label = "Download",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #123F34;")
                                                                                                 ))
                                                               ),fluidRow(br())),
                                                               tabPanel("Raw Data", dataTableOutput(outputId = "raw_dataTable_output_4"),
                                                                        downloadButton(outputId = "raw_download_button_4", label = "Download raw data",class = paste("btn btn-default shiny-download-link",'mybutton'), style = "color: #17414E;"))
                                                               
                             ) ))
                  )
      )
    }
    
  })
  # volcano plots
  output$volcano1 <- renderPlot({
     result <- limmaall( experiment_title = reac$exp1, control_group = input$control_group, treatment_group = input$treatment_group)
     result <- cbind(rownames(result),result)
     colnames(result) <- c("hgnc_symbol","logFC", "AveExpr", "t value", "p value", "p adj value", "B value")
   
    EnhancedVolcano(result,
                    lab = rownames(result),
                    x = "logFC" ,
                    y =  "p value", 
                    pCutoff = input$pCutoff,
                    FCcutoff =input$FCcutoff ,
                    pointSize = input$pointsize,
                    labSize = input$labsize)
  })
  output$volcano2 <- renderPlot({
    result <- limmaall( experiment_title = reac$exp2, control_group = input$control_group_2, treatment_group = input$treatment_group_2)
    
    result <- cbind(rownames(result),result)
    colnames(result) <- c("hgnc_symbol","logFC", "AveExpr", "t value", "p value", "p adj value", "B value")
    EnhancedVolcano(result,
                    lab = rownames(result),
                    x = "logFC" ,
                    y =  "p value", 
                    pCutoff = input$pCutoff2,
                    FCcutoff =input$FCcutoff2 ,
                    pointSize = input$pointsize2,
                    labSize = input$labsize2)
  })
  output$volcano3 <- renderPlot({
    result <- limmaall( experiment_title = reac$exp3, control_group = input$control_group_3, treatment_group = input$treatment_group_3)
    
    result <- cbind(rownames(result),result)
    colnames(result) <- c("hgnc_symbol","logFC", "AveExpr", "t value", "p value", "p adj value", "B value")
    EnhancedVolcano(result,
                    lab = rownames(result),
                    x = "logFC" ,
                    y =  "p value", 
                    pCutoff = input$pCutoff3,
                    FCcutoff =input$FCcutoff3 ,
                    pointSize = input$pointsize3,
                    labSize = input$labsize3)
  })
  output$volcano4 <- renderPlot({
    result <-  limmaall( experiment_title = reac$exp4, control_group = input$control_group_4, treatment_group = input$treatment_group_4)
    result <- cbind(rownames(result),result)
    colnames(result) <- c("hgnc_symbol","logFC", "AveExpr", "t value", "p value", "p adj value", "B value")
    EnhancedVolcano(result,
                    lab = rownames(result),
                    x = "logFC" ,
                    y =  "p value", 
                    pCutoff = input$pCutoff4,
                    FCcutoff =input$FCcutoff4 ,
                    pointSize = input$pointsize4,
                    labSize = input$labsize4)
  })
  
  
  

  
  ###- ---------- -------------- --------------  Raw data, library size plots, MDS in "DATA SUMMARY" tab-----------------------

 
  # raw data object that is used for the table output in "Data summary" tab
  # if there are two or more comparisons, and they are from different experiments, merge the two datasets and show them in the table
  raw_data_1 <- eventReactive(input$analyze, {
    raw_data <- rawdatafunction(input$experiment_title_1)
  })
  raw_data_2 <- eventReactive(input$analyze, {
    raw_data <- rawdatafunction(input$experiment_title_2)
  })
  raw_data_3 <- eventReactive(input$analyze, {
    raw_data <- rawdatafunction(input$experiment_title_3)
  })
  raw_data_4 <- eventReactive(input$analyze, {
    raw_data <- rawdatafunction(input$experiment_title_4)
  })
  # library_size_barplot
  b_plot1<- eventReactive(input$analyze, {
    if (input$experiment_title_1 %in% onlyrnaseq){
      count_data <- raw_data_1()
      library_size <- as.vector(colSums(count_data[, -c(1:2)]))
      par(mar=c(10,4,4,4))
      barplot(library_size, names = colnames(count_data[, -c(1:2)]), las = 2, col = "purple", main = "Library Size Bar Plot")
    }else {
      par(mar=c(10,4,4,4))
      boxplot(raw_data_1()[,-1],las = 2, boxwex=0.6, notch=T, outline=FALSE, main="Box-and-whisker plot for selected array",col = "purple")
    }
  })
  output$library_size_barplot_output_1 <- renderPlot({input$analyze
    b_plot1()
   
  })
  output$library_size_barplot_output_2 <- renderPlot({
    if (input$experiment_title_2 %in% onlyrnaseq){
      count_data <- raw_data_2()
      library_size <- as.vector(colSums(count_data[, -c(1:2)]))
      par(mar=c(10,4,4,4))
      barplot(library_size, names = colnames(count_data[, -c(1:2)]), las = 2, col = "purple", main = "Library Size Bar Plot")
    }else {
      par(mar=c(10,4,4,4))
      boxplot(raw_data_2()[,-1],las = 2, boxwex=0.6, notch=T, outline=FALSE, main="Box-and-whisker plot for selected array",col = "purple")
    }
    
  })
  output$library_size_barplot_output_3 <- renderPlot({
    if (input$experiment_title_3 %in% onlyrnaseq){
      count_data <- raw_data_3()
      library_size <- as.vector(colSums(count_data[, -c(1:2)]))
      par(mar=c(10,4,4,4))
      barplot(library_size, names = colnames(count_data[, -c(1:2)]), las = 2, col = "purple", main = "Library Size Bar Plot")
    }else {
      par(mar=c(10,4,4,4))
      boxplot(raw_data_3()[,-1],las = 2, boxwex=0.6, notch=T, outline=FALSE, main="Box-and-whisker plot for selected array",col = "purple")
    }
  })
  output$library_size_barplot_output_4 <- renderPlot({
    
    if (input$experiment_title_4 %in% onlyrnaseq){
      count_data <- raw_data_4()
      library_size <- as.vector(colSums(count_data[, -c(1:2)]))
      par(mar=c(10,4,4,4))
      barplot(library_size, names = colnames(count_data[, -c(1:2)]), las = 2, col = "purple", main = "Library Size Bar Plot")
    }else {
      par(mar=c(10,4,4,4))
      boxplot(raw_data_4()[,-1],las = 2, boxwex=0.6, notch=T, outline=FALSE, main="Box-and-whisker plot for selected array",col = "purple")
    }
  })
  
  output$pca_plot_output_1 <- renderPlot({
   req( input$analyze)
  
    if (input$experiment_title_1 %in% onlyrnaseq){
      rawdata <- raw_data_1()[,  -c(1:2)]
      design_grouping=data.frame(row.names = colnames(rawdata), condition=colnames(rawdata))
      
      dds=DESeqDataSetFromMatrix(rawdata, colData = design_grouping, design = ~condition)
      # Take regularized logarithm
      rld=rlog(dds, blind = TRUE)
      pcaData <- plotPCA(rld, returnData = TRUE) 
      ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
        geom_point(size =4) +
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }else {
      rawdata <- raw_data_1()[,  -1]
      PCA_g <- prcomp(t(rawdata))
      pcadata<-as.data.frame(PCA_g$x)
      pcadata[,"name"] <- colnames(rawdata)
      samples=factor(colnames(rawdata))
      ggplot(as.data.frame(pcadata), aes(x = PC1, y = PC2, color = samples)) + 
        geom_point(size =4)+
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }
    
    
    
  })
  output$pca_plot_output_2 <- renderPlot({req(input$analyze)
    

    if (input$experiment_title_2 %in% onlyrnaseq){
      rawdata <- raw_data_2()[,  -c(1:2)]
      design_grouping=data.frame(row.names = colnames(rawdata), condition=colnames(rawdata))
      
      dds=DESeqDataSetFromMatrix(rawdata, colData = design_grouping, design = ~condition)
      # Take regularized logarithm
      rld=rlog(dds, blind = TRUE)
      pcaData <- plotPCA(rld, returnData = TRUE) 
      ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
        geom_point(size =4) +
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }else {
      rawdata <- raw_data_2()[,  -1]
      PCA_g <- prcomp(t(rawdata))
      pcadata<-as.data.frame(PCA_g$x)
      pcadata[,"name"] <- colnames(rawdata)
      samples=factor(colnames(rawdata))
      ggplot(as.data.frame(pcadata), aes(x = PC1, y = PC2, color = samples)) + 
        geom_point(size =4)+
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }
    

  })
  output$pca_plot_output_3 <- renderPlot({req(input$analyze)

    if (input$experiment_title_3 %in% onlyrnaseq){
      rawdata <- raw_data_3()[,  -c(1:2)]
      design_grouping=data.frame(row.names = colnames(rawdata), condition=colnames(rawdata))
      
      dds=DESeqDataSetFromMatrix(rawdata, colData = design_grouping, design = ~condition)
      # Take regularized logarithm
      rld=rlog(dds, blind = TRUE)
      pcaData <- plotPCA(rld, returnData = TRUE) 
      ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
        geom_point(size =4) +
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }else {
      rawdata <- raw_data_3()[,  -1]
      PCA_g <- prcomp(t(rawdata))
      pcadata<-as.data.frame(PCA_g$x)
      pcadata[,"name"] <- colnames(rawdata)
      samples=factor(colnames(rawdata))
      ggplot(as.data.frame(pcadata), aes(x = PC1, y = PC2, color = samples)) + 
        geom_point(size =4)+
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }
    
  })
  output$pca_plot_output_4 <- renderPlot({req(input$analyze)
 
    if (input$experiment_title_4 %in% onlyrnaseq){
      rawdata <- raw_data_4()[,  -c(1:2)]
      design_grouping=data.frame(row.names = colnames(rawdata), condition=colnames(rawdata))
      
      dds=DESeqDataSetFromMatrix(rawdata, colData = design_grouping, design = ~condition)
      # Take regularized logarithm
      rld=rlog(dds, blind = TRUE)
      pcaData <- plotPCA(rld, returnData = TRUE) 
      ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(condition))) + 
        geom_point(size =4) +
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }else {
      rawdata <- raw_data_4()[,  -1]
      PCA_g <- prcomp(t(rawdata))
      pcadata<-as.data.frame(PCA_g$x)
      pcadata[,"name"] <- colnames(rawdata)
      samples=factor(colnames(rawdata))
      ggplot(as.data.frame(pcadata), aes(x = PC1, y = PC2, color = samples)) + 
        geom_point(size =4)+
        geom_label_repel(aes(label = name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') 
    }
    
  })
  
 
  ###- ---------- --------------  OUTPUT in "DATA SUMMARY" in TEXT and TABLES -----------------------
  
  output$raw_dataTable_output_1 <- DT::renderDataTable(raw_data_1(),caption = paste("Raw data from the ", input$experiment_title_1, ".", sep = ""), 
                                                       options = list(scrollX = TRUE,initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color':'#cedbd3', 'color': '164036'});",
                                                         "}")),rownames = FALSE)
  output$raw_download_button_1 <- shiny::downloadHandler(
    filename = function() {
      paste(input$experiment_title, " raw data", '.csv', sep='')
    },
    content = function(con) {
      write.csv(raw_data_1()[input$raw_dataTable_output_1_rows_all, ], con)
    }
  )
  
  
  output$raw_dataTable_output_2 <- DT::renderDataTable(raw_data_2(),caption=paste("Raw data from the ", input$experiment_title_2, ".", sep = ""), 
                                                       options = list(scrollX = TRUE,initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color':'#cedbd3', 'color': '164036'});",
                                                         "}")),rownames = FALSE)
  output$raw_download_button_2 <- shiny::downloadHandler(
    filename = function() {
      paste(input$experiment_title_2, " raw data", '.csv', sep='')
    },
    content = function(con) {
      write.csv(raw_data_2()[input$raw_dataTable_output_2_rows_all, ], con)
    }
  )
  
  output$raw_dataTable_output_3 <- DT::renderDataTable(raw_data_3(), caption=paste("Raw data from the ", input$experiment_title_3, ".", sep = ""),
                                                       options = list(scrollX = TRUE,initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color':'#cedbd3', 'color': '164036'});",
                                                         "}")),rownames = FALSE)
  output$raw_download_button_3 <- shiny::downloadHandler(
    filename = function() {
      paste(input$experiment_title_3, " raw data", '.csv', sep='')
    },
    content = function(con) {
      write.csv(raw_data_3()[input$raw_dataTable_output_3_rows_all, ], con)
    }
  )
  
  output$raw_dataTable_output_4 <- DT::renderDataTable(raw_data_4(), caption=paste("Raw data from the ", input$experiment_title_4, ".", sep = ""),
                                                       options = list(scrollX = TRUE,initComplete = JS(
                                                         "function(settings, json) {",
                                                         "$(this.api().table().header()).css({'background-color':'#cedbd3', 'color': '164036'});",
                                                         "}")),rownames = FALSE)
  output$raw_download_button_4 <- shiny::downloadHandler(
    filename = function() {
      paste(input$experiment_title_4, " raw data", '.csv', sep='')
    },
    content = function(con) {
      write.csv(raw_data_3()[input$raw_dataTable_output_4_rows_all, ], con)
    }
  )
  
  
}
 
