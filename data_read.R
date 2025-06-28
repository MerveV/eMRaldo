# exp selection is shorten with numbers
# conditions <- reactive({
#   req(input$analyze)
#   if (reac$exp2== "Select"&&reac$exp3 == "Select"&& reac$exp4 != "Select"){ 
#     conditions <- "14"
#   }
#   else if (reac$exp2 == "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
#     conditions <- "13"
#   }    # 1 ve 2
#   else if(reac$exp3 == "Select"&& reac$exp4 == "Select"&& reac$exp2 != "Select"){
#     conditions <- "12"
#   }
#   else if(reac$exp2 != "Select"&& reac$exp4 == "Select"&& reac$exp3 != "Select"){
#     conditions <- "123"
#   } # 1,3,4
#   else if(reac$exp3 != "Select"&& reac$exp2 == "Select"&& reac$exp4 != "Select"){
#     conditions <- "134"
#   }
#   else if(reac$exp2 != "Select"&& reac$exp3 == "Select"&& reac$exp4 != "Select"){
#     conditions <- "124"
#   }
#   else if(reac$exp2 != "Select"&& reac$exp3 != "Select"&& reac$exp4 != "Select"){
#     conditions <- "1234"
#   }
#   else if(reac$exp2 == "Select"&& reac$exp3 == "Select"&& (reac$exp4 == "Select"|| reac$userfile==TRUE)){
#     conditions <- "1"
#   }
# })


# reading files in global
## -----------READING DATASETS -----------

ST.rawdata <- read.csv("data_files/ST_Rawdata.csv", header=TRUE)
ST_allfiltered <- read.csv("data_files/ST_MCF7_ZR751_Filtered.csv", header=TRUE)
ST_onlyMCF.filtered <- read.csv("data_files/ST_only_MCF7_Filtered_sep.csv", header=TRUE)
ST_onlyZR75.filtered <- read.csv("data_files/ST_only_ZR751_Filtered_sep.csv", header=TRUE)
GSE172478.rawdata <- read.csv("data_files/GSE172478_entrez_hgnc_rawdata.csv", header=TRUE)
GSE172478.filtered <- read.csv("data_files/GSE172478_Filtered.csv", header=TRUE)

gene.list <- read.csv("data_files/genename_ensembleID_63682genes.csv",header = T)
rownames(gene.list) <- gene.list$genename
#GSE151321 <- read.csv("data_files/GSE151321_kidney_Filtered.csv", header=TRUE)



##### Analysis of RNAseq 

observeEvent(input$analyze,
             {
               reac$analyze_type <- input$analyze_type
               
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
ST.group.MCF7=list( "MCF7_EV"  ,   "MCF7_EVA",    "MCF7_EVAS" ,  "MCF7_OV"   ,  "MCF7_OVA" ,"MCF7_OVAS")
ST.group.ZR75=list( "ZR751_EV",  "ZR751_EVA" ,"ZR751_OV",  "ZR751_OVA")

GSE172478_groups=list("Etoh","Methanol","Aldo","Cortisol","Cortisol_RU486")
#GSE151321_kidneyl=list("Control","Aldo","Finerenone","Spironolactone","Aldo_Spiro","Aldo_Fine")

switc <- function(exp){
  group_names <- switch(exp,
                        "ST_MCF7_ZR751"=ST.Group.all,
                        "ST_MCF7"=ST.group.MCF7,
                        "ST_ZR751"=ST.group.ZR75,
                        "GSE172478_RPEcells"=GSE172478_groups
                        
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
  if (reac$exp2== "Select"&&reac$exp3 == "Select"){ 
    conditions <- "1"
  }
  else if (reac$exp2 == "Select"&& reac$exp3 != "Select"){
    conditions <- "13"
  }    # 1 ve 2
  else if(reac$exp3 == "Select"&&  reac$exp2 != "Select"){
    conditions <- "12"
  }
  else if(reac$exp2 != "Select"&& reac$exp3 != "Select"){
    conditions <- "123"
  } # 1,3,4
  
  
})

uploaded.data <- reactive({
  req(input$analyze)
  if(!is.null(input$userfile)){
    print("girdi")
    tryCatch(
      {
        uploaded.data <- read.csv(input$userfile$datapath,header = TRUE)
        colnames(uploaded.data)[1]="Gene"
        print(head(uploaded.data))
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    if(base::ncol(uploaded.data)<5){
      return(uploaded.data)
    }else{shiny::validate("Please make sure that your data contains 4 columns")}
    mergedfiles=merge(uploaded.data,gene.list,by.x="Gene",by.y="genename")
    if(base::nrow(mergedfiles)>0){
      return(uploaded.data)
    }else{shiny::validate("Please make sure that first column contains gene names!")}
  }else {shiny::validate("Please upload your data containing gene names in the first colum, logFC in the second column, p-value a in the third colum and p-adj value in fourth column.")}
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

# create limma or DeSeq analyzed tables here
observeEvent(input$analyze, {
  
  if (reac$analyze_type=="Limma-Voom"){ columnames <- input$colname } else { columnames <- input$colname2 }
  
  
  if(reac$analyze_type=="Limma-Voom"){
    withProgress(message = 'Analyzing the selected data with Limma-Voom and merging them...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      print(conditions())
      if (conditions()=="1"){
        incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom...', detail = "First data")
        l<- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.6)
        limma1  <- cbind(rownames(  l ), l  )
        rownames(limma1) <- 1:dim(limma1  )[1]
        colnames(limma1 )[1] <- "hgnc_symbol"
        
        if(reac$uploadoption==TRUE){
          
          uploaded.data=uploaded.data()
          merged=merge(limma1, uploaded.data,by.x= "hgnc_symbol", by.y="Gene") # merge uploaded data and limma result
          reac$limma<-  metaanalysis(data.frame(merged),reac$metod,columnames ) # apply metaanalysis
          
        }else if(reac$uploadoption==FALSE && reac$exp4 != "Select"){
          
          limma4 <- limmaall( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
          limma4 <-cbind(rownames(limma4),limma4) 
          colnames(limma4)[1] <- "hgnc_symbol"
          reac$limma <-base::merge(limma1,limma4, by="hgnc_symbol")
          
          colnames(reac$limma)[1] <- "hgnc_symbol"
          reac$limma<-  metaanalysis(data.frame(reac$limma),reac$metod,columnames )
          print(head( reac$limma))
        }
        else{
          print("sadece 1")
          reac$limma <- limma1
          print(head(reac$limma))}
        
        
      }
      else if (conditions()=="12"){ 
        incProgress(0.2, message = 'Analyzing the selected data with Limma-Voom  and merging them...', detail = "First data")
        limma1 <- limmaall( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with Limma-Voom and merging them...', detail = "Second data")
        limma2 <- limmaall( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
        reac$limma <-base::merge(limma1,limma2, by="row.names")
        colnames(reac$limma)[1] <- "hgnc_symbol"
        if(reac$uploadoption==TRUE){
          uploaded.data=uploaded.data()
          reac$limma=merge(reac$limma, uploaded.data,by.x= "hgnc_symbol", by.y="Gene")
        }
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
  }else {
    withProgress(message = 'Analyzing the selected data with DESeq2 and merging them...', style = "notification", value = 0.1, {
      Sys.sleep(0.25)
      
      if (conditions()=="1"){
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2...', detail = "First data")
        l<- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.6)
        reac$limma  <- cbind(rownames(  l ), l  )
        rownames(reac$limma ) <- 1:dim(reac$limma  )[1]
        colnames(reac$limma )[1] <- "hgnc_symbol"
      }
      else if (conditions()=="12"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma2 <- Deseq2( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
        reac$limma <-base::merge(limma1,limma2, by="row.names")
        colnames(reac$limma)[1] <- "hgnc_symbol"
      }
      else if (conditions()=="14"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma4 <- Deseq2( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
        reac$limma <-base::merge(limma1,limma4, by="row.names")
        colnames(reac$limma)[1] <- "hgnc_symbol"
      }
      else if (conditions()=="13"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma3 <- Deseq2( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
        incProgress(0.7, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        reac$limma <-base::merge(limma1,limma3, by="row.names")
        colnames(reac$limma)[1] <- "hgnc_symbol"
      }
      else if (conditions()=="123"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma2 <- Deseq2( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
        incProgress(0.7, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        limma3 <- Deseq2( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
        incProgress(0.9, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        l12  <- merge(limma1,limma2, by="row.names")
        colnames(  l12)[1] <- "hgnc_symbol"
        limma3 <-cbind(rownames(limma3),limma3) 
        colnames(limma3)[1] <- "hgnc_symbol"
        reac$limma <- merge(l12 ,limma3, by="hgnc_symbol")
      }
      else if (conditions()=="124"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma2 <- Deseq2( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
        incProgress(0.7, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        limma4 <- Deseq2( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
        incProgress(0.9, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        l12  <- merge(limma1,limma2, by="row.names")
        colnames(  l12)[1] <- "hgnc_symbol"
        limma4 <-cbind(rownames(limma4),limma4) 
        colnames(limma4)[1] <- "hgnc_symbol"
        reac$limma <- merge(l12 ,limma4, by="hgnc_symbol")
      }
      else if (conditions()=="134"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma2 <- Deseq2( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
        incProgress(0.7, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        limma4 <- Deseq2( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
        incProgress(0.9, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        l12  <- merge(limma1,limma2, by="row.names")
        colnames(  l12)[1] <- "hgnc_symbol"
        limma4 <-cbind(rownames(limma4),limma4) 
        colnames(limma4)[1] <- "hgnc_symbol"
        reac$limma <- merge(l12 ,limma4, by="hgnc_symbol")
      }
      else if (conditions()=="1234"){ 
        incProgress(0.2, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "First data")
        limma1 <- Deseq2( experiment_title = reac$exp1, control_group = reac$control_group, treatment_group = reac$treatment_group)
        incProgress(0.4, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Second data")
        limma2 <- Deseq2( experiment_title = reac$exp2, control_group = reac$control_group_2, treatment_group = reac$treatment_group_2)
        incProgress(0.5, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Third data")
        limma3 <- Deseq2( experiment_title = reac$exp3, control_group = reac$control_group_3, treatment_group = reac$treatment_group_3)
        incProgress(0.7, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Fourth data")
        limma4 <- Deseq2( experiment_title = reac$exp4, control_group = reac$control_group_4, treatment_group = reac$treatment_group_4)
        incProgress(0.9, message = 'Analyzing the selected data with DESeq2 and merging them...', detail = "Fourth data")
        l12<- merge(limma1,limma2, by="row.names")
        colnames(l12 )[1] <- "hgnc_symbol"
        l3 <-cbind(rownames(limma3),limma3) 
        colnames(l3)[1] <- "hgnc_symbol"
        l123  <- merge(  l12 ,l3, by="hgnc_symbol")
        l4<-cbind(rownames(limma4),limma4) 
        colnames(l4)[1] <- "hgnc_symbol"
        reac$limma  <- merge(  l123 ,l4, by="hgnc_symbol")
      }
      
      shiny::setProgress(1)
    })
  }
  
  
})