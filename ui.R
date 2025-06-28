#options(repos = BiocManager::repositories())

#install.packages("remotes")
#remotes::install_github("AnalytixWare/ShinySky")

#setwd("~/Desktop/HRonBC")

# GSE196876 24 hours
# GSE184378 1 hour
# GSE99626 24 hour 
#GSE68358 3 hrs


# RNAseq data
l <- list("ST_MCF7_ZR751","GSE172478_RPE","GSE151321_HK","GSE113659_PAE","GSE160497_MCF7","GSE28009_EA","GSE84992_SM","GSE70822_PM")
# "ST_ZR751"



ui <-fluidPage(  #  setSliderColor(c(rep("#00838f",600)), c(1:600)),
  shinyWidgets::chooseSliderSkin("Flat", color = "#00838f"),
  use_theme(create_theme(
    theme = "default",
    bs_vars_global(
      body_bg = "#ffffff",
      text_color = "#17414E"
      
    ),
    
    bs_vars_navbar(
      padding_horizontal = "15px",
      default_bg = "#0CA7B2", 
      default_color = "#00838f",
      default_link_color = "#FFFFFF", 
      default_link_active_color = "#FFFFFF",
      default_border ="#FFFFFF"
    ),
    bs_vars_color(
      gray_base = "#17414E", # 
      gray_darker = "#17414E",
      gray_dark = "#17414E",
      gray = "#cedbd3",
      gray_light = "#17414E",# tablonun caption text rengi 
      gray_lighter = "#cedbd3", # tab i hover edince
      brand_primary = "#17414E",# text color in tab
      brand_success = "#D95F02", #
      brand_info = "#233F29", # 
      brand_warning = "#FCAA1C", # 
      brand_danger = "#2b6777" # 
      
    ),
    bs_vars_panel(
      border_radius = "0px",
      default_text = "#FFF",
      default_heading_bg = "#FFF",
      default_border = "#FFF",
      primary_heading_bg = "#2b6777", # functional profile panel colors, heading
      primary_border = "#2b6777",
      success_heading_bg = "#D95F02",
      success_border = "#D95F02",
      success_text = "#FFF",
      danger_heading_bg = "#7570B3",
      danger_border = "#7570B3",
      danger_text = "#FFF",
      info_text = "#FFF",
      info_border = "#2b6777",
      info_heading_bg = "#2b6777"
    ),
    bs_vars_state(
      success_text = "#FFF",
      success_bg = "#28a745",
      success_border = "#28a745",
      info_text = "#FFF",
      info_bg = "#233F29",
      info_border = "#233F29", 
      danger_text = "#FFF",
      danger_bg = "#00838f",
      danger_border = "#00838f" ,
      warning_text = NULL,
      warning_bg = NULL,
      warning_border = NULL
    ),
    # bs_vars_badge(
    #   color = "firebrick",
    #   bg = "steelblue"
    # ),
    # 164036 koyu yeşil text rengi fln
    bs_vars_tabs(
      border_color = "#17414E", # 
      active_link_hover_bg = "#cedbd3",
      active_link_hover_color = "#17414E", # 
      active_link_hover_border_color = "#17414E",
      link_hover_border_color = "#17414E",
      justified_link_border_color = "#7570B3",
      justified_active_link_border_color = "#7570B3"
    ),
    bs_vars_progress(
      bg ="15px",
      bar_color = "#00838f",
      border_radius = NULL,
      bar_bg = "#00838f",
      bar_success_bg = "#267",
      bar_warning_bg = "#00838f",
      bar_danger_bg = "#00838f",
      bar_info_bg = "#00838f"
    ),
    bs_vars_modal(
      
      backdrop_opacity = 0,
      header_border_color = "#AFDBB3",
      footer_border_color = "#FFF"
    ),
    bs_vars_button(
      font_weight = 850,
      border_radius_base = "10px",
      default_color = "#000",
      default_border = "#123F34",
      primary_color = "#FFF",
      primary_bg = "#FCAA1C",
      primary_border = "#123F34",
      success_color = "#123F34",
      success_bg = "#EFF8CD",
      success_border = "#0A5D32",
      info_color = "#123F34",
      info_bg = "#cedbd3", # bspover style=info
      info_border = "#2b6777",
      ### disgenet kısmında dropdown
      warning_color = "#17414E",
      warning_bg = "#cedbd3",
      warning_border = "#AFDBB3",
      
    )
    
    ,
    bs_vars_pills(
      border_radius = "30px",
      active_link_hover_bg = "#00838f",
      active_link_hover_color = "#FFF"
    ),
    bs_vars_wells(
      bg = "#cedbd3" # sidebar color and wellPanel
    ),
    # cedbd3
    output_file = NULL
  )),
  # This part arranges the navbar head (sanırım)
  tags$head(tags$style(HTML(".navbar-brand {height: 75px;width: 280px; font-size: 28pt;}
  .navbar-nav {height: 75px; font-size: 20pt;  }
  .navbar-nav li { height: 75px;width: 160px;display: table;}
  .navbar-nav li a {height: 75px; display: table-cell;vertical-align: middle;text-align: center;}
  .nav-tabs {font-size: 12pt;}

 "))
  ),
  
  
  tags$head(
    HTML(
      "  <script>
          var socket_timeout_interval
          var n = 0
          $(document).on('shiny:connected', function(event) {
          socket_timeout_interval = setInterval(function(){
          Shiny.onInputChange('count', n++)
          }, 15000)
          });
          $(document).on('shiny:disconnected', function(event) {
          clearInterval(socket_timeout_interval)
          });
          </script>
          "
    )
  ), 
  tags$style(type = "text/css",
             "label { font-size: 16px; }"
  ),
  useShinydashboard(),  shinydisconnect::disconnectMessage(
    text = "Your session timed out, reload the application!",
    refresh = "Reload now",
    background = "#00838f",
    colour = "white",
    overlayColour = "grey",
    overlayOpacity = 0.75,
    top = 250,
    refreshColour = "brown"
  ),   tags$head(tags$style("
      .myRow1{height:80px;color:#FFFFFF;background-color:#00838f;width:1000px;}
      .myRow2{height:600px;background-color: white;}
      .myRow3{height:600px;background-color: #EFF8CD;}
                                            " )
  ),  tags$head(
    tags$style(type = 'text/css',".myclass1 {background-color: #EFF8CD;}"),
    tags$style(type = 'text/css',".myclass2 {height:100px;color:#17414E;background-color: #EFF8CD;}"),
    tags$style(type = 'text/css',".myclass3 {width:360px;color:#17414E;background-color: #EFF8CD;}")
  ), 
  # download buttons
  htmltools::tags$head(htmltools::tags$style(".mybutton{color:black;background-color:#FCAA1C;} .skin-red .sidebar .mybutton{color: #114109;}.mybutton:hover {background-color:#BD7F13;transition: 0.01s;  }") ),
  
  tags$style(HTML(".tooltip > .tooltip-inner {width: 400px;color: #17414E;background-color: #cedbd3;}")),  
  tags$head(tags$style(HTML("
      .table.dataTable tbody td.active, .table.dataTable tbody tr.active td {
            background-color: #AFDBB3!important;}
      "))),
  tags$head(tags$style(HTML("table.dataTable.hover tbody tr:hover, table.dataTable.display tbody tr:hover {background-color: #cedbd3 !important;}"))),
  tags$style(HTML(".dataTables_wrapper .dataTables_length,.dataTables_wrapper .dataTables_buttons, .dataTables_wrapper 
                  .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button {color: #123F34 !important;}
                  
                  #imgid{
  position: relative;
  left: 0px;
  top: -26px;
   bottom: 10;
}

")),
  
  navbarPage(
    
    # Application title.
    
    title = div(div(id="imgid",
      img(
        src = "LogoW.svg",
        height = 100,
        width = 180,
       # style = "background-color: white;"
      ))
    ),
      #tags$div(class = "text-left", style = "margin-bottom:0px;margin-top:0px",img(src="logotitle.png",align="top", height = '50', width ='170')),
   
      #div(span("eMRald", style = "font-size:36px;color:white; position: relative; top: 60%; transform: translateY(-60%);"),windowTitle = "Welcome to eMRald!"),
    id = "tabs",
 
       
     
    
    
    shiny::tabPanel("Home",
                                     tabBox(title = "",width = NULL,
                                            shiny::tabPanel(title = "Welcome",icon = icon("info"),fluidRow(column(includeMarkdown("mdfiles/welcome.md"),width =11,offset = 1))),  
                                            shiny::tabPanel(title = "Data Explanation",icon = icon("dna"),fluidRow(column(includeMarkdown("mdfiles/dataexplanationtab.md"),width =12,offset = 0))),
                                            shiny::tabPanel(
                                              title = "User Guide",icon = icon("book"),
                                              navlistPanel( 
                                               # id = "tabset",
                                               widths = c(2, 10),
                                                
                                                tabPanel("Data Selection", h3("Data Selection"),
                                                         tags$hr(),
                                                         column(
                                                          includeMarkdown("mdfiles/dataselection.md"),
                                                           width =12,
                                                           offset = 0)),
                                                tabPanel("Gene Selection",h3("Gene Selection"),
                                                         tags$hr(),
                                                         column(
                                                          includeMarkdown("mdfiles/geneselection.md"),
                                                           width =12,
                                                           offset = 0)),
                                               tabPanel("Scatter Plot",h3("Scatter Plot"),
                                                        tags$hr(),
                                                        column( includeMarkdown("mdfiles/scatter.md")  ,
                                                                width =12,
                                                                offset = 0)),
                                                tabPanel("Heatmap", h3("Heatmap"),
                                                         tags$hr(),
                                                         column(
                                                           includeMarkdown("mdfiles/heatmap.md"),
                                                           width =12,
                                                           offset = 0)),
                                                tabPanel("Functional Profile", h3("Functional Profile"),
                                                         tags$hr(),
                                                         column(
                                                         includeMarkdown("mdfiles/functionalprofile.md"),
                                                           width =12,
                                                           offset = 0)),
                                               tabPanel("Transcription Factors", h3("Transcription Factors"),
                                                        tags$hr(),
                                                        column(
                                                           includeMarkdown("mdfiles/TFpart.md"),
                                                          width =12,
                                                          offset = 0)),
                                                tabPanel("MSigDB", h3("MSigDB"),
                                                         tags$hr(),
                                                         column(
                                                           includeMarkdown("mdfiles/msigdb.md"),
                                                           width =12,
                                                           offset = 0)),
                                                tabPanel("Compare Clusters", h3("Compare Clusters"),
                                                         tags$hr(),
                                                         column(
                                                           includeMarkdown("mdfiles/comparecluster.md"),
                                                           width =12,
                                                           offset = 0))
                                                
                                              ),
                                            )
                                     )
                    ,textOutput("keepAlive") ),#tabPanel ends
    
    shiny::tabPanel("Analysis",
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(width = 250,
                                          fluidRow(
                                            column(width = 3, 
                                                   wellPanel(
                                                     fluidRow( selectInput(inputId = "experiment_title_1",label = "Data 1:",choices = l)),
                                                     fluidRow( column(width = 6, uiOutput( "select_control_output")),column(width = 6, uiOutput( "select_treatment_output"))))
                                            ),
                                            column(width = 3,
                                                   wellPanel(
                                                     fluidRow(selectInput(inputId = "experiment_title_2",label = "Data 2:",choices = c("Select",l))),
                                                     fluidRow(column(width = 6, uiOutput(outputId = "select_control_output_2")),column(width = 6, uiOutput(outputId = "select_treatment_output_2"))))
                                            ),
                                            column(width = 3, 
                                                   wellPanel(
                                                     fluidRow( selectInput(inputId = "experiment_title_3",label = "Data 3:",choices = c("Select",l))),
                                                     fluidRow(column(width = 6, uiOutput(outputId = "select_control_output_3")),column(width = 6, uiOutput(outputId = "select_treatment_output_3"))))
                                            ),
                                            
                                            column(3,
                                                   wellPanel(
                                                     fluidRow( selectInput(inputId = "experiment_title_4",label = "Data 4:",choices = c("Select",l))),
                                                     fluidRow(column(width = 6, uiOutput(outputId = "select_control_output_4")),column(width = 6, uiOutput(outputId = "select_treatment_output_4")))
                                                   )
                                            )  ),# fluidrow,
                                          
                                          
                                          
                                          
                                          fluidRow( column(3,selectInput("metapvalue","Method for the meta–analysis of p–values", 
                                                                         c("Fisher’s method"="fisher","Stouffer’s method"="stouffer","Inverse chi-square method"="invchisq","Binomial test"="binomtest","Bonferroni method"="bonferroni","Tippett’s method"="tippett"),
                                                                         selected="fisher")),
                                                    
                                                    column(width = 2,offset = 1, actionBttn(inputId = "analyze", label = "Analyze", icon = icon("play"), color = "warning" )),
                                                    column(3,  selectInput(inputId = "colname",label = "Select columns to display", 
                                                                                    choices = c("LogFC"="LogFC", "Average Exp"="AveExp", "t value"="t.value","p value"="p.value","p adj value"="p.adj.value", "B value"="B.value"),
                                                                                    #inline = TRUE,
                                                                           selected = c("LogFC", "p.value"),multiple = T,
                                                                           selectize = TRUE 
                                                    ) ))
                      ), # sidebarPanel ends
                      shiny::mainPanel(width = 250,
                                       
                                       tabsetPanel(type = "tabs",
                                                   tabPanel( icon =icon("info-circle"), title ="Data Summary",
                                                            conditionalPanel(condition="input.analyze != 0",
                                                                             panel(class = "myclass2",style="font-size: 16px; line-height: 16px;",
                                                                                   HTML(paste0("<p>'Volcano plot' tab shows the size of change (fold change) in x-axis with statistical significance (P value) in y-axis.<p>",
                                                                                               " <p>'Plots' include bar/box and PCA plot. For RNASeq data, bar plots shows library size, and PCA plot uses rlog. 
                                                                                                  For Microarray data, Box-and-whisker plot shows the expressions in each sample by removing outliers, and PCA plot uses prcomp .<p>",
                                                                                               "<p>'Raw Data' tab shows raw count data.<p>"
                                                                                   ))
                                                                                   
                                                                             )   ),
                                                            
                                                            
                                                            
                                                            uiOutput(outputId = "data_summary_ui_output")),
                                                   tabPanel( title = "Results", icon = icon("th", class = NULL, lib = "font-awesome"),
                                                             uiOutput(outputId = "results_ui_output"))
                                       )#tabsetPanel ends
                      )# mainPanel ends
                    )
                    
    )
    
    
    
  )
)

