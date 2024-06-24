 # This is a Shiny web application for visualizing HSML profiles  --------


# Load libraries ----------------------------------------------------------

library(shiny)
library(dplyr)
library(stringr)
library(plotly)
library(htmlwidgets)
library(shinyHugePlot)
library(bslib)
library(duckdb)
library(arrow)

options(shiny.maxRequestSize=50*1024^2) #max dim for input files 




# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)
addResourcePath(prefix = 'www', directoryPath = '/app/www')
shinyOptions(cache = cachem::cache_mem(max_size = 500e6))


# Read variants annotations data ---------------------------------------------------

Annotations_37 <- open_dataset("/app/37")
Annotations_38 <- open_dataset("/app/38")

# Read genes annotation data ----------------------------------------------

genes_annotation_38 <-readRDS("/app/genes_annotation_38.rds")
genes_annotation_37 <-readRDS("/app/genes_annotation_37.rds")



# Define UI for application that plots the output of EXCAVATOR2 ---------


ui  <- page_sidebar(
  theme = bs_theme(
  bootswatch = "cosmo"),

 # Application title
  title ="ReViewCNV",
    sidebar = sidebar(width = "35%",
      fileInput("Bed", "Load the bed file of annotated targeted regions"),
      fileInput("HSLM_1", "Select the HSLMResults_* file"),
      fileInput("FastCall_Results_1", "Load the FastCall_Result from Results folder"),
        selectInput("Genome", "Select the Genome version", selected = NULL,
        choices = c("","GRCh37", "GRCh38")),
      checkboxInput("GenomeBrowser", "Show genes annotations", value = FALSE, width = NULL),
      checkboxInput("ShareAxes", "Share x axis", value = TRUE, width = NULL),
      checkboxInput("Set_y_axis", "Set the same Log2R range", value = FALSE, width = NULL),
      uiOutput("Button_1"), 
      uiOutput("Second_Individual"),
      uiOutput("Button_2"), 
      uiOutput("Third_Individual"),
        selectInput("chr", "Select the Chr", selected = "chr1",
          choices = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
            "chr7", "chr8", "chr9", "chr10", "chr11",
           "chr12", "chr13", "chr14", "chr15", "chr16",
           "chr17", "chr18", "chr19", "chr20", "chr21",
            "chr22", "chrX")),
        selectInput("Prob", "Select Calls", selected = "All",
          choices = c("All" = "-1", "Prob Call \u2265  0.5" = "0.5", "Prob Call \u2265  0.6" = "0.6", 
            "Prob Call \u2265  0.7" = "0.7", "Prob Call \u2265  0.8" = "0.8",
            "Prob Call \u2265  0.9" = "0.9")),
        h4("Choose CNVs datasets"),
        checkboxInput("AnnotSV", "AnnotSV", value = TRUE),
        checkboxInput("DGV_Merge", "DGV_Merge", value = FALSE),
        checkboxInput("DGV_Gold", "DGV_Gold", value = FALSE),
        checkboxInput("GnomAD_Genome", "GnomAD_Genome", value = FALSE),
        checkboxInput("GnomAD_Exome", "GnomAD_Exome (only GRCh38)", value = FALSE),
        selectInput("Type", "Select CNVs type", selected = "All",
                  choices = c("All" = "%","Gain" = "Gain", 
                              "Loss" = "Loss")),
        selectInput("Freq", "Select CNVs frequency", selected = "All",
            choices = c("All" = 0, "Freq \u003E  0.01" = 0.01, 
              "Freq \u003E  0.02" = 0.02,"Freq \u003E  0.05" = 0.05, 
              "Freq \u003E  0.1" = 0.1, "Freq \u003E  0.2" = 0.2,
              "Freq \u003E  0.3" = 0.3)),
        uiOutput("limits"),
        uiOutput("Annotations_limits"),
        HTML('<img src = "www/Logo.png" width = "60%" hight = "auto" >')
        ),
        plotlyOutput("zoomPlot"),
        uiOutput("Download")
        )


# Server ------------------------------------------------------------------


server <- function(input, output, session) {
  bs_themer()
  session_store <- reactiveValues()
  options(warn = -1)

  
  
  

# Importing Bed -----------------------------------------------------------

  
  target_data<- reactive({
    if (is.null(input$Bed)) {
      return(NULL)
     }
    else{
     target_data <-read.table(input$Bed$datapath, quote="\"",fill=T,
      col.names= c("Chromosome","Start","End", "Exon"))
    }
  })
  

  
# File and target data first individual -----------------------------------

    file_data_1 <- reactive({
      if (is.null(input$HSLM_1)) {
      return(NULL)
      }
    else{
     read.table(input$HSLM_1$datapath, fill=T, quote="\"", header = T) |> 
        mutate(Log2R = round(Log2R,2))
    }
    })
      
    fast_call_1 <- reactive({
      if (is.null(input$FastCall_Results_1)) {
        return(NULL)
        }
      else{
        read.table(input$FastCall_Results_1$datapath, header = T,
          fill=T, quote="\"") |> 
        mutate(Mutation = case_when(
          CN == 0 ~ "2-DEL",
          CN ==1 ~ "DEL",
          CN == 3 ~ "AMP",
          TRUE ~ "2-AMP")) |> 
        filter(ProbCall >= as.numeric(input$Prob)) |> 
        mutate(mid = (Start + End)/2 )
        }
    })
  
# File and FastCall data second individual ----------------------------------

  file_data_2 <- reactive({
    if (is.null(input$HSLM_2)) {
      return(NULL)
      }
    else{
      read.table(input$HSLM_2$datapath, fill=T, quote="\"", header = T)  |> 
      mutate(Log2R = round(Log2R,2))
      }
    })

  fast_call_2 <- reactive({
    if (is.null(input$FastCall_Results_2)) {
      return(NULL)
    }
    else{
      read.table(input$FastCall_Results_2$datapath,
      fill=T, quote="\"", header = T) |> 
      mutate(Mutation = case_when(
        CN == 0 ~ "2-DEL",
        CN ==1 ~ "DEL",
        CN == 3 ~ "AMP",
        TRUE ~ "2-AMP")) |> 
      filter(ProbCall >= as.numeric(input$Prob)) |> 
      mutate(group=seq_along(Start)) |> 
      mutate(mid = (Start + End)/2 )
      }
    })

# File and FastCall data third individual -----------------------------------

  file_data_3 <- reactive({
    if (is.null(input$HSLM_3)) {
      return(NULL)
      }
    else{
      read.table(input$HSLM_3$datapath, fill=T, quote="\"", header = T) |> 
      mutate(Log2R = round(Log2R,2))
      }
    })
    
  fast_call_3 <- reactive({
    if (is.null(input$FastCall_Results_3)) {
      return(NULL)
      }
    else{
      read.table(input$FastCall_Results_3$datapath,
        fill=T, quote="\"", header = T) |> 
      mutate(Mutation = case_when(
        CN == 0 ~ "2-DEL",
        CN ==1 ~ "DEL",
        CN == 3 ~ "AMP",
        TRUE ~ "2-AMP")) |> 
      filter(ProbCall >= as.numeric(input$Prob)) |> 
      mutate(group=seq_along(Start)) |> 
      mutate(mid = (Start + End)/2 )
      }
    })
# Merging targets and  file data----------------------------------------------------


  
  merge_data_1 <- reactive({
    left_join(file_data_1(),target_data(),by=c("Chromosome","Start", "End"))
   })

  merge_data_2 <- reactive({
    left_join(file_data_2(),target_data(),by=c("Chromosome","Start", "End"))
    })
    
  merge_data_3 <- reactive({
    left_join(file_data_3(),target_data(),by=c("Chromosome","Start", "End"))
    })
  
# Subsetting merge data--------------------------------------------------------------

#First individual  

  subset_data_1 <- reactive({
   merge_data_1()[merge_data_1()$Chromosome == input$chr,]  |>
    mutate(Text = paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon))
    })
  
  subset_data_1_range <- reactive({subset_data_1() |> 
    filter(Start >= input$slider[1] & End <= input$slider[2]) 
    })
  
  
#Second individual 
  
  subset_data_2 <- reactive({
    merge_data_2()[merge_data_2()$Chromosome == input$chr,] |>
    mutate(Text = paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon))
    })
  
  subset_data_2_range <- reactive({subset_data_2() |> 
    filter(Start >= input$slider[1] & End <= input$slider[2])
  })
  

#Third individual  
   
  subset_data_3 <- reactive({
    merge_data_3()[merge_data_3()$Chromosome == input$chr,] |>
    mutate(Text = paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon))
    })
  
  subset_data_3_range <- reactive({subset_data_3() |> 
    filter(Start >= input$slider[1] & End <= input$slider[2])
  })  
  

# Subsetting FastCall data ------------------------------------------------

#First individual
  
  rects_1 <- reactive({fast_call_1() |> filter(Chromosome == input$chr)})
  
  rects_1_db <- reactive({ rects_1() |> to_duckdb()})
  
  rects_1_range <- reactive ({
    if (is.null(input$slider[1])) {
      return(NULL)
     }
    else{rects_1() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })

#Second individual   
  
  rects_2 <- reactive({fast_call_2() |> filter(Chromosome == input$chr)})
  
  rects_2_range <- reactive ({
    if (is.null(input$slider[1])) {
      return(NULL)
     }
    else { rects_2() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })
   
#Third individual
  
  rects_3 <- reactive({fast_call_3() |> filter(Chromosome == input$chr)})
  
  rects_3_range <- reactive ({
    if (is.null(input$slider[1])) {
     return(NULL)
     }
   else{rects_3() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })
  
  
# Subsetting variants annotations data ---------------------------------------------
 
  Annot_SV_D <- reactive({
    if(input$AnnotSV){"AnnotSV BenignSV (v. 3.4)"}
    else{return(NULL)}
  })

  DGV_Merge_D <- reactive({
    if(input$DGV_Merge){
      "dgvMerged (last updated 2020-02-25)"}
    else{return(NULL)}
  })

  DGV_Gold_D <- reactive({
      if(input$DGV_Gold){
       "dgvGold (last updated 2016-05-15)"}
    else{return(NULL)}
  })

  GnomAD_D <- reactive({
    if(input$GnomAD_Genome){
      c("gnomad_v2.1_sv.controls_only.site", "gnomad.v4.1.sv.non_neuro_controls.sites")}
    else{return(NULL)}
  })
  
  GnomAD_Exome_D <- reactive({
    if(input$GnomAD_Exome){
      GnomAD_Exome_D <-  "gnomad.v4.1.cnv.non_neuro_controls"}
    else{return(NULL)}
  })
  

  Annotations_list <-reactive({
    Annotations_list<-c(Annot_SV_D(), DGV_Gold_D(), DGV_Merge_D(), GnomAD_D(),  GnomAD_Exome_D())
    Annotations_list
  })
  
  
  
  Annotations_subset <-reactive({
    if(input$Genome == "GRCh37"){
      
      Annotations_37 |>  
        to_duckdb() |>
        filter(Chromosome %in% !!fast_call_1()$Chromosome) |>  
        filter (Chromosome == input$chr) |> 
        filter(calls %like% input$Type)|> 
        filter(Frequency > input$Freq) |> 
        filter (Database %in%  !!Annotations_list() ) |> 
        inner_join(rects_1_db() , join_by(overlaps(Start, End, Start, End)))|>
        rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |> 
        select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |> 
        arrange(End_FastCall, End) |> 
        collect()}


    else if(input$Genome == "GRCh38"){
      Annotations_38 |>
        to_duckdb() |>
        filter(Chromosome %in% !!fast_call_1()$Chromosome) |>  
        filter (Chromosome == input$chr) |> 
        filter(calls %like% input$Type) |>  
        filter(Frequency > input$Freq) |> 
        filter (Database %in% !!Annotations_list() ) |> 
        inner_join(rects_1_db(), join_by(overlaps(Start, End, Start, End)))|>
        rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |> 
        select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |> 
        arrange(End_FastCall, End) |> 
        collect()}
    
    else {return(NULL)}
    
  }) %>% bindCache (input$Genome, input$chr, input$Type, input$Freq, Annotations_list(), fast_call_1())
  
  
  CNV1 <- reactive ({
    if(!is.null(Annotations_subset())){
    if(input$AnnotSV){
      Annotations_subset() |>  filter( AnnotSV_Present == "No")  }
    else {Annotations_subset()}}
    else{return(NULL)}
  })

  
  CNV_FC <- reactive({
    if(!is.null(Annotations_subset())){
      CNV1() |> 
      select(Start_FastCall, End_FastCall) |> 
      distinct(Start_FastCall, End_FastCall)}
    else{return(NULL)}
      
      })
  

  CNV2 <-reactive({
    if(!is.null(CNV1())){
      v <- c(1)
      k = 1
      if(dim(CNV1())[1]>0){
        if(dim(CNV1())[1]>1){
        for (i in (2:dim(CNV1())[1])){
          if(CNV1()[i,]$Start < CNV_FC()[k,]$End_FastCall)  {
            v <-append(v,k)}
          else {k = k +1
            v<-append(v,k)}}} 

       CNV2 <- CNV1()|> 
        mutate(n_overlap = v) 
       
       for (i in (2:dim(CNV2)[1])){
         if (CNV2[i,]$Unique_ID  %in% CNV2[-i,]$Unique_ID){
           CNV2[i,]$n_overlap <-0}}
      
           
      CNV2 <- CNV2 |> group_by(Unique_ID) |> 
             slice(1) |> 
             ungroup() |> 
             arrange(n_overlap, desc(round(log10(Length),0)), desc(Frequency))
      CNV2
       
      }}
    else{return(NULL)}
    }) |> bindCache(CNV1(), CNV_FC())
  
  CNV3 <- reactive({
    if(!is.null(CNV2())){
      v2 <- c(1)
      k2 = 1
     if(dim(CNV2())[1]>0){
        if(dim(CNV2())[1]>1){
          for (i in (2:dim(CNV2())[1])){
            if(CNV2()[i,]$n_overlap == CNV2()[i-1,]$n_overlap) {
              k2 = k2+1
              v2<-append(v2,k2)}
            else {k2 = 1
              v2<-append(v2,k2)}
            }}
  
        CNV3 <- CNV2()|> 
        mutate(level = v2) 
        

      if(dim(CNV2())[1] >1){
        i = 2
        while (i <= dim(CNV3)[1]){
          if(CNV3[i,]$Start <  max(CNV3[1:i-1,] |> filter(level == CNV3[i,]$level) |> select(End))){
            CNV3[i,]$level <-  CNV3[i,]$level +1 
            i = i}
          else{i = i+1}
          }}
    
        
      CNV3
      }}
   else{return(NULL)}
   }) |> bindCache(CNV2())
  

# Subsetting genes annotations data ---------------------------------------

genes_annotations <-reactive ({
  if(input$GenomeBrowser){
  
    if (input$Genome == "GRCh37"){
      genes_annotation_37 |>  
      filter (Chr == input$chr)}
    else if (input$Genome == "GRCh38"){
      genes_annotation_38 |>  
      filter (Chr == input$chr)
    }
  }
  else(return(NULL))
  }) |> bindCache(input$Genome, input$chr)
  

# First button -------------------------------------------------------

   x <- reactiveValues(val = 0)
   
   observeEvent(input$Button_1, {
     x$val = x$val +1 
   })
   
   output$Button_1  =renderUI({
     if(x$val == 0){
      actionButton("Button_1", "Add new sample")
       }
     else {return(NULL)}
   })
   
   
    output$Second_Individual=renderUI({
      if(x$val == 1){list(
        fileInput("HSLM_2", "Select the HSLMResults_* file"),
         fileInput("FastCall_Results_2", "Load the FastCall_Result from Results folder"))
        }
      else{return(NULL)}
    
     })



# Second button -----------------------------------------------------------

    z <- reactiveValues(val = 0)
    
    observeEvent(input$Button_1, {
      z$val = z$val +1 
    })
    
    observeEvent(input$Button_2, {
      z$val = z$val +1 
    })
    
    output$Button_2  =renderUI({
      if(z$val == 1){
        actionButton("Button_2", "Add new sample")
      }
      else {return(NULL)}
    })
    
    
    output$Third_Individual=renderUI({
      if(z$val == 2){list(
        fileInput("HSLM_3", "Select the HSLMResults_* file"),
        fileInput("FastCall_Results_3", "Load the FastCall_Result from Results folder"))
      }
      else{return(NULL)}
      
    })
     
    




# Setting the range slider for coordinates------------------------------------------------
   
   output$limits=renderUI({
     if (!(is.null(input$HSLM_1)|is.null(input$Bed))){
       
       sliderInput('slider','Choose the chromosome coordinates for the download',
                   min=min(subset_data_1()$Start),
                   max=max(subset_data_1()$End),
                   value=c(min(subset_data_1()$Start),
                           max(subset_data_1()$End)))}
     
     else {return(NULL)}
   })
   


# Setting the range slider for annotations --------------------------------

    
    output$Annotations_limits=renderUI({
      if (!(is.null(input$HSLM_1)|is.null(input$Bed)|is.null(input$FastCall_Results_1))){
        if(!(is.null(CNV3()))){
          if (dim(CNV3())[1] >0){
            sliderInput('slider_Annotations',
            'Choose the maximum number of overlapping CNVs to visualize',
            min=0,
            max=max(CNV3()$level),
            value=c(0,
            min(20,max(CNV3()$level))))}
          }}
      
    
      else {return(NULL)}
    })
    
    



# Download ----------------------------------------------------------------


output$Download = renderUI({
  if (!(is.null(input$HSLM_1)|is.null(input$Bed)|is.null(input$FastCall_Results_1))){
    if (input$Genome  != ""){
      downloadButton("downloadplot", "Download HTML")}
  }
  else{return(NULL)}
  })
 
     

# Plot variants annotations --------------------------------------------------------

    output$zoomPlot <- renderPlotly({
      
     
     
if(!(is.null(input$HSLM_1)|is.null(input$Bed) | is.null(input$FastCall_Results_1) | is.null(input$slider_Annotations[1])|
     is.null(input$slider_Annotations[2]) )){
  if(!(is.null(CNV3()))){
    if(dim(CNV3())[1]>0){
      
    
    
     CNV_Gain <- CNV3() |>
       filter(calls == "Gain") 
  

     CNV_Loss <- CNV3() |>
       filter(calls == "Loss")
    
     

  
  rect_1_CNV_Gain <- list(
    type ="polygon",
    fillcolor = "blue",
    line = list( color = "blue" ))
  
  
  
  rect_CNV_Gain <- list()
  for (i in c(1:dim(CNV_Gain)[1])) {
    rect_1_CNV_Gain[["x0"]] <- CNV_Gain[i,]$Start
    rect_1_CNV_Gain[["x1"]] <- CNV_Gain[i,]$End
    rect_1_CNV_Gain[["y0"]] <- CNV_Gain[i,]$level -0.6
    rect_1_CNV_Gain[["y1"]] <- CNV_Gain[i,]$level -0.3
    rect_CNV_Gain <- c(rect_CNV_Gain, list(rect_1_CNV_Gain))
  }
  
  

  
  
  rect_1_CNV_Loss <- list(
    type ="polygon",
    fillcolor = "orange",
    line = list( color = "orange" ))
  
  
  rect_CNV_Loss <- list()
  for (i in c(1:dim(CNV_Loss)[1])) {
    rect_1_CNV_Loss[["x0"]] <- CNV_Loss[i,]$Start
    rect_1_CNV_Loss[["x1"]] <- CNV_Loss[i,]$End
    rect_1_CNV_Loss[["y0"]] <- CNV_Loss[i,]$level - 0.6
    rect_1_CNV_Loss[["y1"]] <- CNV_Loss[i,]$level - 0.3
    rect_1_CNV_Loss[["name"]] <- CNV_Loss[i,]$ID 
    rect_CNV_Loss <- c(rect_CNV_Loss, list(rect_1_CNV_Loss))
  }
  
  
  
 
  rect <- c()
  
  if (dim(CNV_Loss)[1] >0){
    rect <- append(rect, rect_CNV_Loss)
  }
  
  if (dim(CNV_Gain)[1] >0){
    rect <- append(rect, rect_CNV_Gain)
  }
 fig <- plot_ly() |> 
   layout(shapes = rect, xaxis = list(title = "Chromosome coordinates",
          range = list(min(subset_data_1_range()$Start),
          max(subset_data_1_range()$End))),
          yaxis = list(title = "CNVs",range =list((input$slider_Annotations[1]),
          (input$slider_Annotations[2])), tickformat=',d'), type = "scatter")|>
    add_trace(data = CNV_Gain, type = 'scatter', mode = 'markers', x =  ~middle, y = ~level -0.45, color = I("blue"),
                        text = ~paste("Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>", 
                        "End:", End, "<br>", "Length: ", Length, "<br>",
                        "Reported frequency: ", Frequency),
                        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                        hoverlabel = list(bgcolor = "blue", align = "left"))|>
    add_trace(data = CNV_Loss, type = 'scatter', mode = 'markers', x = ~middle, y = ~level -0.45,color = I("orange"),
                        text = ~paste("Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>", 
                        "End:", End, "<br>", "Length: ", Length,"<br>",
                        "Reported frequency: ", Frequency),
                        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                        hoverlabel = list(bgcolor ="orange", align = "left")) |> 
     partial_bundle()  |>  toWebGL()
  }
  }
  else{
    fig <- plot_ly(type = 'scatter', mode = "markers") |> 
      layout( yaxis = list(range =list(0, 20))) |> 
      partial_bundle()  |>  toWebGL()
  }
}

else{return(NULL)}

      

# Plot genes annotations --------------------------------------------------
if(input$GenomeBrowser){
  

  
  rect_1_genes_annotation <- list(
    type ="line",
    line = list( color = "#008000" ))
      
      
      
  rect_genes_annotation <- list()
    for (i in c(1:dim(genes_annotations())[1])) {
      rect_1_genes_annotation[["x0"]] <- genes_annotations()[i,]$Start
      rect_1_genes_annotation[["x1"]] <- genes_annotations()[i,]$End
      rect_1_genes_annotation[c("y0", "y1")] <- genes_annotations()[i,]$level 
      rect_genes_annotation <- c(rect_genes_annotation, list(rect_1_genes_annotation))
      }
      
      
      
  fig2 <- plot_ly(data = genes_annotations(), type = 'scatter', mode = 'markers', x =  ~middle_gene, y = ~level, color = I("#008000"),
          text = ~paste("RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol, 
          "<br>","Start:", Start, "<br>", "End:", End),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#008000" ,align = "left")) |> 
        layout(shapes = rect_genes_annotation, xaxis = list(title = "Chromosome coordinates"),
          xaxis = list(range =list((min(subset_data_1_range()$Start)),
          max(subset_data_1_range()$End))), 
          yaxis = list( showline=  F, tickvals = c(1, 2), title = "Direction",
          ticktext = c('+', '\u2212'),  tickfont = list(size = 20, family = "Arial black"))) |> 
         partial_bundle()  |>  toWebGL()
  }
 
# Plot for first individual -----------------------------------------------
   

     if(!(is.null(input$HSLM_1)|is.null(input$Bed) | is.null(input$FastCall_Results_1))){
         
       pal <- c("blue", "lightblue")
       
       pal <- setNames(pal, c("IN", "OUT"))
       
       
       
       rects_2DEL <- 
         rects_1_range() |> 
         filter(Mutation == "2-DEL")
       
       
       rects_DEL <- 
         rects_1_range() |> 
         filter(Mutation == "DEL")
       
       
       rects_AMP <- 
         rects_1_range() |> 
         filter(Mutation == "AMP")
       
       
       rects_2AMP <- 
         rects_1_range() |> 
         filter(Mutation == "2-AMP")
       
       
       
       
       
       rect_1D <- list(
         type ="rect",
         fillcolor =   "#e69f00",
         line = list(color = "#e69f00"),
         opacity = 0.6
       )
       
       rect_D <- list()
       
       
       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- floor(min(subset_data_1()$Log2R) - 0.5)
         rect_1D[["y1"]] <- ceiling(max(subset_data_1()$Log2R) + 0.5)
         rect_D <- c(rect_D, list(rect_1D))
       }
       
       
       rect_1A <- list(
         type ="rect",
         fillcolor = "#56B4E9",
         line = list(color =  "#56B4E9"),
         opacity = 0.6
       )
       
       rect_A <-  list()
       for (i in c(1:dim(rects_AMP)[1])) {
         rect_1A[["x0"]] <- rects_AMP[i,]$Start
         rect_1A[["x1"]] <- rects_AMP[i,]$End
         rect_1A[["y0"]] <- floor(min(subset_data_1()$Log2R) - 0.5)
         rect_1A[["y1"]] <- ceiling(max(subset_data_1()$Log2R)+0.5)
         rect_A <- c(rect_A, list(rect_1A))
       }
       
       
       
       rect_1_2A <- list(
         type ="rect",
         fillcolor = "#0072B2",
         line = list(color = "#0072B2"),
         opacity = 0.6
       )
       
       rect_2A <- list()
       for (i in c(1:dim(rects_2AMP)[1])) {
         rect_1_2A[["x0"]] <- rects_2AMP[i,]$Start
         rect_1_2A[["x1"]] <- rects_2AMP[i,]$End
         rect_1_2A[["y0"]] <- floor(min(subset_data_1()$Log2R) - 0.5)
         rect_1_2A[["y1"]] <- ceiling(max(subset_data_1()$Log2R)+0.5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }
       
       
       
       
       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#d55e00",
         line = list(color = "#d55E00"),
         opacity = 0.6
       )
       
       rect_2D <- list()
       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- floor(min(subset_data_1()$Log2R) - 0.5)
         rect_1_2D[["y1"]] <- ceiling(max(subset_data_1()$Log2R)+0.5)
         rect_2D <- c(rect_2D, list(rect_1_2D))}
       
       
       
       rect<-c()
       
       if (dim(rects_AMP)[1] >0){
         rect <- append(rect, rect_A )
       }
       
       if (dim(rects_DEL)[1] >0){
         rect <- append(rect, rect_D) 
       }
       
       if (dim(rects_2AMP)[1] >0){
         rect <- append(rect, rect_2A) 
       }
       
       if (dim(rects_2DEL)[1] >0){
         rect <- append(rect, rect_2D) 
       }   
       
       
       pl_1<- plot_ly(data = subset_data_1_range(), x = ~Position) |> 
       add_trace(y = ~Log2R, color = ~Class, colors =pal, type = 'scatter', mode = "markers",
        text = ~paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon),
        hoverinfo = "text", 
        marker =list(size = 4)) |>
      add_trace(y = ~SegMean, color = "Seg", type = 'scatter', mode = "lines", line  = list(color = 'red')) |> 
      add_trace(type = 'bar', 
         x = min(subset_data_1_range()$Position), 
         y = 0, 
         name = "DEL", 
         color= I("#E69F00"),
         opacity = 0.6) |> 
      add_trace(type = 'bar', 
         x = min(subset_data_1_range()$Position), 
         y = 0, 
         name = "AMP", 
         color = I("#56B4E9"),
        opacity = 0.6)|>
      add_trace(type = 'bar', 
         x = min(subset_data_1_range()$Position), 
         y = 0, 
         name = "2_DEL", 
         color = I("#D55E00"),
         opacity = 0.6)|> 
      add_trace(type = 'bar', 
        x = min(subset_data_1_range()$Position), 
        y = 0, 
        name = "2_AMP", 
        color= I("#0072B2"),
        opacity = 0.6) |> 
      add_trace(type = 'bar', 
        x = min(subset_data_1_range()$Position), 
        y = 0, 
       name = "Gain CNVs from selected databases", 
       color= I("blue"))|>
    add_trace(type = 'bar', 
      x = min(subset_data_1_range()$Position), 
      y = 0, 
     name = "Loss CNVs from selected databases", 
     color= I("orange")
       ) |> 
    layout(shapes = rect, legend=list(title=list(text='Window')),
                       yaxis = list(title = 'Log2Ratio'), xaxis = list(title = "Chromosome coordinates"))

       if (dim(rects_2DEL)[1] >0){
       pl_1<- pl_1|> add_trace(data = rects_2DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_1()$Log2R)+0.5),   color = I("#d55e00"),
                             text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                            "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                            "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                            "ProbCall: ", round(ProbCall,2)),
                             hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                             hoverlabel = list(bgcolor = "#d55e00"))|>
                 add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y = ceiling(max(subset_data_1()$Log2R)+0.5),   color = I("#d55e00"),
                            text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                            "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                            "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                            "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#d55e00") ) |> 
                add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_1()$Log2R)-0.5),   color = I("#d55e00"),
                            text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                            "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                            "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                            "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#d55e00")) |> 
                add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_1()$Log2R)-0.5),   color = I("#d55e00"),
                              text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                              "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                              "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                              "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#d55e00")) |> 
               add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#d55e00"),
                               text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                               "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#d55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_1()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00") ) |> 
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_1()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                     add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~Start, y = floor(min(subset_data_1()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |>  
                  add_trace(data = rects_DEL, type = 'scatter', mode = 'markers',  x =  ~Start, y =  ceiling(max(subset_data_1()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                               "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                  add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                               "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) }
       if (dim(rects_AMP)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_1()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                       add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  floor(min(subset_data_1()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                      add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_1()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                      add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  ceiling(max(subset_data_1()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                    add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
        if (dim(rects_2AMP)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_1()$Log2R)+0.5), color = I("#0072B2"),
                                text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                         add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  floor(min(subset_data_1()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                          add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_1()$Log2R)+0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                          add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_1()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                           add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                   "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}
       
       if(input$Set_y_axis) {
         if (is.null(input$HSLM_2) | is.null(input$FastCall_Results_2)){
           pl_1 <- pl_1|> layout(
             yaxis = list(range =c(floor(min(subset_data_1()$Log2R) -0.5),
                                      ceiling(max(subset_data_1()$Log2R)+0.5))))}
         else if (is.null(input$HSLM_3) | is.null(input$FastCall_Results_3)){
           pl_1<- pl_1|> layout(
             yaxis = list(range = c(
               min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5)), 
               max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5))
             )))}
         else {
           pl_1 <- pl_1|> layout(
             yaxis = list(range = c(
               min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5, min(subset_data_3()$Log2R) - 0.5 )), 
               max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5))))
           )}
       }
       else{
         pl_1<- pl_1|> layout(
           yaxis = list(range =c(floor(min(subset_data_1()$Log2R) -0.5),
                                    ceiling(max(subset_data_1()$Log2R)+0.5))))}
       
       if(input$GenomeBrowser){
         if (input$Genome == "GRCh38"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = min(subset_data_1_range()$Position),
                                   y = 0,
                                   name = "Genes annotation from MANE_RefSEq v 1.3",
                                   color= I("#008000"))}
         
         if (input$Genome == "GRCh37"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = min(subset_data_1_range()$Position),
                                   y = 0,
                                   name = "Genes annotation from ncbi RefSeq Select \n(last updated 2022-03-16)",
                                   color= I("#008000"))}
       }
       
      pl_1 <- pl_1 |>  partial_bundle()  |>  toWebGL()
        
        if(input$GenomeBrowser){
      
        pl<-subplot(plotly_build_light(fig2),plotly_build_light(pl_1),plotly_build_light(fig), nrows =3,heights = c(1/6,2/6,3/6),shareX = TRUE, titleY  = TRUE)}
        
        else{ pl<-subplot(plotly_build_light(pl_1),plotly_build_light(fig), nrows =2,heights = c(1/2,1/2),shareX = TRUE, titleY  = TRUE)}
      }
  

       
     

# Plot for the second individual ------------------------------------------

    

     
     
     if(!(is.null(input$HSLM_2) | is.null(input$FastCall_Results_2))){
         
       pal <- c("blue", "lightblue")
       
       pal <- setNames(pal, c("IN", "OUT"))
       
       
       
       rects_2DEL <- 
         rects_2_range() |> 
         filter(Mutation == "2-DEL")
       
       
       rects_DEL <- 
         rects_2_range() |> 
         filter(Mutation == "DEL")
       
       
       rects_AMP <- 
         rects_2_range() |> 
         filter(Mutation == "AMP")
       
       
       rects_2AMP <- 
         rects_2_range() |> 
         filter(Mutation == "2-AMP")
       
       
       
       
       rect_1D <- list(
         type ="rect",
         fillcolor =   "#e69f00",
         line = list(color = "#e69f00"),
         opacity = 0.6
       )
       
       rect_D <- list()
       
       
       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- floor(min(subset_data_2()$Log2R) - 0.5)
         rect_1D[["y1"]] <- ceiling(max(subset_data_2()$Log2R)+0.5)
         rect_D <- c(rect_D, list(rect_1D))
       }
       
       
       rect_1A <- list(
         type ="rect",
         fillcolor = "#56B4E9",
         line = list(color =  "#56B4E9"),
         opacity = 0.6
       )
       
       rect_A <-  list()
       
       for (i in c(1:dim(rects_AMP)[1])) {
         rect_1A[["x0"]] <- rects_AMP[i,]$Start
         rect_1A[["x1"]] <- rects_AMP[i,]$End
         rect_1A[["y0"]] <- floor(min(subset_data_2()$Log2R) - 0.5)
         rect_1A[["y1"]] <- ceiling(max(subset_data_2()$Log2R)+0.5)
         rect_A <- c(rect_A, list(rect_1A))
       }
       
       
       
       rect_1_2A <- list(
         type ="rect",
         fillcolor = "#0072B2",
         line = list(color = "#0072B2"),
         opacity = 0.6
       )
       
       rect_2A <- list()
       
       for (i in c(1:dim(rects_2AMP)[1])) {
         rect_1_2A[["x0"]] <- rects_2AMP[i,]$Start
         rect_1_2A[["x1"]] <- rects_2AMP[i,]$End
         rect_1_2A[["y0"]] <- floor(min(subset_data_2()$Log2R) - 0.5)
         rect_1_2A[["y1"]] <- ceiling(max(subset_data_2()$Log2R)+0.5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }
       
       
       
       
       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#d55e00",
         line = list(color = "#d55E00"),
         opacity = 0.6
       )
       
       rect_2D <- list()
       
       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- floor(min(subset_data_2()$Log2R) - 0.5)
         rect_1_2D[["y1"]] <- ceiling(max(subset_data_2()$Log2R)+0.5)
         rect_2D <- c(rect_2D, list(rect_1_2D))}
       
       
       
       rect<-c()
       
       if (dim(rects_AMP)[1] >0){
         rect <- append(rect, rect_A )
       }
       
       if (dim(rects_DEL)[1] >0){
         rect <- append(rect, rect_D) 
       }
       
       if (dim(rects_2AMP)[1] >0){
         rect <- append(rect, rect_2A) 
       }
       
       if (dim(rects_2DEL)[1] >0){
         rect <- append(rect, rect_2D) 
       }   
       
       
       
       
       
       
       
       pl_2<- plot_ly(data = subset_data_2_range(), x = ~Position) |> 
        add_trace(type = 'scatter', mode = 'markers', y = ~Log2R, color = ~Class, colors =pal, 
                  text = ~paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon),
                  hoverinfo = "text",
                   marker =list(size = 4), showlegend = FALSE)|>
         add_trace(type = 'scatter', mode = 'lines', y = ~SegMean, color = "Seg", line  = list(color = 'red'), showlegend = FALSE) |> 
          layout(shapes = rect,  yaxis = list(title = 'Log2Ratio'), xaxis = list(title = "Chromosome coordinates"))

       if (dim(rects_2DEL)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_2()$Log2R)+0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                    add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_2()$Log2R)+0.5),   color = I("#d55e00"),
                                text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                     add_trace(data = rects_2DEL,type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_2()$Log2R)-0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                     add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers',x =  ~End, y =  floor(min(subset_data_2()$Log2R)-0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                               "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00"))|> 
                    add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00"))}
     if (dim(rects_DEL)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers',x =  ~End, y =  ceiling(max(subset_data_2()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_2()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_2()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_2()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                     add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) }

       if (dim(rects_AMP)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  ceiling(max(subset_data_2()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  floor(min(subset_data_2()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  floor(min(subset_data_2()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_2()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
       if (dim(rects_2AMP)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  ceiling(max(subset_data_2()$Log2R)+0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                       add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  floor(min(subset_data_2()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))
                       add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_2()$Log2R)+0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                      add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_2()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                      add_trace(data = rects_2AMP,type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}
       if(input$Set_y_axis){
         if (is.null(input$HSLM_3) | is.null(input$FastCall_Results_3)){
           pl_2 <- pl_2|> layout(
             yaxis = list(range = c(
               min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5)), 
               max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5)))
             ))}
         else {
           pl_2 <- pl_2|> layout(
             yaxis = list(range = c(
               min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5, min(subset_data_3()$Log2R) - 0.5 )), 
               max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5))))
           )}
       }
       else{
         pl_2<- pl_2|> layout(
           yaxis = list(range = c(floor(min(subset_data_2()$Log2R) -0.5),
                                    ceiling(max(subset_data_2()$Log2R)+0.5))))}
       pl_2 <- pl_2 |>  partial_bundle()  |>  toWebGL()}
    

# Plot for third individual -----------------------------------------------

     if(!(is.null(input$HSLM_3)|is.null(input$Bed) | is.null(input$FastCall_Results_3))) {
      
       pal <- c("blue", "lightblue")
       
       pal <- setNames(pal, c("IN", "OUT"))
       
       
       
       rects_2DEL <- 
         rects_2_range() |> 
         filter(Mutation == "2-DEL")
       
       
       rects_DEL <- 
         rects_2_range() |> 
         filter(Mutation == "DEL")
       
       
       rects_AMP <- 
         rects_2_range() |> 
         filter(Mutation == "AMP")
       
       
       rects_2AMP <- 
         rects_2_range() |> 
         filter(Mutation == "2-AMP")
       
       
       
       
       
       rect_1D <- list(
         type ="rect",
         fillcolor =   "#e69f00",
         line = list(color = "#e69f00"),
         opacity = 0.6
       )
       
       rect_D <- list()
       
       
       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- floor(min(subset_data_3()$Log2R) - 0.5)
         rect_1D[["y1"]] <- ceiling(max(subset_data_3()$Log2R)+0.5)
         rect_D <- c(rect_D, list(rect_1D))
       }
       
       
       rect_1A <- list(
         type ="rect",
         fillcolor = "#56B4E9",
         line = list(color =  "#56B4E9"),
         opacity = 0.6
       )
       
       rect_A <-  list()
       
       for (i in c(1:dim(rects_AMP)[1])) {
         rect_1A[["x0"]] <- rects_AMP[i,]$Start
         rect_1A[["x1"]] <- rects_AMP[i,]$End
         rect_1A[["y0"]] <- floor(min(subset_data_3()$Log2R) - 0.5)
         rect_1A[["y1"]] <- ceiling(max(subset_data_3()$Log2R)+0.5)
         rect_A <- c(rect_A, list(rect_1A))
       }
       
       
       
       rect_1_2A <- list(
         type ="rect",
         fillcolor = "#0072B2",
         line = list(color = "#0072B2"),
         opacity = 0.6
       )
       
       rect_2A <- list()
       
       for (i in c(1:dim(rects_2AMP)[1])) {
         rect_1_2A[["x0"]] <- rects_2AMP[i,]$Start
         rect_1_2A[["x1"]] <- rects_2AMP[i,]$End
         rect_1_2A[["y0"]] <- floor(min(subset_data_3()$Log2R) - 0.5)
         rect_1_2A[["y1"]] <- ceiling(max(subset_data_3()$Log2R)+0.5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }
       

       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#d55e00",
         line = list(color = "#d55E00"),
         opacity = 0.6
       )
       
       rect_2D <- list()
       
       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- floor(min(subset_data_3()$Log2R) - 0.5)
         rect_1_2D[["y1"]] <- ceiling(max(subset_data_3()$Log2R)+0.5)
         rect_2D <- c(rect_2D, list(rect_1_2D))}
       
       
       
       rect<-c()
       
       if (dim(rects_AMP)[1] >0){
         rect <- append(rect, rect_A )
       }
       
       if (dim(rects_DEL)[1] >0){
         rect <- append(rect, rect_D) 
       }
       
       if (dim(rects_2AMP)[1] >0){
         rect <- append(rect, rect_2A) 
       }
       
       if (dim(rects_2DEL)[1] >0){
         rect <- append(rect, rect_2D) 
       }   
       
       
       
       
       
       
       
       pl_3<- plot_ly(data = subset_data_3_range(), x = ~Position)|>
        add_trace(y = ~Log2R, type = 'scatter', mode = 'markers', color = ~Class, colors =pal, 
          text = ~paste("Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>", "Exon: ", Exon),
         hoverinfo = "text",marker =list(size = 4), showlegend = FALSE)|>
        add_trace(type = 'scatter', mode = 'lines',y = ~SegMean, color = "Seg", line  = list(color = 'red'), showlegend = FALSE) |> 
        layout(shapes = rect,  yaxis = list(title = 'Log2Ratio'), xaxis = list(title = "Chromosome coordinates"))

  
  
       if(input$Set_y_axis) {
             pl_3 <- pl_3|> layout(
               yaxis = list(range = c(
                 min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5, min(subset_data_3()$Log2R) - 0.5 )), 
                 max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5))))
             )}
       else{
         pl_3 <- pl_3|> layout(
           yaxis = list(range = c(floor(min(subset_data_3()$Log2R) -0.5),
                                    ceiling(max(subset_data_3()$Log2R)+0.5))))}

       if (dim(rects_2DEL)[1] >0){
         pl_3<- pl_3|> add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_3()$Log2R)+0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_3()$Log2R)+0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                      add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_3()$Log2R)-0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_3()$Log2R)-0.5),   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00")) |> 
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#d55e00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#d55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_3 <- pl_3|> add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  ceiling(max(subset_data_3()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                      add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_3()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_3()$Log2R)-0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_3()$Log2R)+0.5),  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00")) |> 
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#e69f00"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#e69f00"))}

       if (dim(rects_AMP)[1] >0){
         pl_3 <- pl_3|> add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  ceiling(max(subset_data_3()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                   "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_3()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  floor(min(subset_data_3()$Log2R)-0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))|>
                  add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  ceiling(max(subset_data_3()$Log2R)+0.5),  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |> 
                   add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y = 0,  color = I("#56B4E9"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
       
       if (dim(rects_2AMP)[1] >0){
         pl_3<- pl_3|> add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  ceiling(max(subset_data_3()$Log2R)+0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                  "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                    add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  floor(min(subset_data_3()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                   "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                    add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  ceiling(max(subset_data_3()$Log2R)+0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                 "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                  "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                    add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  floor(min(subset_data_3()$Log2R)-0.5), color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                 "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |> 
                     add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~mid, y =  0, color = I("#0072B2"),
                                 text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                                 "Segment: " , round(Segment,2), "<br>", "CNF: ", round(CNF, 2), "<br>",
                                  "CN: ", round(CN, 2), "<br>", "Call: ", round(Call,2), "<br>",
                                   "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}
        if(input$Set_y_axis) {
          pl_3 <- pl_3|> layout(
          yaxis = list(
            min=floor(min(min(subset_data_1()$Log2R) - 0.5, min(subset_data_2()$Log2R) - 0.5, min(subset_data_3()$Log2R) - 0.5 )), 
            max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5)))
        )}
        else{
          pl_3 <- pl_3|> layout(
          yaxis = list(range =list(floor(min(subset_data_3()$Log2R) -0.5),
                                   ceiling(max(subset_data_3()$Log2R)+0.5))))}
       
       pl_3 <- pl_3 |>  partial_bundle()  |>  toWebGL()}
      
       
    if(is.null(input$HSLM_1)|is.null(input$Bed) | is.null(input$FastCall_Results_1)){return(NULL)}   
       
    else  if(is.null(input$HSLM_2) | is.null(input$FastCall_Results_2) | is.null(input$ShareAxes)){
        plt <-(plotly_build_light(pl))}
    
    else if (is.null(input$HSLM_3) | is.null(input$FastCall_Results_3)){
      if(input$GenomeBrowser){
        if(input$ShareAxes){
          plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), nrows =2, heights = c(3/4,1/4), shareX = TRUE, titleY  = TRUE)}
        else {plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), nrows =2, heights = c(3/4,1/4),titleY  = TRUE)}
      }
      else{
        if(input$ShareAxes){
          plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), nrows =2, heights = c(2/3,1/3), shareX = TRUE, titleY  = TRUE)}
        else {plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), nrows =2, heights = c(2/3,1/3),titleY  = TRUE)}
      }
       }
      
     else {
       if(input$GenomeBrowser){
          if(input$ShareAxes){
             plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), plotly_build_light(pl_3),
                          nrows =3, heights = c(3/5,1/5,1/5), shareX = TRUE,  titleY  = TRUE)}
          else {plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), plotly_build_light(pl_3),nrows =3,  heights = c(3/5,1/5, 1/5),  titleY  = TRUE)}
       }
       else{
         if(input$ShareAxes){
           plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), plotly_build_light(pl_3),
                          nrows =3, heights = c(2/4,1/4,1/4), shareX = TRUE,  titleY  = TRUE)}
         else {plt <- subplot(plotly_build_light(pl), plotly_build_light(pl_2), plotly_build_light(pl_3),nrows =3,  heights = c(2/4,1/4, 1/4),  titleY  = TRUE)}
         
       }}
    plt <- plt |> partial_bundle()  |>  toWebGL()
    
      session_store$plt <- plt
      session_store$plt

   })  %>% bindCache(input$chr, input$Prob, input$Type, input$Freq, input$slider, input$slider_Annotations, input$GenomeBrowser,
                     input$HSLM_1, input$FastCall_Results_1,input$HSLM_2, input$FastCall_Results_2, input$HSLM_3, 
                     input$FastCall_Results_3, input$AnnotSV, input$DGV_Merge,input$ShareAxes,
                     input$DGV_Gold, input$GnomAD_Genome, input$GnomAD_Exome, input$Genome, input$Set_y_axis, session
                     )
    
    

    
    output$downloadplot <- downloadHandler(
      filename = function() {
        paste(Sys.Date(), " ", input$chr," Coordinates", " ", input$slider[1],"-", input$slider[2],  ".html", sep = "")
      },
      content = function(file) {
        # export plotly html widget as a temp file to download.
        saveWidget(as_widget(session_store$plt), file, selfcontained = TRUE)
      }
    )

}  

     

# Run the application  ----------------------------------------------------

shinyApp(ui, server)
