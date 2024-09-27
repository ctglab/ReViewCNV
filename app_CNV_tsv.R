library(shiny)
library(dplyr)
library(stringr)
library(plotly)
library(htmlwidgets)
library(shinyHugePlot)
library(bslib)
library(arrow)


options(shiny.maxRequestSize=50*1024^2) #max dim for input files 


# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 3838)
addResourcePath(prefix = 'www', directoryPath = 'www')
shinyOptions(cache = cachem::cache_mem(max_size = 500e6))
options(warn = -1)


# Read variants annotations data ---------------------------------------------------

Annotations_37 <- open_dataset("37") |> collect()
Annotations_38 <- open_dataset("38") |> collect()

# Read genes annotation data ----------------------------------------------

genes_annotation_38 <-readRDS("genes_annotation_38.rds")
genes_annotation_37 <-readRDS("genes_annotation_37.rds")

# Read exons annotations --------------------------------------------------

exons_annotation_37 <- readRDS("Exons_37.rds")
exons_annotation_38 <- readRDS("Exons_38.rds")


# Read chromosome coordinates --------------------------------------------

hg37_Chromosome_coordinates <- readRDS("hg37_Coordinates.rds")
hg38_Chromosome_coordinates <- readRDS("hg38_Coordinates.rds")



# ui ----------------------------------------------------------------------


ui  <- page_sidebar(
  theme = bs_theme(
  bootswatch = "cosmo"),

 # Application title
  title ="ReViewCNV",
  tags$style(".recalculating { opacity: inherit !important; }"),
    sidebar = sidebar(width = "35%",
      fileInput("CNVs1", "Load the tsv file with the CNVs"),
        selectInput("Genome", "Select the Genome version", selected = NULL,
        choices = c("","GRCh37", "GRCh38")),
      checkboxInput("GenomeBrowser", "Show genes annotations", value = FALSE, width = NULL),
      checkboxInput("ShareAxes", "Share x axis", value = TRUE, width = NULL),
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
        h4("Choose CNVs datasets"),
        checkboxInput("AnnotSV", "AnnotSV", value = TRUE),
        checkboxInput("DGV_Merge", "DGV_Merge", value = FALSE),
        checkboxInput("DGV_Gold", "DGV_Gold", value = FALSE),
        checkboxInput("GnomAD_Genome", "GnomAD_Genome", value = FALSE),
        checkboxInput("GnomAD_Exome", "GnomAD_Exome (only GRCh38)", value = FALSE),
        selectInput("Type", "Select CNVs type", selected = "All",
                  choices = c("All" = ".*","Gain" = "Gain", 
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
  # Storing the session for the Download Handler
  session_store <- reactiveValues()
  # This variable is updated every time there is a Download and is cached, this allows to download an unpadated plot
  rv <- reactiveValues(download_flag = 0)
  options(warn = -1)

  
# Coordinates -----------------------------------
  

Coordinates <- reactive({
  if (is.null(input$CNVs1)| input$Genome == "") {
    return(NULL)}
  else{
    if(input$Genome == "GRCh37"){
      Coordinates <- hg37_Chromosome_coordinates |> 
      filter(Chromosome == input$chr)
      return(Coordinates)
    }
    if(input$Genome == "GRCh38"){
      Coordinates <- hg38_Chromosome_coordinates |> 
        filter(Chromosome == input$chr)
      return(Coordinates)}
  else{return(NULL)} }
})


# CNVs  ----------------------------------------------


CNVs_1 <- reactive({
      if (is.null(input$CNVs1)| input$Genome == "") {
        return(NULL)
        }
      else{
     table <-read.table(input$CNVs1$datapath)
if(dim(table)[2] == 3){
  names(table)  = c("Chromosome", "Start", "End")
  table$Category = NA}
  
if(dim(table)[2] == 4){
  names(table)  = c("Chromosome", "Start", "End", "Category")}
table <- table|> 
mutate(mid = (Start + End)/2 )
      
        }
    })
  
  
CNVs_2 <- reactive({
    if (is.null(input$CNVs2)| input$Genome == "") {
      return(NULL)
    }
    else{
      table <-read.table(input$CNVs2$datapath)
      if(dim(table)[2] == 3){
        names(table)  = c("Chromosome", "Start", "End")
        table$Category = NA}
      
      if(dim(table)[2] == 4){
        names(table)  = c("Chromosome", "Start", "End", "Category")}
      table <- table|> 
        mutate(mid = (Start + End)/2 )
      
    }
  })


CNVs_3 <- reactive({
  if (is.null(input$CNVs3)| input$Genome == "") {
    return(NULL)
  }
  else{
    table <-read.table(input$CNVs3$datapath)
    if(dim(table)[2] == 3){
      names(table)  = c("Chromosome", "Start", "End")
      table$Category = NA}
    
    if(dim(table)[2] == 4){
      names(table)  = c("Chromosome", "Start", "End", "Category")}
    table <- table|> 
      mutate(mid = (Start + End)/2 )
    
  }
})
  

# Subsetting CNVs data ------------------------------------------------

  
rects_1 <- reactive({
    CNVs_1() |> filter(Chromosome == input$chr)})
  
  
rects_1_range <- reactive ({
    if (is.null(input$slider)| input$Genome == "") {
      return(NULL)
     }
    else{rects_1_range <- rects_1() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])
    return(rects_1_range)}
   })




rects_2 <- reactive({ 
  if(is.null(CNVs_2())){return(NULL)}
  else{
  CNVs_2() |> filter(Chromosome == input$chr)}
    })


rects_2_range <- reactive ({
  if (is.null(input$slider)| input$Genome == "" | is.null(CNVs_2())) {
    return(NULL)
  }
  else{rects_2_range <- rects_2() |>
    filter(Start >= input$slider[1] & End <= input$slider[2])
  return(rects_2_range)}
})



rects_3 <- reactive({
  CNVs_3() |> filter(Chromosome == input$chr)})


rects_3_range <- reactive ({
  if (is.null(input$slider)| input$Genome == "" | is.null(CNVs_3())) {
    return(NULL)
  }
  else{rects_3_range <- rects_3() |>
    filter(Start >= input$slider[1] & End <= input$slider[2])
  return(rects_3_range)}
})



# Setting the range slider for coordinates------------------------------------------------
   
output$limits=renderUI({
     if (is.null(input$CNVs1) | input$Genome == ""){
       return(NULL)}
    else{
       
       sliderInput('slider','Choose the chromosome coordinates for the download',
                   min=1,
                   max=Coordinates()$End,
                   value=c(1,
                   Coordinates()$End))}
   }) 
   

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
})|> bindCache(input$Genome, input$chr)


rect_1_genes_annotation <- list(
  type ="line",
  line = list( color = "#008000" ))



rect_genes_annotation <- reactive({
  rect_genes_annotation_1 <- list()
  for (i in c(1:dim(genes_annotations())[1])) {
    rect_1_genes_annotation[["x0"]] <- genes_annotations()[i,]$Start
    rect_1_genes_annotation[["x1"]] <- genes_annotations()[i,]$End
    rect_1_genes_annotation[c("y0", "y1")] <- genes_annotations()[i,]$level 
    rect_genes_annotation_1 <- c(rect_genes_annotation_1, list(rect_1_genes_annotation))
  }
  return(rect_genes_annotation_1)
  
})|> bindCache(genes_annotations())


# Subsetting exons annotations --------------------------------------------

exons_annotations <-reactive ({
  if(input$GenomeBrowser){
    
    if (input$Genome == "GRCh37"){
      if (dim(rects_1())[1]>0){
        
        genes <- genes_annotation_37 |>  
          filter (Chr == input$chr) |> 
          inner_join(rects_1() , join_by(overlaps(Start, End, Start, End))) 
        
        exons <- exons_annotation_37 |> 
          filter(RefSeq_ID %in% genes$RefSeq_ID)
        
        return(exons)
        
      }
      else(return(NULL))}
    
    else if (input$Genome == "GRCh38"){
      if (dim(rects_1())[1]>0){
        
        genes <- genes_annotation_38 |>  
          filter (Chr == input$chr) |> 
          inner_join(rects_1() , join_by(overlaps(Start, End, Start, End))) 
        
        exons <- exons_annotation_38 |> 
          filter(RefSeq_ID %in% genes$RefSeq_ID)
        
        return(exons)
        
      }
      else(return(NULL))}
  }
  else(return(NULL))
}) |> bindCache(input$Genome, input$chr)


rect_1_exons_annotation <- list(
  type ="rect",
  fillcolor = "#008000",
  opacity =0.4,
  line = list( color = "#008000" ))



rect_exons_annotation <- reactive({
  if(!is.null(exons_annotations())){
    if(dim(exons_annotations())[1] >0){
      
      rect_exons_annotation_1<- list()
      for (i in c(1:dim(exons_annotations())[1])) {
        rect_1_exons_annotation[["x0"]] <- exons_annotations()[i,]$Start
        rect_1_exons_annotation[["x1"]] <- exons_annotations()[i,]$End
        rect_1_exons_annotation[["y0"]] <-  exons_annotations()[i,]$level + 0.3
        rect_1_exons_annotation[["y1"]] <-  exons_annotations()[i,]$level - 0.3
        rect_exons_annotation_1 <- c(rect_exons_annotation_1, list(rect_1_exons_annotation))
      }
      return(rect_exons_annotation_1)
    }
    else{return(NULL)}}
  else{return(NULL)}
}) |> bindCache(exons_annotations())



shapes <- reactive({
  c(rect_genes_annotation(), rect_exons_annotation())
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
      filter (Chromosome == input$chr) |> 
      filter(str_detect(calls, input$Type))|> 
      filter(Frequency > input$Freq) |> 
      filter (Database %in%  Annotations_list() ) |> 
      inner_join(rects_1() , join_by(overlaps(Start, End, Start, End)))|>
      rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |> 
      select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |> 
      arrange(End_FastCall, End)}
  
  
  else if(input$Genome == "GRCh38"){
    Annotations_38 |>
      filter (Chromosome == input$chr) |> 
      filter(str_detect(calls, input$Type))|>  
      filter(Frequency > input$Freq) |> 
      filter (Database %in% Annotations_list() ) |> 
      inner_join(rects_1(), join_by(overlaps(Start, End, Start, End)))|>
      rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |> 
      select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |> 
      arrange(End_FastCall, End)}
  
  else {return(NULL)}
  
})  %>% bindCache (input$Genome, input$chr, input$Type, input$Freq, Annotations_list())


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
})|> bindCache(CNV2())



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
    
    
# Setting the range slider for annotations --------------------------------
    
    
    output$Annotations_limits=renderUI({
      if (!(is.null(input$CNVs1))){
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
  if(x$val == 1){
    fileInput("CNVs2", "Load the tsv file with the CNVs")
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
  if(z$val == 2){
    fileInput("CNVs3", "Load the tsv file with the CNVs")
  }
  else{return(NULL)}
  
})


# Download button----------------------------------------------------------------


output$Download = renderUI({
  if (!(is.null(input$CNVs1))){
    if (input$Genome  != ""){
      downloadButton("downloadplot", "Download HTML")}
  }
  else{return(NULL)}
})




# Plot --------------------------------------------------------------------



output$zoomPlot <- renderPlotly({
  
  

# Plot variants -----------------------------------------------------------

if(!(is.null(input$CNVs1)| is.null(input$slider_Annotations[1])|
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
            range = c(min=input$slider[1], 
              max=input$slider[2])),
            yaxis = list(title = "Polymorphisms",  range =list((input$slider_Annotations[1]),
              (input$slider_Annotations[2]))))|>
          add_trace(data = CNV_Gain, type = 'scatter', mode = 'markers', x =  ~End, y = ~level -0.3, color = I("blue"),
            text = ~paste("Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>", 
              "End:", End, "<br>", "Length: ", Length, "<br>",
              "Reported frequency: ", Frequency),
            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
            hoverlabel = list(bgcolor = "blue", align = "left"))|>
          add_trace(data = CNV_Loss, type = 'scatter', mode = 'markers', x = ~End, y = ~level -0.3,color = I("orange"),
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
        layout( yaxis = list(title = "Polymorphisms", range =list(0, 20), showticklabels=FALSE)) |> 
        partial_bundle()  |>  toWebGL()
    }
  }
  
  else{return(NULL)}
  
# Plot genes annotations --------------------------------------------------
  if(input$GenomeBrowser){
    
fig2 <- plot_ly(data = genes_annotations(), type = 'scatter', mode = 'markers', x =  ~middle_gene, y = ~level, color = I("#008000"),
      text = ~paste("RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol, 
        "<br>","Start:", Start, "<br>", "End:", End),
      hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
      hoverlabel = list(bgcolor = "#008000" ,align = "left")) |> 
      layout(shapes = shapes(), xaxis = list(title = "Chromosome coordinates",
        range = c(min=input$slider[1], 
          max=input$slider[2])),
        yaxis = list( showline=  F, tickvals = c(1, 2), title = "Genes",
          ticktext = c('+', '\u2212'),  tickfont = list(size = 20, family = "Arial black"))) 
    
    if(!is.null(exons_annotations())){
      fig2 <- fig2 |>  add_trace(data = exons_annotations(), type = 'scatter', mode = 'markers', x =  ~Start, y = ~level + 0.3, color = I("#008000"),
        opacity =0.4,  text = ~paste("RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol, "<br>",
          "Exon:",  Exon_number, "<br>","Exon Start:", Start, "<br>", "Exon End:", End, "<br>", "Direction:", Direction),
        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
        hoverlabel = list(bgcolor = "#008000", opacity = 0.4, align = "left"))
    }
fig2 |>  partial_bundle()  |>  toWebGL()
  }
    
  
# Plot first individual ---------------------------------------------------


if(!(is.null(rects_1_range()))){
         
       rect  <- list(
         type ="rect",
         fillcolor =   "#e69f00",
         line = list(color = "#e69f00"),
         opacity = 0.6
       )
       
       rect_plot <- list()
       
       
       for (i in c(1:dim(rects_1_range())[1])) {
         rect[["x0"]] <- rects_1_range()[i,]$Start
         rect[["x1"]] <- rects_1_range()[i,]$End
         rect[["y0"]] <- 0
         rect[["y1"]] <- 1
         rect_plot <- c(rect_plot, list(rect))
       }
       
       rect_plot2 <-c()
       if(dim(rects_1_range())[1]>0){
         rect_plot2 <- rect_plot}
      
       
       pl_1<- plot_ly(data = rects_1_range(),  x =  ~End ) |> 
         layout(shapes = rect_plot2, legend=list(title=list(text='Window')),
           yaxis = list(title = 'CNVs',  showticklabels=FALSE), xaxis = list(title = "Chromosome coordinates",  range = c(min=input$slider[1], 
             max=input$slider[2]))) |> 
      add_trace(type = 'bar', 
         x = 0, 
         y = 0, 
         name = "CNV", 
         color= I("#E69F00"),
        opacity = 0.6) |> 
      add_trace(type = 'bar', 
        x = 0, 
        y = 0, 
       name = "Gain CNVs from selected databases", 
       color= I("blue"))|>
    add_trace(type = 'bar', 
      x = 0, 
      y = 0, 
     name = "Loss CNVs from selected databases", 
     color= I("orange")
       ) 
       
  
      
  

       if (dim(rects_1_range())[1] >0){
       pl_1 <- pl_1|> add_trace(data = rects_1_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  1,   color = I("#e69f00"),
                             text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
                             hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                             hoverlabel = list(bgcolor = "#e69f00"))|>
         add_trace(data = rects_1_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  1,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|> 
         add_trace(data = rects_1_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  0,   color = I("#e69f00"),
           text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
           hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
           hoverlabel = list(bgcolor = "#e69f00"))|>
         add_trace(data = rects_1_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  0,   color = I("#e69f00"),
           text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
           hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
           hoverlabel = list(bgcolor = "#e69f00"))|>
         add_trace(data = rects_1_range(),  type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#e69f00"),
           text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
           hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
           hoverlabel = list(bgcolor = "#e69f00"))}

       

       
       if(input$GenomeBrowser){
         if (input$Genome == "GRCh38"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = input$slider[1],
                                   y = 0,
                                   name = "Genes annotation from MANE_RefSEq v 1.3",
                                   color= I("#008000"))}
         
         if (input$Genome == "GRCh37"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = input$slider[1],
                                   y = 0,
                                   name = "Genes annotation from ncbi RefSeq Select \n(last updated 2022-03-16)",
                                   color= I("#008000"))}
       }
       
       
       
      pl_1  <- pl_1 |>  partial_bundle()  |>  toWebGL()
      
     if(input$GenomeBrowser){
        
        pl<-subplot(plotly_build_light(fig2),plotly_build_light(pl_1),plotly_build_light(fig), nrows =3,heights = c(1/6,2/6,3/6),shareX = TRUE, titleY  = TRUE)}
      
      else{ pl<-subplot(plotly_build_light(pl_1),plotly_build_light(fig), nrows =2,heights = c(1/2,1/2),shareX = TRUE, titleY  = TRUE)}
}
  
    


  

# Plot second individual --------------------------------------------------

if(!(is.null(rects_2_range()))){
    
    rect  <- list(
      type ="rect",
      fillcolor =   "#e69f00",
      line = list(color = "#e69f00"),
      opacity = 0.6
    )
    
    rect_plot <- list()
    
    
    for (i in c(1:dim(rects_2_range())[1])) {
      rect[["x0"]] <- rects_2_range()[i,]$Start
      rect[["x1"]] <- rects_2_range()[i,]$End
      rect[["y0"]] <- 0
      rect[["y1"]] <- 1
      rect_plot <- c(rect_plot, list(rect))
    }
    
    rect_plot2 <-c()
    if(dim(rects_2_range())[1]>0){
      rect_plot2 <- rect_plot}
      
    
    
    pl_2 <- plot_ly(data = rects_2_range(),  x =  ~End ) |> 
      layout(shapes = rect_plot2, legend=list(title=list(text='Window')),
        yaxis = list(title = 'CNVs',  showticklabels=FALSE), xaxis = list(title = "Chromosome coordinates",  range = c( min=input$slider[1], 
          max=input$slider[2]))) 
    
    
    if (dim(rects_2_range())[1] >0){
      pl_2 <- pl_2|> add_trace(data = rects_2_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  1,   color = I("#e69f00"),
        text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
        hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_2_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  1,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|> 
        add_trace(data = rects_2_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_2_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_2_range(),  type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))}
    
    
    
    
  pl2  <- pl_2 |>  partial_bundle()  |>  toWebGL()}
 
  

# Plot third individual ---------------------------------------------------

if(!(is.null(rects_3_range()))){
    
    rect  <- list(
      type ="rect",
      fillcolor =   "#e69f00",
      line = list(color = "#e69f00"),
      opacity = 0.6
    )
    
    rect_plot <- list()
    
    
    for (i in c(1:dim(rects_3_range())[1])) {
      rect[["x0"]] <- rects_3_range()[i,]$Start
      rect[["x1"]] <- rects_3_range()[i,]$End
      rect[["y0"]] <- 0
      rect[["y1"]] <- 1
      rect_plot <- c(rect_plot, list(rect))
    }
    
    rect_plot2 <-c()
    if(dim(rects_3_range())[1]>0){
      rect_plot2 <- rect_plot} 
      
    
    
    pl_3 <- plot_ly(data = rects_3_range(),  x =  ~End ) |> 
      layout(shapes = rect_plot2, legend=list(title=list(text='Window')),
        yaxis = list(title = 'CNVs',  showticklabels=FALSE), xaxis = list(title = "Chromosome coordinates",  range = c( min=input$slider[1], 
          max=input$slider[2]))) 
    
    
    if (dim(rects_3_range())[1] >0){
      pl_3 <- pl_3|> add_trace(data = rects_3_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  1,   color = I("#e69f00"),
        text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
        hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_3_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  1,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|> 
        add_trace(data = rects_3_range(),  type = 'scatter', mode = 'markers', x =  ~Start, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_3_range(),  type = 'scatter', mode = 'markers', x =  ~End, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))|>
        add_trace(data = rects_3_range(),  type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#e69f00"),
          text = ~paste("Start: ", Start, "<br>", "End: ", End, "<br>", "Category: ", Category),
          hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
          hoverlabel = list(bgcolor = "#e69f00"))}
    
    
  
    pl3  <- pl_3 |>  partial_bundle()  |>  toWebGL()}


# Final plot --------------------------------------------------------------

  
  if(is.null(input$CNVs1)){return(NULL)}   
  
  else  if(is.null(input$CNVs2) |  is.null(input$ShareAxes)){
    plt <-plotly_build_light(pl)}
  
  else if (is.null(input$CNVs3)){
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
  
})  %>% bindCache(input$chr, input$Type, input$Freq, input$slider, input$slider_Annotations, input$GenomeBrowser,
  input$CNVs1, input$CNVs2, input$CNVs3, input$AnnotSV, input$DGV_Merge,input$ShareAxes,
  input$DGV_Gold, input$GnomAD_Genome, input$GnomAD_Exome, input$Genome, input$Set_y_axis, rv$download_flag)
    
    output$downloadplot <- downloadHandler(
      filename = function() {
        paste(Sys.Date(), " ", input$chr," Coordinates", " ", input$slider[1],"-", input$slider[2],  ".html", sep = "")
      },
      content = function(file) {
        # export plotly html widget as a temp file to download.
        saveWidget(as_widget(isolate(session_store$plt)), file, selfcontained = TRUE)
        rv$download_flag <- rv$download_flag + 1
      }
    )

}  

shinyApp(ui, server)

