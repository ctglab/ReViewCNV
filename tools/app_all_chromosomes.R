library(shiny)
library(dplyr)
library(stringr)
library(plotly)
library(htmlwidgets)
library(bslib)
library(arrow)


# Arguments ---------------------------------------------------------------




# Specify the application port
options(shiny.maxRequestSize=50*1024^2) #max dim for input files 
options(shiny.host = "0.0.0.0")
options(shiny.port = 6868)
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

hg37_Chromosomes_Coordinates <- readRDS("hg37_Coordinates.rds")
hg38_Chromosomes_Coordinates <- readRDS("hg38_Coordinates.rds")



# ui ----------------------------------------------------------------------


ui  <- page_sidebar(
  theme = bs_theme(
  bootswatch = "cosmo"),

 # Application title
  title ="ReViewCNV",
  #tags$style(".recalculating { opacity: inherit !important; }"),
    sidebar = sidebar(width = "35%",
      selectInput("Genome", "Select the Genome version", selected = NULL,
        choices = c("","GRCh37", "GRCh38")),
      fileInput("FastCall_Results_1", "Load the CNV calls file"),
      fileInput("HSLM_1", "If available load the HSLM/TR level CN estimation file"),
      fileInput("Bed", "If available load the bed file of annotated targeted regions"),
      checkboxInput("GenomeBrowser", "Show genes annotations", value = FALSE, width = NULL),
      checkboxInput("ShareAxes", "Share x axis", value = FALSE, width = NULL),
      checkboxInput("Set_y_axis", "Set the same Log2R range", value = FALSE, width = NULL),
      uiOutput("Button_1"), 
      uiOutput("Second_Individual"),
      uiOutput("Button_2"), 
      uiOutput("Third_Individual"),
        selectInput("chr", "Select the Chr", selected = "All",
          choices = c("All")),
      uiOutput("Button_3"), 
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
        uiOutput("Plot"),
        uiOutput("Download")
        )


# Server ------------------------------------------------------------------


server <- function(input, output, session) {
#Find pandoc
  
  if (!rmarkdown::pandoc_available()) {
    alt_pandoc <- "~/local/pandoc/bin"
    alt_pandoc_bin <- file.path(alt_pandoc, "pandoc")
    
    if (file.exists(alt_pandoc_bin)) {
      message("Using user-installed Pandoc at: ", alt_pandoc)
      Sys.setenv(RSTUDIO_PANDOC = alt_pandoc)
    } else {
      warning("Pandoc not available and no fallback found. Some features may not work.")
    }
  }
  
  
  
  
  
  
  bs_themer()
  # Storing the session for the Download Handler
  session_store <- reactiveValues()
  # This variable is updated every time there is a Download and is cached, this allows to download an unpadated plot
  rv <- reactiveValues(download_flag = 0)
  options(warn = -1)

  
# Coordinates -----------------------------------
  
  Coordinates <- reactive({
    if ( is.null(input$FastCall_Results_1)|input$Genome == "") {
      return(NULL)}
    else{
      if(input$Genome == "GRCh37"){
        Coordinates <- hg37_Chromosomes_Coordinates|>
          filter(Chromosome == input$chr)
        return(Coordinates)
      }
      if(input$Genome == "GRCh38"){
        Coordinates <- hg38_Chromosomes_Coordinates|>
          filter(Chromosome == input$chr)
        return(Coordinates)}
      else{return(NULL)} }
  }) |> bindCache(input$chr)
  
  
Coordinates_all <- reactive({
    if ( is.null(input$FastCall_Results_1)|input$Genome == "") {
      return(NULL)}
    else{
      if(input$Genome == "GRCh37"){
        Coordinates <- hg37_Chromosomes_Coordinates
        return(Coordinates)
      }
      if(input$Genome == "GRCh38"){
        Coordinates <- hg38_Chromosomes_Coordinates
        return(Coordinates)}
      else{return(NULL)} }
  })
  

# Coordinates for all chromosome introductory map -------------------------

  
Chromosomes_Coordinates <-  reactive({
  if ( is.null(input$FastCall_Results_1)|input$Genome == "") {
    return(NULL)}
  else{
  Coordinates_all()|> 
    mutate(level = rev(seq(0,2*dim(Coordinates_all())[1]-2, 2))) |> 
    arrange(level) |> 
    mutate(Chromoseome_n = str_remove(Chromosome, "chr"))}
  })
  
  
  
# Bed file -----------------------------------------------------------
  
target_data <- reactive({
  if (is.null(input$Bed)) {
    return(NULL)
  }
  else{
    target_data <-read.table(input$Bed$datapath, quote="\"",fill=T,
    col.names= c("Chromosome","Start","End", "Exon"))
  }
})

# File, target data first individual and CNV all Chromosomes-----------------------------------
  
  file_data_1_pre <- reactive({
    if (is.null(input$HSLM_1)) {
      return(NULL)
    }
    else{
      file_data_1 <- read.table(input$HSLM_1$datapath, fill=T, quote="\"", sep="\t", h = T) |> 
        select(any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |> 
        mutate(across(where(is.double), ~ round(.,2)))
      
      if(is.null(file_data_1$Position)) {
        file_data_1 <- file_data_1 |> 
          mutate(Position = (Start + End)/2)}
      
      if(is.null(file_data_1$Chromosome)) {
        file_data_1 <- file_data_1 |> 
          rename("Chromosome" = "Chr")}
      
      if(is.null(file_data_1$Log2R)) {
        file_data_1 <- file_data_1 |> 
          rename("Log2R" = "NRC_poolNorm") }  
      
      if(is.null(file_data_1$Class)) {
        file_data_1 <- file_data_1|> 
          mutate(Class = "Yes") }
      
      if(is.null(file_data_1$Exon)) {
        file_data_1$Exon = NA }
      
      return(file_data_1)
    }
  })


file_data_1 <- reactive({
  if(!is.null(input$Bed) & ! is.null(input$HSLM_1)){
    file_data_1_merged <- 
      file_data_1_pre() |> 
      select(!Exon) |> 
      left_join(target_data(),by=c("Chromosome","Start", "End"))
  }
  else{file_data_1_merged <- file_data_1_pre()}
  return(file_data_1_merged)
})

  

      
fast_call_1 <- reactive({
      if (is.null(input$FastCall_Results_1)| input$Genome == "") {
        return(NULL)
        }
      else{
        fast_call_1 <- read.table(input$FastCall_Results_1$datapath, header = T,
        fill=T, quote="\"") 
        
        fast_call_1 <- fast_call_1 |> 
        select(any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
        if(is.null(fast_call_1$Mutation)){
        fast_call_1 <- fast_call_1 |> 
          mutate(Mutation = case_when(
          Call == -2 ~ "2-DEL",
          Call == -1 ~ "DEL",
          Call == 1 ~ "AMP",
          TRUE ~ "2-AMP"))}
        
        if(is.null(fast_call_1$ProbCall)){fast_call_1$ProbCall = 1}
        if(is.null(fast_call_1$CN)){fast_call_1$CN = "NA"}
        if(is.null(fast_call_1$Call)){fast_call_1$Call = "NA"} 
        
        fast_call_1 <- fast_call_1 |> 
        filter(ProbCall >= as.numeric(input$Prob)) |> 
        mutate(mid = (Start + End)/2 )
        
        names(fast_call_1) <- sub("Chr$", "Chromosome", names(fast_call_1))
        
        return(fast_call_1)
      }
    }) 

CNV_all_Chromosomes <- reactive({ 
  if (is.null(input$FastCall_Results_1)| input$Genome == "") {
    return(NULL)
  }
  else{
    fast_call_1() |> 
      left_join(Chromosomes_Coordinates(), join_by(Chromosome))}
})
    

 # File and FastCall data second individual ----------------------------------

    file_data_2_pre <- reactive({
      if (is.null(input$HSLM_2)) {
        return(NULL)
      }
      else{
        file_data_2 <- read.table(input$HSLM_2$datapath, fill=T, quote="\"", sep="\t", h = T) |> 
          select(any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |> 
          mutate(across(where(is.double), ~ round(.,2)))
        
        if(is.null(file_data_2$Position)) {
          file_data_2 <- file_data_2 |> 
            mutate(Position = (Start + End)/2)}
        
        if(is.null(file_data_2$Chromosome)) {
          file_data_2 <- file_data_2 |> 
            rename("Chromosome" = "Chr")}
        
        if(is.null(file_data_2$Log2R)) {
          file_data_2 <- file_data_2 |> 
            rename("Log2R" = "NRC_poolNorm") }  
        
        if(is.null(file_data_2$Class)) {
          file_data_2 <- file_data_2 |> 
            mutate(Class = "Yes") }
        
        if(is.null(file_data_2$Exon)) {
          file_data_2$Exon = NA }
        
        return(file_data_2)
      }
    })

file_data_2 <- reactive({
      if(!is.null(input$Bed) & ! is.null(input$HSLM_2)){
        file_data_2_merged <- 
          file_data_2_pre() |> 
          select(!Exon) |> 
          left_join(target_data(),by=c("Chromosome","Start", "End"))
      }
      else{file_data_2_merged <- file_data_2_pre()}
      return(file_data_2_merged)
    })
    
    
    
fast_call_2 <- reactive({
      if (is.null(input$FastCall_Results_2)| input$Genome == "" | x$val == 1) {
        return(NULL)
      }
      else{
        fast_call_2 <- read.table(input$FastCall_Results_2$datapath, header = T,
          fill=T, quote="\"") 
        
        fast_call_2 <- fast_call_2 |> 
          select(any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
        if(is.null(fast_call_2$Mutation)){
          fast_call_2 <- fast_call_2 |> 
            mutate(Mutation = case_when(
              Call == -2 ~ "2-DEL",
              Call == -1 ~ "DEL",
              Call == 1 ~ "AMP",
              TRUE ~ "2-AMP"))}
        
        if(is.null(fast_call_2$ProbCall)){fast_call_2$ProbCall = 1}
        if(is.null(fast_call_2$CN)){fast_call_2$CN = "NA"}
        if(is.null(fast_call_2$Call)){fast_call_2$Call = "NA"} 
        
        fast_call_2 <- fast_call_2 |> 
          filter(ProbCall >= as.numeric(input$Prob)) |> 
          mutate(mid = (Start + End)/2 )
        
        names(fast_call_2) <- sub("Chr$", "Chromosome", names(fast_call_2))
        
        return(fast_call_2)
      }
    }) 
    
    

# File and FastCall data third individual -----------------------------------
    
    file_data_3_pre <- reactive({
      if (is.null(input$HSLM_3)) {
        return(NULL)
      }
      else{
        file_data_3 <- read.table(input$HSLM_3$datapath, fill=T, quote="\"", sep="\t", h = T) |> 
          select(any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |> 
          mutate(across(where(is.double), ~ round(.,2)))
        
        if(is.null(file_data_3$Position)) {
          file_data_3 <- file_data_3 |> 
            mutate(Position = (Start + End)/2)}
        
        if(is.null(file_data_3$Chromosome)) {
          file_data_3 <- file_data_3 |> 
            rename("Chromosome" = "Chr")}
        
        if(is.null(file_data_3$Log2R)) {
          file_data_3 <- file_data_3 |> 
            rename("Log2R" = "NRC_poolNorm") }  
        
        if(is.null(file_data_3$Class)) {
          file_data_3 <- file_data_3 |> 
            mutate(Class = "Yes") }
        
        if(is.null(file_data_3$Exon)) {
          file_data_3$Exon = NA }
        
        return(file_data_3)
      }
    })

file_data_3 <- reactive({
  if(!is.null(input$Bed) & ! is.null(input$HSLM_3)){
    file_data_3_merged <- 
      file_data_3_pre() |> 
      select(!Exon) |> 
      left_join(target_data(),by=c("Chromosome","Start", "End"))
  }
  else{file_data_3_merged <- file_data_3_pre()}
  return(file_data_3_merged)
})


    fast_call_3 <- reactive({
      if (is.null(input$FastCall_Results_3)| input$Genome == "") {
        return(NULL)
      }
      else{
        fast_call_3 <- read.table(input$FastCall_Results_3$datapath, header = T,
          fill=T, quote="\"") 
        
        fast_call_3 <- fast_call_3 |> 
          select(any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
        if(is.null(fast_call_3$Mutation)){
          fast_call_3 <- fast_call_3 |> 
            mutate(Mutation = case_when(
              Call == -2 ~ "2-DEL",
              Call == -1 ~ "DEL",
              Call == 1 ~ "AMP",
              TRUE ~ "2-AMP"))}
        
        if(is.null(fast_call_3$ProbCall)){fast_call_3$ProbCall = 1}
        if(is.null(fast_call_3$CN)){fast_call_3$CN = "NA"}
        if(is.null(fast_call_3$Call)){fast_call_3$Call = "NA"} 
        
        fast_call_3 <- fast_call_3 |> 
          filter(ProbCall >= as.numeric(input$Prob)) |> 
          mutate(mid = (Start + End)/2 )
        
        names(fast_call_3) <- sub("Chr$", "Chromosome", names(fast_call_3))
        
        return(fast_call_3)
      }
    }) 
    
  
    
#  Select plot ------------------------------------------------------------
    
    
output$Plot =renderUI({
        if(input$chr == "All"){
          plotlyOutput("Plot_all_chr", width = "100%", height = "100%")}
        else{plotlyOutput("Plot_single_chr",  width = "100%", height = "100%")}
    })


# Plot all chromosomes ----------------------------------------------------

    

output$Plot_all_chr <- renderPlotly({

      if (is.null(input$FastCall_Results_1)|input$Genome == "" ) {
        return(NULL)}
      else{


        rect_1_Chromosome <- list(
          type ="polygon",
          fillcolor = "green",
          line = list( color = "green" ))


        rect_Chromosome <- list()
        for (i in c(1:dim(Chromosomes_Coordinates())[1])) {
          rect_1_Chromosome[["x0"]] <- Chromosomes_Coordinates()[i,]$Start
          rect_1_Chromosome[["x1"]] <- Chromosomes_Coordinates()[i,]$End
          rect_1_Chromosome[["y0"]] <- Chromosomes_Coordinates()[i,]$level - 0.8
          rect_1_Chromosome[["y1"]] <- Chromosomes_Coordinates()[i,]$level + 0.8
          rect_Chromosome <- c(rect_Chromosome, list(rect_1_Chromosome))
        }



        CNV <- CNV_all_Chromosomes() |>
          left_join(Chromosomes_Coordinates(), join_by(Chromosome))


        
        
        
        rect_CNV_2DEL <-
          CNV |>
          filter(Mutation == "2-DEL")
        
        
        rect_CNV_DEL <-
          CNV |>
          filter(Mutation == "DEL")
        
        
        rect_CNV_AMP <-
          CNV |>
          filter(Mutation == "AMP")
        
        
        rect_CNV_2AMP <-
          CNV |>
          filter(Mutation == "2-AMP")
        
        
        rect_2_CNV_2DEL <- list(
          type ="polygon",
          fillcolor = "yellow",
          line = list( color = "yellow" ))
        
        
        rect_2DEL <- list()
        for (i in c(1:dim(rect_CNV_2DEL)[1])) {
          rect_2_CNV_2DEL[["x0"]] <- rect_CNV_2DEL[i,]$Start.x
          rect_2_CNV_2DEL[["x1"]] <- rect_CNV_2DEL[i,]$End.x
          rect_2_CNV_2DEL[["y0"]] <- rect_CNV_2DEL[i,]$level.x - 0.8
          rect_2_CNV_2DEL[["y1"]] <- rect_CNV_2DEL[i,]$level.x + 0.8
          rect_2DEL <- c(rect_2DEL, list(rect_2_CNV_2DEL))
        }
        
        
        
        rect_2_CNV_DEL <- list(
          type ="polygon",
          fillcolor = "yellow",
          line = list( color = "yellow" ))
        
        
        rect_DEL <- list()
        for (i in c(1:dim(rect_CNV_DEL)[1])) {
          rect_2_CNV_DEL[["x0"]] <- rect_CNV_DEL[i,]$Start.x
          rect_2_CNV_DEL[["x1"]] <- rect_CNV_DEL[i,]$End.x
          rect_2_CNV_DEL[["y0"]] <- rect_CNV_DEL[i,]$level.x - 0.8
          rect_2_CNV_DEL[["y1"]] <- rect_CNV_DEL[i,]$level.x + 0.8
          rect_DEL <- c(rect_DEL, list(rect_2_CNV_DEL))
        }
        
        
        
        rect_2_CNV_AMP <- list(
          type ="polygon",
          fillcolor = "yellow",
          line = list( color = "yellow" ))
        
        
        rect_AMP <- list()
        for (i in c(1:dim(rect_CNV_AMP)[1])) {
          rect_2_CNV_AMP[["x0"]] <- rect_CNV_AMP[i,]$Start.x
          rect_2_CNV_AMP[["x1"]] <- rect_CNV_AMP[i,]$End.x
          rect_2_CNV_AMP[["y0"]] <- rect_CNV_AMP[i,]$level.x - 0.8
          rect_2_CNV_AMP[["y1"]] <- rect_CNV_AMP[i,]$level.x + 0.8
          rect_AMP <- c(rect_AMP, list(rect_2_CNV_AMP))
        }
        
        
        rect_2_CNV_2AMP <- list(
          type ="polygon",
          fillcolor = "yellow",
          line = list( color = "yellow" ))
        
        
        rect_2AMP <- list()
        for (i in c(1:dim(rect_CNV_2AMP)[1])) {
          rect_2_CNV_2AMP[["x0"]] <- rect_CNV_2AMP[i,]$Start.x
          rect_2_CNV_2AMP[["x1"]] <- rect_CNV_2AMP[i,]$End.x
          rect_2_CNV_2AMP[["y0"]] <- rect_CNV_2AMP[i,]$level.x - 0.8
          rect_2_CNV_2AMP[["y1"]] <- rect_CNV_2AMP[i,]$level.x + 0.8
          rect_2AMP <- c(rect_2AMP, list(rect_2_CNV_2AMP))
        }
        
        rect<-c(rect_Chromosome)
        
        if (dim(rect_CNV_2AMP)[1] >0){
          rect <- append(rect, rect_2AMP )
        }
        
        if (dim(rect_CNV_AMP)[1] >0){
          rect <- append(rect, rect_AMP) 
        }
        
        if (dim(rect_CNV_DEL)[1] >0){
          rect <- append(rect, rect_DEL) 
        }
        
        if (dim(rect_CNV_2DEL)[1] >0){
          rect <- append(rect, rect_2DEL) 
        }  

        



        fig_all <- plot_ly() |>
          layout(shapes = rect, 
            xaxis =list(title = "Chromosome coordinate"),
            yaxis = list(
            title = "Chromosome",
            range = c(-1,46.8),
            zeroline = FALSE,
            showline = FALSE,
            showgrid = FALSE,
            tickmode = "array",
            tickvals = c(Chromosomes_Coordinates()$level),
            ticktext = c(Chromosomes_Coordinates()$Chromoseome_n))
           )

fig_all2 <- fig_all|>
  add_trace(data = CNV, type = 'scatter', mode = 'markers', x =  ~End.x, y =  ~level.x + 0.8,
    color ="red",
    marker =list(size = 2), showlegend = FALSE,
    text = ~paste(" Chromosome: ", Chromosome,  "<br>",
      "Start: ", Start.x, "<br>", "End: ", End.x,
      "<br>", "Mutation: ", Mutation ),
    hoverinfo = "text",
    hoverlabel = list(bgcolor = "yellow"))


plt <-fig_all2 |> partial_bundle()  |>  toWebGL()

session_store$plt <- plt
session_store$plt
}
    }) |>  bindCache(input$Prob,  input$slider, input$chr,
     input$FastCall_Results_1,  rv$download_flag)
     

# Observe click event -----------------------------------------------------

    hover_reactive <- reactiveVal()   
    
    
    observe({
      hover_data <- event_data("plotly_click")
      if (!is.null(hover_data))
        hover_reactive(hover_data)  
      else{return(NULL)}
    })
    
    observe({ 
      if(!is.null(hover_reactive()) & h$val == 1){
        updateSelectInput(session,'chr',
          selected =  chromosome(),
            choices = c("All", chromosome()))
        
        updateSliderInput(session,
          'slider', value  =c(Start() -2000000, Start() + 1000000))
      }
    })
    
    observe({
      if (input$chr != "All"){
        hover_data <- NULL
        hover_reactive(hover_data)
        updateCheckboxInput(session,"GenomeBrowser", value = TRUE)
        }
    })
    
    

    
    
    chromosome <- reactive({
      if(!is.null(hover_reactive()) & input$chr == "All" & h$val == 1){
      
      b <- abs((hover_reactive()$y - 8.8)/2 -20)
      if (b == "23"){b = "X"}
      if (b == "24"){b = "Y"}
      c <- paste0("chr",b)
      c}
      else{return(NULL)}
    })
    
    Start <- reactive({
      hover_reactive()$x 
    })
    
    observe({
      if(h$val == -1) {
        updateSelectInput(session,'chr',
        selected =  input$chr,
        choices = c( "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11",
        "chr12", "chr13", "chr14", "chr15", "chr16",
        "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22", "chrX"))
        }
      
      })
    
    observe({
      if(h$val == 1) {
        updateSelectInput(session,'chr',
          selected =  "All",
          choices = c( "All"))
        }
    })
    
    observe({
      if(input$chr == "All"){
        updateCheckboxInput(session,"GenomeBrowser", value = FALSE)}  
      })
        
    

# Button3 -----------------------------------------------------------------

    
    h <- reactiveValues(val = 1) 
    
output$Button_3  =renderUI({
  if(input$chr != "All"){
      if(h$val == 1){
        actionButton("Button_3", "Go to single chromosome visualization")
      }
      else {actionButton("Button_3", "Go to multiple chromosomes visualization")}
  }
  else{return(NULL)}
    })

    
observeEvent(input$Button_3, {h$val = h$val*(-1)})


    
# Subsetting file data--------------------------------------------------------------

#First individual  

  subset_data_1 <- reactive({
    
  if(is.null(file_data_1())){return(NULL)}
  else{
    
  subset_data_1 <- file_data_1()[file_data_1()$Chromosome == input$chr,] 

  if(!is.null( subset_data_1$Mappability)){
    subset_data_1 <-  subset_data_1 |>
      mutate(Text = paste(" NRC_poolNorm: ", Log2R, "<br>", "GC content: ", GC_content, "<br>",
      "Start:", Start, "<br>", "End:", End, "<br>",  "Mappability: ", Mappability))}

  else if (!is.null(subset_data_1$SegMean)){
    subset_data_1 <-  subset_data_1 |>
      mutate(Text = paste(" Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>",
        "Exon: ", Exon,"<br>", 
        "Start:", Start, "<br>", "End:", End))}
    
  else {
    subset_data_1 <-  subset_data_1 |>
      mutate(Text = paste( " Log2R: ", Log2R, "<br>", "Start:", Start, "<br>", "End:", End))}
  return(subset_data_1)}
   }) |> bindCache(input$chr, file_data_1())
  
  subset_data_1_range <- reactive({  if(is.null(subset_data_1())){return(NULL)}
    else{
    subset_data_1() |> 
    filter(Start >= input$slider[1] & End <= input$slider[2])}
    })
  
  
#Second individual 
  
  subset_data_2 <- reactive({
    
    if(is.null(file_data_2())){return(NULL)}
    else{
      
      subset_data_2 <- file_data_2()[file_data_2()$Chromosome == input$chr,] 
      
      if(!is.null( subset_data_2$Mappability)){
        subset_data_2 <-  subset_data_2 |>
          mutate(Text = paste(" NRC_poolNorm: ", Log2R, "<br>", "GC content: ", GC_content, "<br>",
            " Start:", Start, "<br>", " End:", End, "<br>",  "Mappability: ", Mappability))}
      
      else if (!is.null(subset_data_2$SegMean)){
        subset_data_2 <-  subset_data_2 |>
          mutate(Text = paste(" Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>",
            "Exon: ", Exon,"<br>", 
            "Start:", Start, "<br>", "End:", End))}
      
      else {
        subset_data_2 <-  subset_data_2 |>
          mutate(Text = paste( " Log2R: ", Log2R, "<br>", "Start:", Start, "<br>", "End:", End))}
      return(subset_data_2)}
  }) |> bindCache(input$chr, file_data_2())
  
  subset_data_2_range <- reactive({  if(is.null(file_data_2())){return(NULL)}
    else{
      subset_data_2() |> 
        filter(Start >= input$slider[1] & End <= input$slider[2])}
  })
  
#Third individual  
   
  subset_data_3 <- reactive({
    
    if(is.null(file_data_3())){return(NULL)}
    else{
      
      subset_data_3 <- file_data_3()[file_data_3()$Chromosome == input$chr,] 
      
      if(!is.null( subset_data_3$Mappability)){
        subset_data_3 <-  subset_data_3 |>
          mutate(Text = paste(" NRC_poolNorm: ", Log2R, "<br>", "GC content: ", GC_content, "<br>",
            "Start:", Start, "<br>", "End:", End, "<br>",  "Mappability: ", Mappability))}
      
      else if (!is.null(subset_data_3$SegMean)){
        subset_data_3 <-  subset_data_3 |>
          mutate(Text = paste(" Class: ", Class, "<br>", "Log2R: ", Log2R, "<br>",
            "Exon: ", Exon,"<br>", 
            "Start:", Start, "<br>", "End:", End))}
      
      else {
        subset_data_3 <-  subset_data_3 |>
          mutate(Text = paste( " Log2R: ", Log2R, "<br>", "Start:", Start, "<br>", "End:", End))}
      return(subset_data_3)}
  }) |> bindCache(input$chr, file_data_3())
  
  subset_data_3_range <- reactive({  if(is.null(file_data_3())){return(NULL)}
    else{
      subset_data_3() |> 
        filter(Start >= input$slider[1] & End <= input$slider[2])}
  })
  
  
  
# Subsetting FastCall data ------------------------------------------------

#First individual
  
rects_1 <- reactive({fast_call_1() |> filter(Chromosome == input$chr)})|>
    bindCache(input$chr,fast_call_1())
  
  
rects_1_range <- reactive ({
    if (is.null(input$slider)| input$Genome == "") {
      return(NULL)
    }
    else{rects_1_range <- rects_1() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])
    return(rects_1_range)}
  })|> bindCache(input$slider, rects_1())


#Second individual   
  
rects_2 <- reactive({fast_call_2() |> filter(Chromosome == input$chr)})|>
  bindCache(input$chr,fast_call_2())
  
  rects_2_range <- reactive ({
    if (is.null(input$slider[1])| input$Genome == "") {
      return(NULL)
     }
    else { rects_2() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })|> bindCache(input$slider, rects_2())
   
  
#Third individual
  
  rects_3 <- reactive({fast_call_3() |> filter(Chromosome == input$chr)})|>
    bindCache(input$chr,fast_call_3())
  
  rects_3_range <- reactive ({
    if (is.null(input$slider[1])| input$Genome == "") {
     return(NULL)
     }
   else{rects_3() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })|> bindCache(input$slider, rects_3())
 
  
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
   })|> bindCache(CNV2())
  

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
  
  
  rect_1_genes_annotation <- list(
    type ="line",
    line = list( color = "#008000" ))
  
  
  
  rect_genes_annotation <- reactive({
  if(!is.null(genes_annotations())){
    rect_genes_annotation_1 <- list()
  for (i in c(1:dim(genes_annotations())[1])) {
    rect_1_genes_annotation[["x0"]] <- genes_annotations()[i,]$Start
    rect_1_genes_annotation[["x1"]] <- genes_annotations()[i,]$End
    rect_1_genes_annotation[c("y0", "y1")] <- genes_annotations()[i,]$level 
    rect_genes_annotation_1 <- c(rect_genes_annotation_1, list(rect_1_genes_annotation))
  }
    return(rect_genes_annotation_1)}
    else{return(NULL)}
    
  })|> bindCache(genes_annotations())


# Subsetting exons annotations --------------------------------------------

  
  
  exons_annotations_1 <-reactive({
    if(input$GenomeBrowser){
      if (input$Genome == "GRCh37"){
        exons <- exons_annotation_37 
        return(exons)
        
      }
      else if (input$Genome == "GRCh38"){
        exons <- exons_annotation_38
        return(exons)
      }
      else(return(NULL))}
    else(return(NULL))
  })|> bindCache(input$Genome)
  
  
  exons_annotations <-reactive ({
    if (dim(rects_1())[1]>0){
      genes <- genes_annotations() |> 
        inner_join(rects_1() , join_by(overlaps(Start, End, Start, End))) 
      
      exons <- exons_annotations_1() |> 
        filter (Chr == input$chr) |> 
        filter(RefSeq_ID %in% genes$RefSeq_ID)
      return(exons)
      
    }
    else(return(NULL))
    
  })|> bindCache(input$chr, genes_annotations(), exons_annotations_1())
  
  
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
  }) |> bindCache(rect_genes_annotation(), rect_exons_annotation())
      
# First button -------------------------------------------------------

   x <- reactiveValues(val = 1)
   
   observeEvent(input$Button_1, {
     x$val = x$val *(-1) 
   })
   
   output$Button_1  =renderUI({
     if(x$val == 1){
      actionButton("Button_1", "Add new sample")
       }
     else { actionButton("Button_1", "Remove sample")}
   })
   
   
    output$Second_Individual=renderUI({
      if(x$val == -1){list(
        
        fileInput("FastCall_Results_2", "Load the CNV calls file"),
        fileInput("HSLM_2", "If available load the HSLM/TR level CN estimation file"))
        
        }
      else{return(NULL)}
    
     })



# Second button -----------------------------------------------------------

    
    
    z <- reactiveValues(val = 1)
    
    observeEvent(input$Button_2, {
      z$val = z$val *(-1) 
    })
    
    
    A <- reactive({      
      if (is.null(input$FastCall_Results_2)){return(NULL)}
      
      else if (is.null(fast_call_2())){return(NULL)}
      
      else if (z$val == 1){
        actionButton("Button_2", "Add new sample") }
      
      else { actionButton("Button_2", "Remove sample")} 
    })
    
    output$Button_2  <- renderUI({
      K = A()
      return(K)
    }) 
    
    
    output$Third_Individual=renderUI({
      if(z$val == -1){list(
        
        fileInput("FastCall_Results_3", "Load the CNV calls file"),
        fileInput("HSLM_3", "If available load the HSLM/TR level CN estimation file"))
      }
      else{return(NULL)}
      
    })
     
    

# Setting the range slider for coordinates------------------------------------------------
   
output$limits=renderUI({
      if (is.null(input$FastCall_Results_1)|input$Genome == "" ) {
        return(NULL)}
        else if (h$val == -1 ){
          if(input$chr != "All"){
            sliderInput('slider','Choose the chromosome coordinates for the download',
            min=1,
            max=Coordinates()$End,
            value=c(1,
             Coordinates()$End))}
          else{return(NULL)}
         }
      else{
        sliderInput('slider','Choose the chromosome coordinates for the download',
        min=1,
        max=max(Coordinates_all()$End),
        value=c(1,
        max(Coordinates_all()$End)))}
    }) 
    
    
    


# Setting the range slider for annotations --------------------------------

    
    output$Annotations_limits=renderUI({
      if (!is.null(input$FastCall_Results_1)){
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
  if (!is.null(input$FastCall_Results_1)){
    if (input$Genome  %in% c("GRCh37", "GRCh38")){
      downloadButton("downloadplot", "Download HTML")}
  }
  else{return(NULL)}
  }) |> bindCache(input$FastCall_Results_1, input$Genome)
 
     

# pal 1---------------------------------------------------------------------

    pal_1 <- reactive({
      if(is.null(subset_data_1_range())){return (NULL)}
      
      else if("IN" %in% subset_data_1_range()$Class){
        pal_1 <- c("blue", "lightblue")
        pal_1 <- setNames(pal_1, c("IN", "OUT"))
      }
      else{pal_1 <- c("blue")
        pal_1 <- setNames(pal_1, c("Yes"))
      }
      return(pal_1)
    }) 

# pal 2---------------------------------------------------------------------
    
    pal_2 <- reactive({
      if(is.null(subset_data_2_range())){return (NULL)}
      
      else if("IN" %in% subset_data_2_range()$Class){
        pal_2 <- c("blue", "lightblue")
        pal_2 <- setNames(pal_2, c("IN", "OUT"))
      }
      else{pal_2 <- c("blue")
      pal_2 <- setNames(pal_2, c("Yes"))
      }
      return(pal_2)
    }) 
    
    
# pal 3---------------------------------------------------------------------
    
    pal_3 <- reactive({
      if(is.null(subset_data_3_range())){return (NULL)}
      
      else if("IN" %in% subset_data_3_range()$Class){
        pal_3 <- c("blue", "lightblue")
        pal_3 <- setNames(pal_3, c("IN", "OUT"))
      }
      else{pal_3 <- c("blue")
      pal_3 <- setNames(pal_3, c("Yes"))
      }
      return(pal_3)
    }) 

    
# Plot variants --------------------------------------------------------

    
  output$Plot_single_chr <- renderPlotly({



if(!(is.null(input$FastCall_Results_1) | is.null(input$slider_Annotations[1])|
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
          yaxis = list(title = "Polymorphisms",range =list((input$slider_Annotations[1]),
          (input$slider_Annotations[2])), tickformat=',d'))|>
    add_trace(data = CNV_Gain, type = 'scatter', mode = 'markers', x =  ~End, y = ~level -0.3, color = I("blue"),
                        text = ~paste(" Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>",
                        "End:", End, "<br>", "Length: ", Length, "<br>",
                        "Reported frequency: ", Frequency),
                        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                        hoverlabel = list(bgcolor = "blue", align = "left"))|>
    add_trace(data = CNV_Loss, type = 'scatter', mode = 'markers', x = ~End, y = ~level -0.3,color = I("orange"),
                        text = ~paste(" Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>",
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
  fig2 <- plot_ly(data = genes_annotations(), type = 'scatter', mode = 'markers', x =  ~middle_gene, y = ~level, color = I("#008000"),
          text = ~paste(" RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol,
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
    opacity =0.4,  text = ~paste(" RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol, "<br>",
      "Exon:",  Exon_number, "<br>","Exon Start:", Start, "<br>", "Exon End:", End, "<br>", "Direction:", Direction),
    hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
    hoverlabel = list(bgcolor = "#008000", opacity = 0.4, align = "left"))
}
        fig2 |>  partial_bundle()  |>  toWebGL()
  }

# Plot for first individual -----------------------------------------------


     if(!is.null(input$FastCall_Results_1)){



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
         fillcolor =   "#E69f00",
         line = list(color = "#E69f00"),
         opacity = 0.6
       )

       rect_D <- list()


       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- min(floor(min(subset_data_1()$Log2R) - 0.5), -5)
         rect_1D[["y1"]] <- max(ceiling(max(subset_data_1()$Log2R)+1.5),5)
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
         rect_1A[["y0"]] <- min(floor(min(subset_data_1()$Log2R) - 0.5), -5)
         rect_1A[["y1"]] <- max(ceiling(max(subset_data_1()$Log2R)+1.5),5)
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
         rect_1_2A[["y0"]] <- min(floor(min(subset_data_1()$Log2R) - 0.5), -5)
         rect_1_2A[["y1"]] <- max(ceiling(max(subset_data_1()$Log2R)+1.5),5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }


       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#D55e00",
         line = list(color = "#D55e00"),
         opacity = 0.6
       )

       rect_2D <- list()
       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- min(floor(min(subset_data_1()$Log2R) - 0.5), -5)
         rect_1_2D[["y1"]] <- max(ceiling(max(subset_data_1()$Log2R)+1.5),5)
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


       pl_1<- plot_ly(data = rects_1_range(),  x =  ~End )


       if(!is.null(subset_data_1_range())){

         pl_1 <- pl_1 |>
       add_trace(data = subset_data_1_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_1(), type = 'scatter', mode = "markers",
        text = ~Text,
        hoverinfo = "text",
        marker =list(size = 4),showlegend = F)}

      pl_1 <- pl_1 |> add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "DEL",
         color= I("#E69f00"),
         opacity = 0.6) |>
      add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "AMP",
         color = I("#56B4E9"),
        opacity = 0.6)|>
      add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "2_DEL",
         color = I("#D55e00"),
         opacity = 0.6)|>
      add_trace(type = 'bar',
        x = 0,
        y = 0,
        name = "2_AMP",
        color= I("#0072B2"),
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
       ) |>
    layout(shapes = rect, legend=list(title=list(text='Window')),
                      xaxis = list(title = "Chromosome coordinates", range = c(min=input$slider[1],
          max=input$slider[2])))



    if(!is.null(subset_data_1_range())){
      if(!is.null(subset_data_1_range()$Mappability)){
        pl_1 <-pl_1 |>  layout(yaxis = list(title = 'NRC poolNorm'))}
      else {
        pl_1 <-pl_1 |>  layout(yaxis = list(title = 'Log2 Ratio'))}
    }
      else {
        pl_1 <-pl_1 |>  layout(yaxis = list(
          title = "Identified CNVs",
          zeroline = FALSE,
          showline = FALSE,
          showticklabels = FALSE,
          showgrid = FALSE
        ))}




       if (dim(rects_2DEL)[1] >0){
       pl_1<- pl_1|> add_trace(data = rects_2DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
           "ProbCall: ", round(ProbCall,2)),
                             hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                             hoverlabel = list(bgcolor = "#D55e00"))|>
                 add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y = max(ceiling(max(subset_data_1()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
                   text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                     "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                     "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#D55e00") ) |>
                add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                  text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                    "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                    "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#D55e00")) |>
                add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                  text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                    "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                    "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#D55e00")) |>
               add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"), opacity = 0.6,
                 text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                   "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                   "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#D55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00") ) |>
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                     add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~Start, y = min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                  add_trace(data = rects_DEL, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                  add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) }
       if (dim(rects_AMP)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                       add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                      add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                      add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
        if (dim(rects_2AMP)[1] >0){
         pl_1<- pl_1|> add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                         add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                          add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                            text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                              "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                              "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                          add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                            text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                              "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                              "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                           add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
                             text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                               "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                               "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}




       if(input$Set_y_axis) {
         if (is.null(input$HSLM_2) | is.null(input$FastCall_Results_2)){
           pl_1 <- pl_1|> layout(
             yaxis = list(range =c(min(floor(min(subset_data_1()$Log2R) -0.5),-5)),
               max(ceiling(max(subset_data_1()$Log2R)+0.5), 5)))}
         else if (is.null(input$HSLM_3) | is.null(input$FastCall_Results_3)){
           pl_1<- pl_1|> layout(
             yaxis = list(range = c(
               min=min(floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5)),-5),
               max=max(ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5)),5)
             )))}
         else {
           pl_1 <- pl_1|> layout(
             yaxis = list(range = c(
               min=min(floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5, min(subset_data_3()$Log2R) -0.5 )),-5),
               max= max(ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5)),5)))
           )}
       }
       else{
         pl_1<- pl_1|> layout(
           yaxis = list(range =c(min(floor(min(subset_data_1()$Log2R) -0.5),-5)),
                                    max(ceiling(max(subset_data_1()$Log2R)+0.5), 5)))}

       if(input$GenomeBrowser){
         if (input$Genome == "GRCh38"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = 0,
                                   y = 0,
                                   name = "Genes annotation from MANE_RefSEq v 1.3",
                                   color= I("#008000"))}

         if (input$Genome == "GRCh37"){
           pl_1<-pl_1|>  add_trace(type = 'bar',
                                   x = 0,
                                   y = 0,
                                   name = "Genes annotation from ncbi RefSeq Select \n(last updated 2022-03-16)",
                                   color= I("#008000"))}
       }


        if(input$GenomeBrowser){

        pl<-subplot(fig2,pl_1, fig, nrows =3,heights = c(1/6,2/6,3/6),shareX = TRUE, titleY  = TRUE)}

        else{ pl<-subplot(pl_1, fig, nrows =2,heights = c(1/2,1/2),shareX = TRUE, titleY  = TRUE)}

      pl_1 <- pl_1 |>  partial_bundle()  |>  toWebGL()
      }





# Plot for the second individual ------------------------------------------


     if((!is.null(input$FastCall_Results_2) & x$val == -1)){



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
         fillcolor =   "#E69f00",
         line = list(color = "#E69f00"),
         opacity = 0.6
       )

       rect_D <- list()


       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- min(floor(min(subset_data_2()$Log2R) - 0.5), -5)
         rect_1D[["y1"]] <- max(ceiling(max(subset_data_2()$Log2R)+1.5),5)
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
         rect_1A[["y0"]] <- min(floor(min(subset_data_2()$Log2R) - 0.5), -5)
         rect_1A[["y1"]] <- max(ceiling(max(subset_data_2()$Log2R)+1.5),5)
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
         rect_1_2A[["y0"]] <- min(floor(min(subset_data_2()$Log2R) - 0.5), -5)
         rect_1_2A[["y1"]] <- max(ceiling(max(subset_data_2()$Log2R)+1.5),5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }




       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#D55e00",
         line = list(color = "#D55e00"),
         opacity = 0.6
       )

       rect_2D <- list()

       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- min(floor(min(subset_data_2()$Log2R) - 0.5), -5)
         rect_1_2D[["y1"]] <- max(ceiling(max(subset_data_2()$Log2R)+1.5),5)
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



       pl_2<- plot_ly(data = rects_2_range(),  x =  ~End )


       if(!is.null(subset_data_2_range())){

         pl_2 <- pl_2 |>
           add_trace(data = subset_data_2_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_2(), type = 'scatter', mode = "markers",
             text = ~Text,
             hoverinfo = "text",
             marker =list(size = 4),showlegend = F)}




       pl_2 <- pl_2|> layout(shapes = rect, xaxis = list(title = "Chromosome coordinates", range = c(min=input$slider[1],
            max=input$slider[2])))

       if(!is.null(subset_data_2_range())){
         if(!is.null(subset_data_2_range()$Mappability)){
           pl_2 <-pl_2 |>  layout(yaxis = list(title = 'NRC poolNorm'))}
         else {
           pl_2 <-pl_2 |>  layout(yaxis = list(title = 'Log2 Ratio'))}
       }
       else {
         pl_2 <-pl_2 |>  layout(yaxis = list(
           title = "Identified CNVs",
           zeroline = FALSE,
           showline = FALSE,
           showticklabels = FALSE,
           showgrid = FALSE
         ))}


       if (dim(rects_2DEL)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),   color = I("#D55e00"),opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                    add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),   color = I("#D55e00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                     add_trace(data = rects_2DEL,type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"),opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                     add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers',x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"),opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))|>
                    add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))}
     if (dim(rects_DEL)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers',x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#E69f00"),opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                     add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) }

       if (dim(rects_AMP)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
       if (dim(rects_2AMP)[1] >0){
         pl_2<- pl_2|> add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                       add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                       add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                      add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                      add_trace(data = rects_2AMP,type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}
       if(input$Set_y_axis){
         if (is.null(input$HSLM_3) | is.null(input$FastCall_Results_3)){
           pl_2 <- pl_2|> layout(
             yaxis = list(range = c(
               min=min(floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5)),-5),
               max=max(ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5)),5)
             )))}
         else {
           pl_2 <- pl_2|> layout(
             yaxis = list(range = c(
               min=min(floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5, min(subset_data_3()$Log2R) -0.5 )),-5),
               max= max(ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5)),5)))
           )}
       }
       else{
         pl_2<- pl_2|> layout(
           yaxis = list(range = c(min(floor(min(subset_data_2()$Log2R) -0.5), -5),
                                    max(ceiling(max(subset_data_2()$Log2R)+1.5),5))))}
       pl_2 <- pl_2 |>  partial_bundle()  |>  toWebGL()}


# Plot for third individual -----------------------------------------------

     if(!is.null(input$FastCall_Results_3) & z$val == -1) {



       rects_2DEL <-
         rects_3_range() |>
         filter(Mutation == "2-DEL")


       rects_DEL <-
         rects_3_range() |>
         filter(Mutation == "DEL")


       rects_AMP <-
         rects_3_range() |>
         filter(Mutation == "AMP")


       rects_2AMP <-
         rects_3_range() |>
         filter(Mutation == "2-AMP")




       rect_1D <- list(
         type ="rect",
         fillcolor =   "#E69f00",
         line = list(color = "#E69f00"),
         opacity = 0.6
       )

       rect_D <- list()


       for (i in c(1:dim(rects_DEL)[1])) {
         rect_1D[["x0"]] <- rects_DEL[i,]$Start
         rect_1D[["x1"]] <- rects_DEL[i,]$End
         rect_1D[["y0"]] <- min(floor(min(subset_data_3()$Log2R) - 0.5), -5)
         rect_1D[["y1"]] <- max(ceiling(max(subset_data_3()$Log2R)+1.5),5)
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
         rect_1A[["y0"]] <- min(floor(min(subset_data_3()$Log2R) - 0.5), -5)
         rect_1A[["y1"]] <- max(ceiling(max(subset_data_3()$Log2R)+1.5),5)
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
         rect_1_2A[["y0"]] <- min(floor(min(subset_data_3()$Log2R) - 0.5), -5)
         rect_1_2A[["y1"]] <- max(ceiling(max(subset_data_3()$Log2R)+1.5),5)
         rect_2A <- c(rect_2A, list(rect_1_2A))
       }


       rect_1_2D <- list(
         type ="rect",
         fillcolor = "#D55e00",
         line = list(color = "#D55e00"),
         opacity = 0.6
       )

       rect_2D <- list()

       for (i in c(1:dim(rects_2DEL)[1])){
         rect_1_2D[["x0"]] <- rects_2DEL[i,]$Start
         rect_1_2D[["x1"]] <- rects_2DEL[i,]$End
         rect_1_2D[["y0"]] <- min(floor(min(subset_data_3()$Log2R) - 0.5), -5)
         rect_1_2D[["y1"]] <- max(ceiling(max(subset_data_3()$Log2R)+1.5),5)
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




       pl_3<- plot_ly(data = rects_3_range(),  x =  ~End )


       if(!is.null(subset_data_3_range())){

         pl_3 <- pl_3 |>
           add_trace(data = subset_data_3_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_3(), type = 'scatter', mode = "markers",
             text = ~Text,
             hoverinfo = "text",
             marker =list(size = 4),showlegend = F)}



        pl_3 <- pl_3 |> layout(shapes = rect,  xaxis = list(title = "Chromosome coordinates", range = c(min=input$slider[1],
          max=input$slider[2])))


        if(!is.null(subset_data_3_range())){
          if(!is.null(subset_data_3_range()$Mappability)){
            pl_3 <-pl_3 |>  layout(yaxis = list(title = 'NRC poolNorm'))}
          else {
            pl_3 <-pl_3 |>  layout(yaxis = list(title = 'Log2 Ratio'))}
        }
        else {
          pl_3 <-pl_3 |>  layout(yaxis = list(
            title = "Identified CNVs",
            zeroline = FALSE,
            showline = FALSE,
            showticklabels = FALSE,
            showgrid = FALSE
          ))}




       if(input$Set_y_axis) {
             pl_3 <- pl_3|> layout(
               yaxis = list(range = c(
                 min=floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5, min(subset_data_3()$Log2R) -0.5 )),
                 max=ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5))))
             )}
       else{
         pl_3 <- pl_3|> layout(
           yaxis = list(range = c(floor(min(subset_data_3()$Log2R) -0.5),
                                    max(ceiling(max(subset_data_3()$Log2R)+1.5),5))))}

       if (dim(rects_2DEL)[1] >0){
         pl_3<- pl_3|> add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                      add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_3 <- pl_3|> add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                      add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00"))}

       if (dim(rects_AMP)[1] >0){
         pl_3 <- pl_3|> add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))|>
                  add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                   add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                     text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                       "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                       "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}

       if (dim(rects_2AMP)[1] >0){
         pl_3<- pl_3|> add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                     add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2"))}
        if(input$Set_y_axis) {
          pl_3 <- pl_3|> layout(
          yaxis = list(
            min=min(floor(min(min(subset_data_1()$Log2R) -0.5, min(subset_data_2()$Log2R) -0.5, min(subset_data_3()$Log2R) -0.5 )),-5),
            max= max(ceiling(max(max(subset_data_1()$Log2R) + 0.5, max(subset_data_2()$Log2R) + 0.5, max(subset_data_3()$Log2R) + 0.5)),5))
        )}
        else{
          pl_3 <- pl_3|> layout(
          yaxis = list(range =list(min(floor(min(subset_data_3()$Log2R) -0.5),5),
                                   max(ceiling(max(subset_data_3()$Log2R)+1.5),5))))}

       pl_3 <- pl_3 |>  partial_bundle()  |>  toWebGL()}



# Final plot --------------------------------------------------------------

    if(is.null(input$FastCall_Results_1)){return(NULL)}

    else  if(is.null(input$FastCall_Results_2) | is.null(input$ShareAxes) | x$val == 1){
        plt <-(pl)}

    else if (is.null(input$FastCall_Results_3) | z$val == 1){
      if(input$GenomeBrowser){
        if(input$ShareAxes){
          plt <- subplot(pl, pl_2, nrows =2, heights = c(3/4,1/4), shareX = TRUE, titleY  = TRUE)}
        else {plt <- subplot(pl, pl_2, nrows =2, heights = c(3/4,1/4),titleY  = TRUE)}
      }
      else{
        if(input$ShareAxes){
          plt <- subplot(pl, pl_2, nrows =2, heights = c(2/3,1/3), shareX = TRUE, titleY  = TRUE)}
        else {plt <- subplot(pl, pl_2, nrows =2, heights = c(2/3,1/3),titleY  = TRUE)}
      }
       }

     else {
       if(input$GenomeBrowser){
          if(input$ShareAxes){
             plt <- subplot(pl, pl_2, pl_3,
                          nrows =3, heights = c(3/5,1/5,1/5), shareX = TRUE,  titleY  = TRUE)}
          else {plt <- subplot(pl, pl_2, pl_3,nrows =3,  heights = c(3/5,1/5, 1/5),  titleY  = TRUE)}
       }
       else{
         if(input$ShareAxes){
           plt <- subplot(pl, pl_2, pl_3,
                          nrows =3, heights = c(2/4,1/4,1/4), shareX = TRUE,  titleY  = TRUE)}
         else {plt <- subplot(pl, pl_2, pl_3,nrows =3,  heights = c(2/4,1/4, 1/4),  titleY  = TRUE)}

       }}
    plt <- plt |> partial_bundle()  |>  toWebGL()

      session_store$plt <- plt
      session_store$plt

   } )  %>% bindCache(input$chr, input$Prob, input$Type, input$Freq, input$slider, input$slider_Annotations, input$GenomeBrowser,
                      input$HSLM_1, input$FastCall_Results_1,input$True_Set, input$HSLM_2, input$FastCall_Results_2,
                      input$HSLM_3, input$FastCall_Results_3, input$AnnotSV, input$DGV_Merge,input$ShareAxes,
                     input$DGV_Gold, input$GnomAD_Genome, input$GnomAD_Exome, input$Genome, input$Set_y_axis, rv$download_flag, z$val, x$val,
                      input$Bed, rects_1(), rects_1_range(), Coordinates(), Chromosomes_Coordinates(), Coordinates_all(), h$val
                     )

    
    
    output$downloadplot <- shiny::downloadHandler(
      filename = function() {
        paste(Sys.Date(), " ", input$chr," Coordinates", " ", input$slider[1],"-", input$slider[2],  ".html", sep = "")
      },
      content = function(file) {
        tryCatch({
          # export plotly html widget as a temp file to download.
          htmlwidgets::saveWidget(plotly::as_widget(shiny::isolate(session_store$plt)), file, selfcontained = TRUE)
          rv$download_flag <- rv$download_flag + 1
        },
          error = function(e) {
            message("Download error: ", e$message)
            writeLines(e$message, "/tmp/render_error.log")
          }) 
      },
      contentType = "text/html"
    )
  
} 

  shinyApp(ui, server)

