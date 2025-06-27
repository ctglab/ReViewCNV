# Arguments ---------------------------------------------------------------

#' ReViewCNV
#'
#' @param ... Additional arguments passed to `shinyApp()`. Currently unused.
#'
#' @returns  Runs a Shiny application in the default web browser.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   Plot_Visualization()
#' }

ReViewCNV <- function(...) {


# Specify the application port
options(shiny.maxRequestSize=50*1024^2) #max dim for input files
options(shiny.host = "0.0.0.0")
options(shiny.port = 6868)
shiny::addResourcePath(prefix = 'www', directoryPath = 'inst/www')
shiny::shinyOptions(cache = cachem::cache_mem(max_size = 500e6))
options(warn = -1)


# Read variants annotations data ---------------------------------------------------

Annotations_37 <-arrow:: open_dataset("inst/37") |> dplyr::collect()
Annotations_38 <-arrow:: open_dataset("inst/38") |> dplyr::collect()



# Read genes annotation data ----------------------------------------------

genes_annotation_38 <-readRDS("inst/genes_annotation_38.rds")
genes_annotation_37 <-readRDS("inst/genes_annotation_37.rds")


# Read exons annotations --------------------------------------------------

exons_annotation_37 <- readRDS("inst/Exons_37.rds")
exons_annotation_38 <- readRDS("inst/Exons_38.rds")

# Read chromosome coordinates --------------------------------------------

hg37_Chromosomes_Coordinates <- readRDS("inst/hg37_Coordinates.rds")
hg38_Chromosomes_Coordinates <- readRDS("inst/hg38_Coordinates.rds")



# ui ----------------------------------------------------------------------


ui  <- bslib::page_sidebar(
  theme = bslib::bs_theme(
  bootswatch = "cosmo"),

 # Application title
  title ="ReViewCNV",
  #tags$style(".recalculating { opacity: inherit !important; }"),
    sidebar = bslib::sidebar(width = "35%",
      shiny:: selectInput("Genome", "Select the Genome version", selected = NULL,
        choices = c("","GRCh37", "GRCh38")),
      shiny::fileInput("FastCall_Results_1", "Load the CNV calls file"),
      shiny::fileInput("HSLM_1", "If available load the HSLM/TR level CN estimation file"),
      shiny::fileInput("Bed", "If available load the bed file of annotated targeted regions"),
      shiny::checkboxInput("GenomeBrowser", "Show genes annotations", value = FALSE, width = NULL),
      shiny::checkboxInput("ShareAxes", "Share x axis", value = FALSE, width = NULL),
      shiny::checkboxInput("Set_y_axis", "Set the same Log2R range", value = FALSE, width = NULL),
      shiny::uiOutput("Button_1"),
      shiny::uiOutput("Second_Individual"),
      shiny::uiOutput("Button_2"),
      shiny::uiOutput("Third_Individual"),
        shiny:: selectInput("chr", "Select the Chr", selected = "All",
          choices = c("All")),
      shiny::uiOutput("Button_3"),
        shiny:: selectInput("Prob", "Select Calls", selected = "All",
          choices = c("All" = "-1", "Prob Call \u2265  0.5" = "0.5", "Prob Call \u2265  0.6" = "0.6",
            "Prob Call \u2265  0.7" = "0.7", "Prob Call \u2265  0.8" = "0.8",
            "Prob Call \u2265  0.9" = "0.9")),
        shiny::h4("Choose CNVs datasets"),
        shiny::checkboxInput("AnnotSV", "AnnotSV", value = TRUE),
        shiny::checkboxInput("DGV_Merge", "DGV_Merge", value = FALSE),
        shiny::checkboxInput("DGV_Gold", "DGV_Gold", value = FALSE),
        shiny::checkboxInput("GnomAD_Genome", "GnomAD_Genome", value = FALSE),
        shiny::checkboxInput("GnomAD_Exome", "GnomAD_Exome (only GRCh38)", value = FALSE),
        shiny:: selectInput("Type", "Select CNVs type", selected = "All",
                  choices = c("All" = ".*","Gain" = "Gain",
                              "Loss" = "Loss")),
        shiny:: selectInput("Freq", "Select CNVs frequency", selected = "All",
            choices = c("All" = 0, "Freq \u003E  0.01" = 0.01,
              "Freq \u003E  0.02" = 0.02,"Freq \u003E  0.05" = 0.05,
              "Freq \u003E  0.1" = 0.1, "Freq \u003E  0.2" = 0.2,
              "Freq \u003E  0.3" = 0.3)),
        shiny::uiOutput("limits"),
        shiny::uiOutput("Annotations_limits"),
        shiny::HTML('<img src = "www/Logo.png" width = "60%" hight = "auto" >')
        ),
        shiny::uiOutput("Plot"),
        shiny::uiOutput("Download")
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






  bslib::bs_themer()
  # Storing the session for the Download Handler
  session_store <- shiny::reactiveValues()
  # This variable is updated every time there is a Download and is cached, this allows to download an unpadated plot
  rv <- shiny::reactiveValues(download_flag = 0)
  options(warn = -1)


# Coordinates -----------------------------------

  Coordinates <- shiny::reactive({
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
  }) |> shiny::bindCache(input$chr)


Coordinates_all <- shiny::reactive({
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


Chromosomes_Coordinates <-  shiny::reactive({
  if ( is.null(input$FastCall_Results_1)|input$Genome == "") {
    return(NULL)}
  else{
  Coordinates_all()|>
    mutate(level = rev(seq(0,2*dim(Coordinates_all())[1]-2, 2))) |>
    arrange(level) |>
    mutate(Chromoseome_n = stringr::str_remove(Chromosome, "chr"))}
  })



# Bed file -----------------------------------------------------------

target_data <- shiny::reactive({
  if (is.null(input$Bed)) {
    return(NULL)
  }
  else{
    target_data <-utils::read.table(input$Bed$datapath, quote="\"",fill=T,
    col.names= c("Chromosome","Start","End", "Exon"))
  }
})

# File, target data first individual and CNV all Chromosomes-----------------------------------

  file_data_1_pre <- shiny::reactive({
    if (is.null(input$HSLM_1)) {
      return(NULL)
    }
    else{
      file_data_1 <- utils::read.table(input$HSLM_1$datapath, fill=T, quote="\"", sep="\t", h = T) |>
        select(dplyr::any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |>
        mutate(dplyr::across(dplyr::where()(is.double), ~ round(.,2)))

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


file_data_1 <- shiny::reactive({
  if(!is.null(input$Bed) & ! is.null(input$HSLM_1)){
    file_data_1_merged <-
      file_data_1_pre() |>
      select(!Exon) |>
      dplyr::left_join(target_data(),by=c("Chromosome","Start", "End"))
  }
  else{file_data_1_merged <- file_data_1_pre()}
  return(file_data_1_merged)
})




fast_call_1 <- shiny::reactive({
      if (is.null(input$FastCall_Results_1)| input$Genome == "") {
        return(NULL)
        }
      else{
        fast_call_1 <- utils::read.table(input$FastCall_Results_1$datapath, header = T,
        fill=T, quote="\"")

        fast_call_1 <- fast_call_1 |>
        select(dplyr::any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
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

CNV_all_Chromosomes <- shiny::reactive({
  if (is.null(input$FastCall_Results_1)| input$Genome == "") {
    return(NULL)
  }
  else{
    fast_call_1() |>
      dplyr::left_join(Chromosomes_Coordinates(), dplyr::join_by(Chromosome))}
})


 # File and FastCall data second individual ----------------------------------

    file_data_2_pre <- shiny::reactive({
      if (is.null(input$HSLM_2)) {
        return(NULL)
      }
      else{
        file_data_2 <- utils::read.table(input$HSLM_2$datapath, fill=T, quote="\"", sep="\t", h = T) |>
          select(dplyr::any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |>
          mutate(dplyr::across(dplyr::where()(is.double), ~ round(.,2)))

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

file_data_2 <- shiny::reactive({
      if(!is.null(input$Bed) & ! is.null(input$HSLM_2)){
        file_data_2_merged <-
          file_data_2_pre() |>
          select(!Exon) |>
          dplyr::left_join(target_data(),by=c("Chromosome","Start", "End"))
      }
      else{file_data_2_merged <- file_data_2_pre()}
      return(file_data_2_merged)
    })



fast_call_2 <- shiny::reactive({
      if (is.null(input$FastCall_Results_2)| input$Genome == "" | x$val == 1) {
        return(NULL)
      }
      else{
        fast_call_2 <- utils::read.table(input$FastCall_Results_2$datapath, header = T,
          fill=T, quote="\"")

        fast_call_2 <- fast_call_2 |>
          select(dplyr::any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
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

    file_data_3_pre <- shiny::reactive({
      if (is.null(input$HSLM_3)) {
        return(NULL)
      }
      else{
        file_data_3 <- utils::read.table(input$HSLM_3$datapath, fill=T, quote="\"", sep="\t", h = T) |>
          select(dplyr::any_of(c("Chr", "Start", "End", "GC_content", "Mappability", "NRC_poolNorm", "Log2R", "SegMean", "Class", "Chromosome", "Position", "Exon"))) |>
          mutate(dplyr::across(dplyr::where()(is.double), ~ round(.,2)))

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

file_data_3 <- shiny::reactive({
  if(!is.null(input$Bed) & ! is.null(input$HSLM_3)){
    file_data_3_merged <-
      file_data_3_pre() |>
      select(!Exon) |>
      dplyr::left_join(target_data(),by=c("Chromosome","Start", "End"))
  }
  else{file_data_3_merged <- file_data_3_pre()}
  return(file_data_3_merged)
})


    fast_call_3 <- shiny::reactive({
      if (is.null(input$FastCall_Results_3)| input$Genome == "") {
        return(NULL)
      }
      else{
        fast_call_3 <- utils::read.table(input$FastCall_Results_3$datapath, header = T,
          fill=T, quote="\"")

        fast_call_3 <- fast_call_3 |>
          select(dplyr::any_of(c("Chr", "Chromosome","Start", "End", "Mutation", "CN", "Call",  "ProbCall")))
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


output$Plot =shiny::renderUI({
        if(input$chr == "All"){
          plotly::plotlyOutput("Plot_all_chr", width = "100%", height = "100%")}
        else{plotly::plotlyOutput("Plot_single_chr",  width = "100%", height = "100%")}
    })


# Plot all chromosomes ----------------------------------------------------



output$Plot_all_chr <- plotly::renderPlotly({

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
          dplyr::left_join(Chromosomes_Coordinates(), dplyr::join_by(Chromosome))





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





        fig_all <- plotly::plot_ly () |>
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
  plotly::add_trace(data = CNV, type = 'scatter', mode = 'markers', x =  ~End.x, y =  ~level.x + 0.8,
    color ="red",
    marker =list(size = 2), showlegend = FALSE,
    text = ~paste(" Chromosome: ", Chromosome,  "<br>",
      "Start: ", Start.x, "<br>", "End: ", End.x,
      "<br>", "Mutation: ", Mutation ),
    hoverinfo = "text",
    hoverlabel = list(bgcolor = "yellow"))


plt <-fig_all2 |>plotly::partial_bundle()  |>  plotly::toWebGL()

session_store$plt <- plt
session_store$plt
}
    }) |>  shiny::bindCache(input$Prob,  input$slider, input$chr,
     input$FastCall_Results_1,  rv$download_flag)


# Observe click event -----------------------------------------------------

    hover_reactive <- shiny::reactiveVal()


    shiny::observe({
      hover_data <- plotly::event_data("plotly_click")
      if (!is.null(hover_data))
        hover_reactive(hover_data)
      else{return(NULL)}
    })

    shiny::observe({
      if(!is.null(hover_reactive()) & h$val == 1){
        shiny::updateSelectInput(session,'chr',
          selected =  chromosome(),
            choices = c("All", chromosome()))

        shiny::updateSliderInput(session,
          'slider', value  =c(Start() -2000000, Start() + 1000000))
      }
    })

    shiny::observe({
      if (input$chr != "All"){
        hover_data <- NULL
        hover_reactive(hover_data)
        shiny::updateCheckboxInput(session,"GenomeBrowser", value = TRUE)
        }
    })





    chromosome <- shiny::reactive({
      if(!is.null(hover_reactive()) & input$chr == "All" & h$val == 1){

      b <- abs((hover_reactive()$y - 8.8)/2 -20)
      if (b == "23"){b = "X"}
      if (b == "24"){b = "Y"}
      c <- paste0("chr",b)
      c}
      else{return(NULL)}
    })

    Start <- shiny::reactive({
      hover_reactive()$x
    })

    shiny::observe({
      if(h$val == -1) {
        shiny::updateSelectInput(session,'chr',
        selected =  input$chr,
        choices = c( "chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
        "chr7", "chr8", "chr9", "chr10", "chr11",
        "chr12", "chr13", "chr14", "chr15", "chr16",
        "chr17", "chr18", "chr19", "chr20", "chr21",
        "chr22", "chrX"))
        }

      })

    shiny::observe({
      if(h$val == 1) {
        shiny::updateSelectInput(session,'chr',
          selected =  "All",
          choices = c( "All"))
        }
    })

    shiny::observe({
      if(input$chr == "All"){
        shiny::updateCheckboxInput(session,"GenomeBrowser", value = FALSE)}
      })



# Button3 -----------------------------------------------------------------


    h <- shiny::reactiveValues(val = 1)

output$Button_3  =shiny::renderUI({
  if(input$chr != "All"){
      if(h$val == 1){
        shiny::actionButton("Button_3", "Go to single chromosome visualization")
      }
      else {shiny::actionButton("Button_3", "Go to multiple chromosomes visualization")}
  }
  else{return(NULL)}
    })


shiny::observeEvent(input$Button_3, {h$val = h$val*(-1)})



# Subsetting file data--------------------------------------------------------------

#First individual

  subset_data_1 <- shiny::reactive({

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
   }) |> shiny::bindCache(input$chr, file_data_1())

  subset_data_1_range <- shiny::reactive({  if(is.null(subset_data_1())){return(NULL)}
    else{
    subset_data_1() |>
    filter(Start >= input$slider[1] & End <= input$slider[2])}
    })


#Second individual

  subset_data_2 <- shiny::reactive({

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
  }) |> shiny::bindCache(input$chr, file_data_2())

  subset_data_2_range <- shiny::reactive({  if(is.null(file_data_2())){return(NULL)}
    else{
      subset_data_2() |>
        filter(Start >= input$slider[1] & End <= input$slider[2])}
  })

#Third individual

  subset_data_3 <- shiny::reactive({

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
  }) |> shiny::bindCache(input$chr, file_data_3())

  subset_data_3_range <- shiny::reactive({  if(is.null(file_data_3())){return(NULL)}
    else{
      subset_data_3() |>
        filter(Start >= input$slider[1] & End <= input$slider[2])}
  })



# Subsetting FastCall data ------------------------------------------------

#First individual

rects_1 <- shiny::reactive({fast_call_1() |> filter(Chromosome == input$chr)})|>
    shiny::bindCache(input$chr,fast_call_1())


rects_1_range <- shiny::reactive ({
    if (is.null(input$slider)| input$Genome == "") {
      return(NULL)
    }
    else{rects_1_range <- rects_1() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])
    return(rects_1_range)}
  })|> shiny::bindCache(input$slider, rects_1())


#Second individual

rects_2 <- shiny::reactive({fast_call_2() |> filter(Chromosome == input$chr)})|>
  shiny::bindCache(input$chr,fast_call_2())

  rects_2_range <- shiny::reactive ({
    if (is.null(input$slider[1])| input$Genome == "") {
      return(NULL)
     }
    else { rects_2() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })|> shiny::bindCache(input$slider, rects_2())


#Third individual

  rects_3 <- shiny::reactive({fast_call_3() |> filter(Chromosome == input$chr)})|>
    shiny::bindCache(input$chr,fast_call_3())

  rects_3_range <- shiny::reactive ({
    if (is.null(input$slider[1])| input$Genome == "") {
     return(NULL)
     }
   else{rects_3() |>
      filter(Start >= input$slider[1] & End <= input$slider[2])}
   })|> shiny::bindCache(input$slider, rects_3())


# Subsetting variants annotations data ---------------------------------------------

  Annot_SV_D <- shiny::reactive({
    if(input$AnnotSV){"AnnotSV BenignSV (v. 3.4)"}
    else{return(NULL)}
  })

  DGV_Merge_D <- shiny::reactive({
    if(input$DGV_Merge){
      "dgvMerged (last updated 2020-02-25)"}
    else{return(NULL)}
  })

  DGV_Gold_D <- shiny::reactive({
      if(input$DGV_Gold){
       "dgvGold (last updated 2016-05-15)"}
    else{return(NULL)}
  })

  GnomAD_D <- shiny::reactive({
    if(input$GnomAD_Genome){
      c("gnomad_v2.1_sv.controls_only.site", "gnomad.v4.1.sv.non_neuro_controls.sites")}
    else{return(NULL)}
  })

  GnomAD_Exome_D <- shiny::reactive({
    if(input$GnomAD_Exome){
      GnomAD_Exome_D <-  "gnomad.v4.1.cnv.non_neuro_controls"}
    else{return(NULL)}
  })


  Annotations_list <-shiny::reactive({
    Annotations_list<-c(Annot_SV_D(), DGV_Gold_D(), DGV_Merge_D(), GnomAD_D(),  GnomAD_Exome_D())
    Annotations_list
  })

  Annotations_subset <-shiny::reactive({
    if(input$Genome == "GRCh37"){

      Annotations_37 |>
        filter (Chromosome == input$chr) |>
        filter(stringr::str_detect(calls, input$Type))|>
        filter(Frequency > input$Freq) |>
        filter (Database %in%  Annotations_list() ) |>
        dplyr::inner_join(rects_1() , dplyr::join_by(overlaps(Start, End, Start, End)))|>
        rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |>
        select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |>
        arrange(End_FastCall, End)}


    else if(input$Genome == "GRCh38"){
      Annotations_38 |>
        filter (Chromosome == input$chr) |>
        filter(stringr::str_detect(calls, input$Type))|>
        filter(Frequency > input$Freq) |>
        filter (Database %in% Annotations_list() ) |>
        dplyr::inner_join(rects_1(), dplyr::join_by(overlaps(Start, End, Start, End)))|>
        rename(chr = Chromosome.x, Start =Start.x, End = End.x, End_FastCall = End.y, Start_FastCall = Start.y) |>
        select("Unique_ID","ID","chr", "Start","Start_FastCall","End_FastCall", "End", "calls", "Length", "Frequency", "Database", "middle", "Frequency", "AnnotSV_Present") |>
        arrange(End_FastCall, End)}

    else {return(NULL)}

  })|> shiny::bindCache (input$Genome, input$chr, input$Type, input$Freq, Annotations_list(), fast_call_1())


  CNV1 <- shiny::reactive ({
    if(!is.null(Annotations_subset())){
    if(input$AnnotSV){
      Annotations_subset() |>  filter( AnnotSV_Present == "No")  }
    else {Annotations_subset()}}
    else{return(NULL)}
  })


  CNV_FC <- shiny::reactive({
    if(!is.null(Annotations_subset())){
      CNV1() |>
      select(Start_FastCall, End_FastCall) |>
      dplyr::distinct(Start_FastCall, End_FastCall)}
    else{return(NULL)}
      })


  CNV2 <-shiny::reactive({
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
             dplyr::slice(1) |>
             ungroup() |>
             arrange(n_overlap, dplyr::desc(round(log10(Length),0)), dplyr::desc(Frequency))
      CNV2

      }}
    else{return(NULL)}
    }) |> shiny::bindCache(CNV1(), CNV_FC())

  CNV3 <- shiny::reactive({
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
   })|> shiny::bindCache(CNV2())


# Subsetting genes annotations data ---------------------------------------

genes_annotations <-shiny::reactive ({
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
  }) |> shiny::bindCache(input$Genome, input$chr)


  rect_1_genes_annotation <- list(
    type ="line",
    line = list( color = "#008000" ))



  rect_genes_annotation <- shiny::reactive({
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

  })|> shiny::bindCache(genes_annotations())


# Subsetting exons annotations --------------------------------------------



  exons_annotations_1 <-shiny::reactive({
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
  })|> shiny::bindCache(input$Genome)


  exons_annotations <-shiny::reactive ({
    if (dim(rects_1())[1]>0){
      genes <- genes_annotations() |>
        dplyr::inner_join(rects_1() , dplyr::join_by(overlaps(Start, End, Start, End)))

      exons <- exons_annotations_1() |>
        filter (Chr == input$chr) |>
        filter(RefSeq_ID %in% genes$RefSeq_ID)
      return(exons)

    }
    else(return(NULL))

  })|> shiny::bindCache(input$chr, genes_annotations(), exons_annotations_1())


      rect_1_exons_annotation <- list(
        type ="rect",
        fillcolor = "#008000",
        opacity =0.4,
        line = list( color = "#008000" ))



      rect_exons_annotation <- shiny::reactive({
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
      }) |> shiny::bindCache(exons_annotations())



shapes <- shiny::reactive({
  c(rect_genes_annotation(), rect_exons_annotation())
  }) |> shiny::bindCache(rect_genes_annotation(), rect_exons_annotation())

# First button -------------------------------------------------------

   x <- shiny::reactiveValues(val = 1)

   shiny::observeEvent(input$Button_1, {
     x$val = x$val *(-1)
   })

   output$Button_1  =shiny::renderUI({
     if(x$val == 1){
      shiny::actionButton("Button_1", "Add new sample")
       }
     else { shiny::actionButton("Button_1", "Remove sample")}
   })


    output$Second_Individual=shiny::renderUI({
      if(x$val == -1){list(

        shiny::fileInput("FastCall_Results_2", "Load the CNV calls file"),
        shiny::fileInput("HSLM_2", "If available load the HSLM/TR level CN estimation file"))

        }
      else{return(NULL)}

     })



# Second button -----------------------------------------------------------



    z <- shiny::reactiveValues(val = 1)

    shiny::observeEvent(input$Button_2, {
      z$val = z$val *(-1)
    })


    A <- shiny::reactive({
      if (is.null(input$FastCall_Results_2)){return(NULL)}

      else if (is.null(fast_call_2())){return(NULL)}

      else if (z$val == 1){
        shiny::actionButton("Button_2", "Add new sample") }

      else { shiny::actionButton("Button_2", "Remove sample")}
    })

    output$Button_2  <-shiny::renderUI({
      K = A()
      return(K)
    })


    output$Third_Individual=shiny::renderUI({
      if(z$val == -1){list(

        shiny::fileInput("FastCall_Results_3", "Load the CNV calls file"),
        shiny::fileInput("HSLM_3", "If available load the HSLM/TR level CN estimation file"))
      }
      else{return(NULL)}

    })



# Setting the range slider for coordinates------------------------------------------------

output$limits=shiny::renderUI({
      if (is.null(input$FastCall_Results_1)|input$Genome == "" ) {
        return(NULL)}
        else if (h$val == -1 ){
          if(input$chr != "All"){
            shiny::sliderInput('slider','Choose the chromosome coordinates for the download',
            min=1,
            max=Coordinates()$End,
            value=c(1,
             Coordinates()$End))}
          else{return(NULL)}
         }
      else{
        shiny::sliderInput('slider','Choose the chromosome coordinates for the download',
        min=1,
        max=max(Coordinates_all()$End),
        value=c(1,
        max(Coordinates_all()$End)))}
    })





# Setting the range slider for annotations --------------------------------


    output$Annotations_limits=shiny::renderUI({
      if (!is.null(input$FastCall_Results_1)){
        if(!(is.null(CNV3()))){
          if (dim(CNV3())[1] >0){
            shiny::sliderInput('slider_Annotations',
            'Choose the maximum number of overlapping CNVs to visualize',
            min=0,
            max=max(CNV3()$level),
            value=c(0,
            min(20,max(CNV3()$level))))}
          }}


      else {return(NULL)}
    })



# Download ----------------------------------------------------------------


output$Download =shiny::renderUI({
  if (!is.null(input$FastCall_Results_1)){
    if (input$Genome  %in% c("GRCh37", "GRCh38")){
      shiny::downloadButton("downloadplot", "Download HTML")}
  }
  else{return(NULL)}
  }) |> shiny::bindCache(input$FastCall_Results_1, input$Genome)



# pal 1---------------------------------------------------------------------

    pal_1 <- shiny::reactive({
      if(is.null(subset_data_1_range())){return (NULL)}

      else if("IN" %in% subset_data_1_range()$Class){
        pal_1 <- c("blue", "lightblue")
        pal_1 <- stats::setNames(pal_1, c("IN", "OUT"))
      }
      else{pal_1 <- c("blue")
        pal_1 <- stats::setNames(pal_1, c("Yes"))
      }
      return(pal_1)
    })

# pal 2---------------------------------------------------------------------

    pal_2 <- shiny::reactive({
      if(is.null(subset_data_2_range())){return (NULL)}

      else if("IN" %in% subset_data_2_range()$Class){
        pal_2 <- c("blue", "lightblue")
        pal_2 <- stats::setNames(pal_2, c("IN", "OUT"))
      }
      else{pal_2 <- c("blue")
      pal_2 <- stats::setNames(pal_2, c("Yes"))
      }
      return(pal_2)
    })


# pal 3---------------------------------------------------------------------

    pal_3 <- shiny::reactive({
      if(is.null(subset_data_3_range())){return (NULL)}

      else if("IN" %in% subset_data_3_range()$Class){
        pal_3 <- c("blue", "lightblue")
        pal_3 <- stats::setNames(pal_3, c("IN", "OUT"))
      }
      else{pal_3 <- c("blue")
      pal_3 <- stats::setNames(pal_3, c("Yes"))
      }
      return(pal_3)
    })


# Plot variants --------------------------------------------------------


  output$Plot_single_chr <- plotly::renderPlotly({



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
 fig <- plotly::plot_ly () |>
   layout(shapes = rect, xaxis = list(title = "Chromosome coordinates",
       range = c(min=input$slider[1],
       max=input$slider[2])),
          yaxis = list(title = "Polymorphisms",range =list((input$slider_Annotations[1]),
          (input$slider_Annotations[2])), tickformat=',d'))|>
    plotly::add_trace(data = CNV_Gain, type = 'scatter', mode = 'markers', x =  ~End, y = ~level -0.3, color = I("blue"),
                        text = ~paste(" Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>",
                        "End:", End, "<br>", "Length: ", Length, "<br>",
                        "Reported frequency: ", Frequency),
                        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                        hoverlabel = list(bgcolor = "blue", align = "left"))|>
    plotly::add_trace(data = CNV_Loss, type = 'scatter', mode = 'markers', x = ~End, y = ~level -0.3,color = I("orange"),
                        text = ~paste(" Database:", Database, "<br>", "ID:", ID, "<br>", "Start:", Start, "<br>",
                        "End:", End, "<br>", "Length: ", Length,"<br>",
                        "Reported frequency: ", Frequency),
                        hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                        hoverlabel = list(bgcolor ="orange", align = "left")) |>
    plotly::partial_bundle()  |>  plotly::toWebGL()
  }
  }
  else{
    fig <- plotly::plot_ly (type = 'scatter', mode = "markers") |>
      layout( yaxis = list(range =list(0, 20))) |>
     plotly::partial_bundle()  |>  plotly::toWebGL()
  }
}

else{return(NULL)}



# Plot genes annotations --------------------------------------------------
if(input$GenomeBrowser){
  fig2 <- plotly::plot_ly (data = genes_annotations(), type = 'scatter', mode = 'markers', x =  ~middle_gene, y = ~level, color = I("#008000"),
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
  fig2 <- fig2 |>  plotly::add_trace(data = exons_annotations(), type = 'scatter', mode = 'markers', x =  ~Start, y = ~level + 0.3, color = I("#008000"),
    opacity =0.4,  text = ~paste(" RefSeq ID:", RefSeq_ID, "<br>", "Gene Symbol:", Gene_Symbol, "<br>",
      "Exon:",  Exon_number, "<br>","Exon Start:", Start, "<br>", "Exon End:", End, "<br>", "Direction:", Direction),
    hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
    hoverlabel = list(bgcolor = "#008000", opacity = 0.4, align = "left"))
}
        fig2 |> plotly::partial_bundle()  |>  plotly::toWebGL()
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


       pl_1<- plotly::plot_ly (data = rects_1_range(),  x =  ~End )


       if(!is.null(subset_data_1_range())){

         pl_1 <- pl_1 |>
       plotly::add_trace(data = subset_data_1_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_1(), type = 'scatter', mode = "markers",
        text = ~Text,
        hoverinfo = "text",
        marker =list(size = 4),showlegend = F)}

      pl_1 <- pl_1 |> plotly::add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "DEL",
         color= I("#E69f00"),
         opacity = 0.6) |>
      plotly::add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "AMP",
         color = I("#56B4E9"),
        opacity = 0.6)|>
      plotly::add_trace(type = 'bar',
         x = 0,
         y = 0,
         name = "2_DEL",
         color = I("#D55e00"),
         opacity = 0.6)|>
      plotly::add_trace(type = 'bar',
        x = 0,
        y = 0,
        name = "2_AMP",
        color= I("#0072B2"),
        opacity = 0.6) |>
      plotly::add_trace(type = 'bar',
        x = 0,
        y = 0,
       name = "Gain CNVs from selected databases",
       color= I("blue"))|>
    plotly::add_trace(type = 'bar',
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
       pl_1<- pl_1|> plotly::add_trace(data = rects_2DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
           "ProbCall: ", round(ProbCall,2)),
                             hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                             hoverlabel = list(bgcolor = "#D55e00"))|>
                 plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y = max(ceiling(max(subset_data_1()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
                   text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                     "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                     "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#D55e00") ) |>
                plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                  text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                    "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                    "ProbCall: ", round(ProbCall,2)),
                            hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                            hoverlabel = list(bgcolor = "#D55e00")) |>
                plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                  text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                    "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                    "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#D55e00")) |>
               plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"), opacity = 0.6,
                 text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                   "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                   "ProbCall: ", round(ProbCall,2)),
                               hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                               hoverlabel = list(bgcolor = "#D55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_1<- pl_1|> plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00") ) |>
                    plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                     plotly::add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~Start, y = min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                  plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                  plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) }
       if (dim(rects_AMP)[1] >0){
         pl_1<- pl_1|> plotly::add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                       plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                      plotly::add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                      plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    plotly::add_trace(data = rects_AMP,  type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
        if (dim(rects_2AMP)[1] >0){
         pl_1<- pl_1|> plotly::add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                         plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                          plotly::add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_1()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                            text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                              "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                              "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                          plotly::add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                            text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                              "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                              "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                           plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
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
           pl_1<-pl_1|>  plotly::add_trace(type = 'bar',
                                   x = 0,
                                   y = 0,
                                   name = "Genes annotation from MANE_RefSEq v 1.3",
                                   color= I("#008000"))}

         if (input$Genome == "GRCh37"){
           pl_1<-pl_1|>  plotly::add_trace(type = 'bar',
                                   x = 0,
                                   y = 0,
                                   name = "Genes annotation from ncbi RefSeq Select \n(last updated 2022-03-16)",
                                   color= I("#008000"))}
       }


        if(input$GenomeBrowser){

        pl<-plotly::subplot(fig2,pl_1, fig, nrows =3,heights = c(1/6,2/6,3/6),shareX = TRUE, titleY  = TRUE)}

        else{ pl<-plotly::subplot(pl_1, fig, nrows =2,heights = c(1/2,1/2),shareX = TRUE, titleY  = TRUE)}

      pl_1 <- pl_1 |> plotly::partial_bundle()  |>  plotly::toWebGL()
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



       pl_2<- plotly::plot_ly (data = rects_2_range(),  x =  ~End )


       if(!is.null(subset_data_2_range())){

         pl_2 <- pl_2 |>
           plotly::add_trace(data = subset_data_2_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_2(), type = 'scatter', mode = "markers",
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
         pl_2<- pl_2|> plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),   color = I("#D55e00"),opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                    plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),   color = I("#D55e00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                     plotly::add_trace(data = rects_2DEL,type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"),opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                     plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers',x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),   color = I("#D55e00"),opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))|>
                    plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))}
     if (dim(rects_DEL)[1] >0){
         pl_2<- pl_2|> plotly::add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers',x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#E69f00"),opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                    plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#E69f00"),opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                     plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) }

       if (dim(rects_AMP)[1] >0){
         pl_2<- pl_2|> plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}
       if (dim(rects_2AMP)[1] >0){
         pl_2<- pl_2|> plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                       plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                       plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_2()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                      plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_1()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                      plotly::add_trace(data = rects_2AMP,type = 'scatter', mode = 'markers',  x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
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
       pl_2 <- pl_2 |> plotly::partial_bundle()  |>  plotly::toWebGL()}


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




       pl_3<- plotly::plot_ly (data = rects_3_range(),  x =  ~End )


       if(!is.null(subset_data_3_range())){

         pl_3 <- pl_3 |>
           plotly::add_trace(data = subset_data_3_range(), x = ~Position, y = ~Log2R,  color = ~Class, colors =pal_3(), type = 'scatter', mode = "markers",
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
         pl_3<- pl_3|> plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                      plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00")) |>
                       plotly::add_trace(data = rects_2DEL, type = 'scatter', mode = 'markers', x =  ~mid, y =  0,   color = I("#D55e00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#D55e00"))}
       if (dim(rects_DEL)[1] >0){
         pl_3 <- pl_3|> plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                      plotly::add_trace(data = rects_DEL,  type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                        text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                          "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                          "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00")) |>
                       plotly::add_trace(data = rects_DEL, type = 'scatter', mode = 'markers', x =  ~mid, y = 0,  color = I("#E69f00"), opacity = 0.6,
                         text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                           "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                           "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#E69f00"))}

       if (dim(rects_AMP)[1] >0){
         pl_3 <- pl_3|> plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                    plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                     plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5),  color = I("#56B4E9"), opacity = 0.6,
                       text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                         "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                         "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))|>
                  plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers', x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5),  color = I("#56B4E9"), opacity = 0.6,
                    text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                      "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                      "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9")) |>
                   plotly::add_trace(data = rects_AMP, type = 'scatter', mode = 'markers',  x =  ~mid, y = 0,  color = I("#56B4E9"), opacity = 0.6,
                     text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                       "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                       "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#56B4E9"))}

       if (dim(rects_2AMP)[1] >0){
         pl_3<- pl_3|> plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~End, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
           text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
             "CN: ", CN, "<br>", "Call: ", Call, "<br>",
             "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    plotly::add_trace(data = rects_2AMP,  type = 'scatter', mode = 'markers', x =  ~End, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  max(ceiling(max(subset_data_3()$Log2R)+1.5),5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                    plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers',  x =  ~Start, y =  min(floor(min(subset_data_3()$Log2R) - 0.5), -5), color = I("#0072B2"), opacity = 0.6,
                      text = ~paste(" Start: ", Start, "<br>", "End: ", End, "<br>", "Mutation: ", Mutation,"<br>",
                        "CN: ", CN, "<br>", "Call: ", Call, "<br>",
                        "ProbCall: ", round(ProbCall,2)),
                                 hoverinfo = "text", marker =list(size = 2), showlegend = FALSE,
                                 hoverlabel = list(bgcolor = "#0072B2")) |>
                     plotly::add_trace(data = rects_2AMP, type = 'scatter', mode = 'markers', x =  ~mid, y =  0, color = I("#0072B2"), opacity = 0.6,
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

       pl_3 <- pl_3 |> plotly::partial_bundle()  |>  plotly::toWebGL()}



# Final plot --------------------------------------------------------------

    if(is.null(input$FastCall_Results_1)){return(NULL)}

    else  if(is.null(input$FastCall_Results_2) | is.null(input$ShareAxes) | x$val == 1){
        plt <-(pl)}

    else if (is.null(input$FastCall_Results_3) | z$val == 1){
      if(input$GenomeBrowser){
        if(input$ShareAxes){
          plt <- plotly::subplot(pl, pl_2, nrows =2, heights = c(3/4,1/4), shareX = TRUE, titleY  = TRUE)}
        else {plt <- plotly::subplot(pl, pl_2, nrows =2, heights = c(3/4,1/4),titleY  = TRUE)}
      }
      else{
        if(input$ShareAxes){
          plt <- plotly::subplot(pl, pl_2, nrows =2, heights = c(2/3,1/3), shareX = TRUE, titleY  = TRUE)}
        else {plt <- plotly::subplot(pl, pl_2, nrows =2, heights = c(2/3,1/3),titleY  = TRUE)}
      }
       }

     else {
       if(input$GenomeBrowser){
          if(input$ShareAxes){
             plt <- plotly::subplot(pl, pl_2, pl_3,
                          nrows =3, heights = c(3/5,1/5,1/5), shareX = TRUE,  titleY  = TRUE)}
          else {plt <- plotly::subplot(pl, pl_2, pl_3,nrows =3,  heights = c(3/5,1/5, 1/5),  titleY  = TRUE)}
       }
       else{
         if(input$ShareAxes){
           plt <- plotly::subplot(pl, pl_2, pl_3,
                          nrows =3, heights = c(2/4,1/4,1/4), shareX = TRUE,  titleY  = TRUE)}
         else {plt <- plotly::subplot(pl, pl_2, pl_3,nrows =3,  heights = c(2/4,1/4, 1/4),  titleY  = TRUE)}

       }}
    plt <- plt |>plotly::partial_bundle()  |>  plotly::toWebGL()

      session_store$plt <- plt
      session_store$plt

   } ) |>shiny::bindCache(input$chr, input$Prob, input$Type, input$Freq, input$slider, input$slider_Annotations, input$GenomeBrowser,
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

  shiny::shinyApp(ui, server)

}
