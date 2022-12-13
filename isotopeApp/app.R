library(shiny)
library(shinyFiles)
library(dplyr)
library(xcms)
library(stringr)
library(gtools)
library(MSnbase)
library(purrr)
library(ggplot2)
library(bslib)
library(future.apply)
library(thematic)
library(shinycssloaders)
source("./helpers.R")
options(warn = -1)

# Define UI ----
ui <- navbarPage(

  title = HTML("v1.0"),
  
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  tabPanel(
    "Isotope Abundance Calculator",
    sidebarLayout(
      # Side Panel ----
      sidebarPanel(
        shinyDirButton('datadir', label = "Choose Data Directory", title = "Please select data path", multiple = FALSE, width = "100%"),
        br(),
        br(),
        uiOutput("centCheck"),
        conditionalPanel(
          condition = "output.showPanel",
          helpText("Load Data when ready"),
          actionButton("loadData", "Load Data")),
        conditionalPanel(
          condition = "input.loadData",
          helpText("Enter compound formula"),
          fluidRow(
            column(6,numericInput("C","C",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("H","H",value = NULL,min=0,max=60,step=1))
          ),
          fluidRow(
            column(6,numericInput("N","N",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("O","O",value = NULL,min=0,max=60,step=1))
          ),
          helpText("Enter mass error in ppm"),
          numericInput("ppm",label = NULL, value = 5, min = 0, max = 50, step = 1, width = "50%"),
          helpText("Choose data polarity"),
          selectInput("polarity", label = NULL, choices = c("+","-"), width = "50%"),
          helpText("Choose compound labeling status"),
          selectInput("unlabeled", label = NULL, choices = c("[13]C and [15]N labeled", "[13]C labeled", "[15]N labeled")),
          helpText("Fill in compound formula to proceed"),
        ),
        conditionalPanel(
          condition = "output.showPlot",
          helpText("Peak Range (sec);"),
          uiOutput("rt_control"),
          helpText("Use intensites greater than what percent of max intensity to compute average?"),
          numericInput("multiplier", label = NULL, value = 1, min = 0.1, max = 1, step = 0.1, width = "50%"),
          checkboxInput("background", label = "Background subtraction (beta)", value = TRUE)
        ),
        conditionalPanel(
          condition = "input.background && output.showPlot",
          helpText("Background Range (sec);"),
          uiOutput("bg_control")
        ),
        conditionalPanel(
          condition = "output.showPlot",
          helpText("Plot chromatogram"),
          actionButton("plotData", "Plot"),
          br(),
          br(),
          helpText("Run when done setting parameters"),
          actionButton("run", "Run")
        ),
        conditionalPanel(
          condition = "output.showDownload",
          br(),
          helpText("Click to save result"),
          downloadButton("download","Save")
        )
      ),
      
      # Main Panel ----
      mainPanel(
        "Files to be processed:",
        verbatimTextOutput("test"),
        br(),
        htmlOutput("filelist"),
        htmlOutput("readDone"),
        conditionalPanel(
          condition = "input.plotData",
          tableOutput("mztable"),
          uiOutput("plotOptions"),
          withSpinner(plotOutput("chr", brush = "chr.brush"), type = 5 ,color = "#18BC9C", size = .65)
        ),
        conditionalPanel(
          condition = "input$run",
          withSpinner(tableOutput("table"), type = 5 ,color = "#18BC9C", size = .65)
        )
      )
    )
  ),
  
  # Help Tab ----
  tabPanel(
    HTML("Help"),
    
    includeHTML("./www/help.html")
  )
)

# Define server logic ----
server <- function(input, output) {
  # Root folder
  roots <- c(data = '../data')
  
  # Connect to root folder
  shinyDirChoose(input, "datadir", roots = roots, allowDirCreate = TRUE)
  
  # Make reactive file list
  files <- reactive({
    req(input$datadir, !isEmpty(input$centroid))
    fls <- list.files(path = parseDirPath(roots, input$datadir), pattern = "^.*\\.mzML$",
                 all.files = TRUE, full.names = TRUE,
                 recursive = FALSE, ignore.case = FALSE,
                 include.dirs = FALSE, no.. = TRUE)  %>%
      mixedsort(decreasing = F)
    if (input$centroid) {
      fls_cent <- fls[str_detect(fls,"cent")]
      if (length(fls_cent)==0) {
        return(fls)
      } else if (length(fls_cent)>0) {
        return(fls_cent)
      }
    } else if (!input$centroid) {
      fls <- fls[!str_detect(fls, "cent")]
      fls
    }
  })
  
  # Render centroid checkbox
  output$centCheck <- renderUI({
    req(input$datadir)
    checkboxInput("centroid", label = "Use centroided data", value = TRUE)
  })
  
  # Print filelist
  output$filelist <- renderUI({
    req(files(), !isEmpty(input$centroid))
    if (input$centroid == TRUE) {
      check <- files() %>% str_detect(.,"cent") %>% sum
      if (check>0) {
        out <- HTML(paste(c(files() %>% basename), collapse = '<br/>' ))
      } else if (check == 0) {
        out <- HTML(paste(c(files() %>% basename,"Centroids will be computed on the fly."), collapse = '<br/>'))
      }
    } else if (input$centroid == FALSE) {
      out <- HTML(paste(files() %>% basename, collapse = '<br/>' ))
    }
    tagList(
      out,
      HTML("<hr/>")
    )
  })
  
  # Condition for showing other controls
  output$showPanel <- reactive({
    req(files())
    TRUE
  })
  outputOptions(output, "showPanel", suspendWhenHidden = FALSE)
  
  # Condition for showing plot button
  output$showPlot <- reactive({
    req(input$C,input$H,input$N,input$O,data())
    TRUE
  })
  outputOptions(output, "showPlot", suspendWhenHidden = FALSE)
  
  # Read reactive data item
  data <- reactive({
    read_data(fls = files(),centroidData = input$centroid)
  }) %>%
    bindCache(files(), input$centroid) %>%
    bindEvent(input$loadData)
  
  dataNames <- reactive({
    req(data())
    paste(fileNames(data())) %>% basename
  })
  
  output$readDone <- renderUI({
    tagList(
      HTML(paste(dataNames(), "reading complete.",collapse = "<br/>")),
      HTML("<hr/>")
    )
  })
  
  # Generate reactive plot options
  pol <- reactive({
    switch(
      input$polarity,
      "+" = 1,
      "-" = -1
    )
  })
  
  mzs <- reactive({
    req(input$C, input$N, input$H, input$O, pol(), data())
    mzs <- get_mz_cn(C=input$C,N=input$N,H=input$H,O=input$O,pol=pol())
    
    mzs <- switch(
      input$unlabeled,
      "[13]C and [15]N labeled" = mzs, 
      "[13]C labeled" = mzs[nrow(mzs),], 
      "[15]N labeled" = select(mzs,ncol(mzs))
    )
    
    return(mzs)
    }) %>%
    bindCache(input$C, input$N, input$H, input$O, pol(), input$unlabeled)
  
  label <- reactive({
    req(mzs())
    sapply(colnames(mzs()), function(x) {
      sapply(rownames(mzs()), function(y) {
        paste(x,y,sep="_")
      })
    })
  }) %>%
    bindCache(mzs())
  
  mzs_vec <- reactive({
    req(mzs())
    as.vector(as.matrix(mzs()))
  }) %>%
    bindCache(mzs())
  
  output$mztable <- renderTable({
    mzs()
  }, rownames = TRUE, digits = 4, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$plotOptions <- renderUI({
    fluidRow(
      column(3,selectInput("whichmz", label = "Choose mz value", choices = round(mzs_vec(),4))),
      column(3, selectInput("whichfile", label = "Choose file to plot", choices = c(dataNames(),"All"))),
      column(6, selectInput("whichrange", label = "Select for", choices = c("Peak", "Background")))
    )
  })
  
  # Generate reactive rt and background controls
  rtRange <- reactive({
    req(data())
    data() %>% filterFile(1) %>% fData %>% select(retentionTime) %>% c(min(.),max(.)) %>% `[`(-1) %>% unlist %>% round(.,3)
  }) %>%
    bindCache(data())
  
  output$rt_control <- renderUI({
    sliderInput("rt_select",label = NULL, min = rtRange()[1], max = rtRange()[2], value = rtRange(), step = 0.001, round = FALSE, ticks = FALSE, dragRange = TRUE, animate = FALSE)
  })
  
  output$bg_control <- renderUI({
    sliderInput("bg_select",label = NULL, min = rtRange()[1], max = rtRange()[2], value = rtRange(), step = 0.001, round = FALSE, ticks = FALSE, dragRange = TRUE, animate = FALSE)
  })
  
  # Generate plot
  thematic::thematic_shiny()
  chr <- reactive({
    req(data(),input$ppm, input$whichmz)
    data() %>% filterMz(mz = get_ppm_range(x = as.numeric(input$whichmz), ppm = input$ppm)) %>% chromatogram(aggregationFun = "max", missing = 0)
  }) %>%
    bindCache(data(),input$ppm, input$whichmz)
  
  output$chr <- renderPlot({
    req(chr(),input$whichfile, dataNames())
    chrom_data <- lapply(chr(),function(x){
      cbind(x@rtime,x@intensity) %>%
        as.data.frame()
    })
    groups <- sub(chr()@phenoData@data[,1], pattern = "\\.mzX?ML", replacement = "")
    
    if (input$whichfile != "All") {
      idx <- which(dataNames() == input$whichfile)
      chrom_data <- chrom_data[idx]
      groups <- groups[idx]
    }
    
    chrom_base <- ggplot() + theme_classic() + xlab("Time(s)") +
      ylab("Intensity") + ggtitle("Chromatogram") +
      theme(axis.line = element_line(colour = "black", size = 0.5,
                                     linetype = "solid"),
            legend.title = element_blank(),
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    chrom_plot <- lapply(seq_along(chrom_data), function(i){
      geom_line(aes(x = chrom_data[[i]][,1], y = chrom_data[[i]][,2], color = groups[i]))
    })
    
    for(i in seq_along(chrom_plot)){
      chrom_base <- chrom_base + chrom_plot[[i]]
    }
    
    return(chrom_base)
  })
  
  # Update rt and bg ranges based on selection
  observe({
    req(input$whichrange, !is.null(input$chr.brush), input$chr.brush)
    minval <- max(round(input$chr.brush$xmin,3),rtRange()[1])
    maxval <- min(round(input$chr.brush$xmax,3),rtRange()[2])
    if (input$whichrange == "Peak") {
      updateSliderInput(inputId = "rt_select",value = c(minval,maxval))
    } else if (input$whichrange == "Background") {
      updateSliderInput(inputId = "bg_select",value = c(minval,maxval))
    }
  })
  
  # Generate final result table
  out <- reactive({
    req(data(),input$ppm,input$rt_select,input$bg_select,mzs(),input$multiplier, !isEmpty(input$background), input$unlabeled)
    
    temp <- get_abundance(data = data(), 
                  ppm = input$ppm, 
                  rt_range = input$rt_select, 
                  bg_range = input$bg_select, 
                  mzs = mzs(), 
                  multiplier = input$multiplier, 
                  background = input$background, 
                  unlabeled = switch(
                    input$unlabeled,
                    "[13]C and [15]N labeled" = NA, 
                    "[13]C labeled" = "N", 
                    "[15]N labeled" = "C"
                  ))
    
    temp <- bind_rows(temp,sapply(temp, rsd))
    rownames(temp)[length(rownames(temp))] <- "RSD"
    
    return(temp)
  }) %>%
    bindCache(data(),input$ppm,input$rt_select,input$bg_select,mzs(),input$multiplier, input$background, input$unlabeled) %>%
    bindEvent(input$run)
  
  output$table <- renderTable({
    req(out())
    out()
  }, rownames = TRUE, digits = 6, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Download button
  output$showDownload <- reactive({
    req(out())
    TRUE
  })
  outputOptions(output, "showDownload", suspendWhenHidden = FALSE)
  
  output$download <- downloadHandler(
    filename = function() {
      "output.csv"
    },
    content = function(file) {
      write.csv(out(), file, row.names = TRUE)
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)