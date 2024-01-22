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

  title = div(img(src="logo.png",height=46,width=46),align="left"),
  
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  windowTitle = "IsoCalc v1.3.3",
  
  tags$head(HTML("<link rel='icon' type='png' href='logo.png'>")),
  
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
          br(),
          helpText("Choose labeling mode"),
          radioButtons("mode",label = NULL, choices = c("C-N Labeled","C-H Labeled","N-H Labeled"),inline = TRUE),
          helpText("Enter compound formula"),
          fluidRow(
            column(6,numericInput("C","C",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("H","H",value = NULL,min=0,max=60,step=1))
          ),
          fluidRow(
            column(6,numericInput("N","N",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("O","O",value = NULL,min=0,max=60,step=1))
          ),
          fluidRow(
            column(6,numericInput("S","S",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("P","P",value = NULL,min=0,max=60,step=1))
          ),
          fluidRow(
            column(6,numericInput("f","F",value = NULL,min=0,max=60,step=1)),
            column(6,numericInput("I","I",value = NULL,min=0,max=60,step=1))
          ),
          uiOutput("formulaInput"),
          verbatimTextOutput("test"),
          helpText("Enter mass error in ppm"),
          numericInput("ppm",label = NULL, value = 5, min = 0, max = 50, step = 1, width = "50%"),
          helpText("Choose data polarity"),
          selectInput("polarity", label = NULL, choices = c("+","-"), width = "50%"),
          helpText("Fill in the compound formula to proceed"),
        ),
        conditionalPanel(
          condition = "output.showPlot",
          helpText("Peak Range (sec);"),
          uiOutput("rt_control"),
          helpText("Use intensites greater than what percent of max intensity to compute average?"),
          numericInput("multiplier", label = NULL, value = 1, min = 0.1, max = 1, step = 0.1, width = "50%"),
          helpText("Choose main peak mz"),
          uiOutput("main_peak"),
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
        "Note: for the sample glycine data, the formula is C2H5NO2 containing 2 [13]Cs and 1 [15]N, fill in 0s for the rest of elements. The data is measured in positive mode. The 79.041 peak should be selected as the main peak (peak with the highest response)",
        hr(),
        "Files to be processed:",
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
  ),
  
  # Changelog
  tabPanel(
    HTML("Change Log"),
    
    includeHTML("./www/changelog.html")
  )
)

# Define server logic ----
server <- function(input, output) {
  # Root folder
  roots <- c(data = 'sampledata/') # demo mode
  # roots <- c(data = '../data/') # change to your data folder
  # Connect to root folder
  shinyDirChoose(input, "datadir", roots = roots, allowDirCreate = F)
  
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
    req(input$C,input$H,input$N,input$O,input$S,input$P,input$f,input$I,input$label1,input$label2,data())
    TRUE
  })
  outputOptions(output, "showPlot", suspendWhenHidden = FALSE)
  
  # Read reactive data item
  data <- reactive({
    req(files(),input$centroid)
    
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
  
  # Generate reactive compound formula input
  output$formulaInput <- renderUI({
    req(input$mode)
    if (input$mode == "C-N Labeled") {
      req(input$C,input$N)
      elem <- tagList(
        helpText("Enter the number of labeled atoms"),
        fluidRow(
          column(6,numericInput("label1","[13]C",value = input$C,min=0,max=input$C,step=1)),
          column(6,numericInput("label2","[15]N",value = input$N,min=0,max=input$N,step=1))
        )
      )
    } else if (input$mode == "C-H Labeled") {
      req(input$C,input$H)
      elem <- tagList(
        helpText("Enter the number of labeled atoms"),
        fluidRow(
          column(6,numericInput("label1","[13]C",value = input$C,min=0,max=input$C,step=1)),
          column(6,numericInput("label2","D",value = input$H,min=0,max=input$H,step=1))
        )
      )
    } else if (input$mode == "N-H Labeled") {
      req(input$N,input$H)
      elem <- tagList(
        helpText("Enter the number of labeled atoms"),
        fluidRow(
          column(6,numericInput("label1","[15]N",value = input$N,min=0,max=input$N,step=1)),
          column(6,numericInput("label2","D",value = input$H,min=0,max=input$H,step=1))
        )
      )
    }
    elem
  })
  
  # Report unlabeled status
  unlabeled <- reactive({
    req(input$mode,input$label1,input$label2)
    if (input$mode == "C-N Labeled") {
      if (input$label1 == 0) {
        temp <- "[13]C"
      } else if (input$label2 == 0) {
        temp <- "[15]N"
      } else {
        temp <- NA
      }
    } else if (input$mode == "C-H Labeled") {
      if (input$label1 == 0) {
        temp <- "[13]C"
      } else if (input$label2 == 0) {
        temp <- "D"
      } else {
        temp <- NA
      }
    } else if (input$mode == "N-H Labeled") {
      if (input$label1 == 0) {
        temp <- "[15]N"
      } else if (input$label2 == 0) {
        temp <- "D"
      } else {
        temp <- NA
      }
    }
    return(temp)
  })
  
  # output$test <- renderText({
  #   input$label2
  # })

  # Generate reactive plot options
  pol <- reactive({
    switch(
      input$polarity,
      "+" = 1,
      "-" = -1
    )
  })
  
  mzs <- reactive({
    req(input$C, input$N, input$H, input$O, input$S, input$P, input$f, input$I, input$label1, input$label2, pol(), data(),input$mode)
    
    if (input$mode == "C-N Labeled") {
    mzs <- get_mz_cn(C=input$C,N=input$N,H=input$H,O=input$O,S=input$S,P=input$P,f=input$f, i=input$I, c13=input$label1,n15=input$label2,pol=pol())
    } else if (input$mode == "C-H Labeled") {
      mzs <- get_mz_ch(C=input$C,N=input$N,H=input$H,O=input$O,S=input$S,P=input$P,f=input$f, i=input$I,c13=input$label1,d=input$label2,pol=pol())
    } else if (input$mode == "N-H Labeled") {
      mzs <- get_mz_nh(C=input$C,N=input$N,H=input$H,O=input$O,S=input$S,P=input$P,f=input$f, i=input$I,n15=input$label1,d=input$label2,pol=pol())
    }

    return(mzs)
    }) %>%
    bindCache(input$C, input$N, input$H, input$O, input$S, input$P, input$f, input$I, input$label1, input$label2, pol())

  output$mztable <- renderTable({
    req(mzs())
    mzs()
  }, rownames = TRUE, digits = 4, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$plotOptions <- renderUI({
    req(mzs(),dataNames(),input$mainPeak)
    mzs_vec <- as.vector(as.matrix(mzs()))
    fluidRow(
      column(3,selectInput("whichmz", label = "Choose mz value", choices = round(mzs_vec,4), selected = input$mainPeak)),
      column(3, selectInput("whichfile", label = "Choose file to plot", choices = c(dataNames(),"All"),selected="All")),
      column(6, selectInput("whichrange", label = "Select for", choices = c("Peak", "Background")))
    )
  }) %>% 
    bindEvent(input$plotData)
  
  # Generate reactive rt and background controls
  rtRange <- reactive({
    req(data())
    data() %>% filterFile(1) %>% fData %>% select(retentionTime) %>% c(min(.),max(.)) %>% `[`(-1) %>% unlist %>% round(.,3)
  }) %>%
    bindCache(data())
  
  output$rt_control <- renderUI({
    sliderInput("rt_select",label = NULL, min = rtRange()[1], max = rtRange()[2], value = rtRange(), step = 0.001, round = FALSE, ticks = FALSE, dragRange = TRUE, animate = FALSE)
  })
  
  output$main_peak <- renderUI({
    # Choose the main peak
    req(mzs)
    mzs_vec <- as.vector(as.matrix(mzs()))
    selectInput("mainPeak",label = NULL, choices=round(mzs_vec,digits=4),multiple=F)
  }) %>% 
    bindCache(mzs())
  
  output$bg_control <- renderUI({
    sliderInput("bg_select",label = NULL, min = rtRange()[1], max = rtRange()[2], value = rtRange(), step = 0.001, round = FALSE, ticks = FALSE, dragRange = TRUE, animate = FALSE)
  })
  
  # Generate plot
  thematic::thematic_shiny()
  chr <- reactive({
    req(data(),input$ppm, input$whichmz)
    data() %>% filterMz(mz = get_ppm_range(x = as.numeric(input$whichmz), ppm = input$ppm)) %>% chromatogram(aggregationFun = "max", missing = 0)
  }) %>% 
    bindCache(data(),input$ppm,input$whichmz)
  
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
    req(data(),input$ppm,input$rt_select,input$bg_select,mzs(),input$multiplier, !isEmpty(input$background), !isEmpty(unlabeled()))
    
    temp <- get_abundance(data = data(), 
                  ppm = input$ppm, 
                  rt_range = input$rt_select, 
                  bg_range = input$bg_select, 
                  mzs = mzs(), 
                  multiplier = input$multiplier, 
                  background = input$background,
                  mainPeak = input$mainPeak,
                  unlabeled = unlabeled())
    
    temp <- bind_rows(temp,sapply(temp, rsd))
    rownames(temp)[length(rownames(temp))] <- "RSD"
    
    return(temp)
  }) %>%
    bindCache(data(),input$ppm,input$rt_select,input$bg_select,mzs(),input$multiplier, input$background, input$mainPeak, unlabeled()) %>%
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