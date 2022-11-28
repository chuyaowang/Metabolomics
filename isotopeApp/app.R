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
library(praise)
library(future.apply)
source("./helpers.R")

# Define UI ----
ui <- navbarPage(
  title = HTML("HRMS"),
  
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  tabPanel(
    "Isotope Abundance Calculator",
    sidebarLayout(
      sidebarPanel(
        shinyDirButton('datadir', label = "Choose Data Directory", title = "Please select data path", multiple = FALSE, width = "60%"),
        conditionalPanel(
          condition = "output.showPanel",
          br(),
          checkboxInput("centroid", label = "Use centroided data"),
          helpText("Enter mass error in ppm"),
          numericInput("ppm",label = NULL, value = 5, min = 0, max = 50, step = 1, width = "50%"),
          helpText("Retention Time Range (sec);"),
          fluidRow(
            column(6, textInput("rt_lower","Lower limit",width= "100%")),
            column(6, textInput("rt_higher","Higher limit", width = "100%"))
          ),
          helpText("Enter compound formula"),
          fluidRow(
            column(3,textInput("C","C")),
            column(3,textInput("H","H")),
            column(3,textInput("N","N")),
            column(3,textInput("O","O"))
          ),
          helpText("Choose compound labeling status"),
          selectInput("unlabeled", label = NULL, choices = c("[13]C and [15]N labeled", "[13]C labeled", "[15]N labeled")),
          helpText("Choose data polarity"),
          selectInput("polarity", label = NULL, choices = c("+","-"), width = "50%"),
          helpText("Minimum percent of max intensity for a spectrum to be used in calculation"),
          numericInput("multiplier", label = NULL, value = 1, min = 0.1, max = 1, step = 0.1, width = "50%"),
          checkboxInput("background", label = "Background subtraction (beta)")
        ),
        conditionalPanel(
          condition = "output.showPlot",
          actionButton("plotData", "Plot")
        ),
        conditionalPanel(
          condition = "input.background",
          helpText("Background Range (sec);"),
          fluidRow(
            column(6, textInput("bg_lower","Lower limit",width= "100%")),
            column(6, textInput("bg_higher","Higher limit", width = "100%"))
          )
        )

        
      ),
      
      mainPanel(
        "Files to be processed:",
        br(),
        htmlOutput("filelist"),
        hr(),
        conditionalPanel(
          condition = "output.showPanel",
          fluidRow(
            column(3, selectInput("mz", label = "Choose mz value", choices = c(1,2,3,4,5,6))),
            column(3, selectInput("whichfile", label = "Choose file to plot", choices = c(1,2,3))),
            column(6, selectInput("whichrange", label = "Range", choices = c("Retention time", "Background")))
          ),
          plotOutput("chr")
        )
      )
    )
  ),
  
  tabPanel(
    HTML("Help &#x1F6C8;")
  )
)

# Define server logic ----
server <- function(input, output) {
  # Choose files ----
  # Root folder
  roots <- c(data = '../data')
  # Connect to root folder
  shinyDirChoose(input, "datadir", roots = roots, allowDirCreate = TRUE)
  # Make reactive file list
  files <- reactive({
    fls <- list.files(path = parseDirPath(roots, input$datadir), pattern = "^.*\\.mzML$",
                 all.files = TRUE, full.names = TRUE,
                 recursive = FALSE, ignore.case = FALSE,
                 include.dirs = FALSE, no.. = TRUE)  %>%
      mixedsort(decreasing = F)
    if (input$centroid == TRUE) {
      fls_cent <- fls[str_detect(fls,"cent")]
      if (length(fls_cent)==0) {
        fls
      } else {
        fls_cent
      }
    } else {
      fls <- fls[!str_detect(fls, "cent")]
      fls
    }
  })
  # Print filelist
  output$filelist <- renderUI({
    if (input$centroid == TRUE) {
      check <- files() %>% str_detect(.,"cent") %>% sum
      if (check>0) {
        HTML(paste(files(), collapse = '<br/>' ))
      } else if (check == 0) {
        HTML(paste(c(files(),"<br/>Centroids will be computed on the fly.")), collapse = '<br/>')
      }
    } else if (input$centroid == FALSE) {
      HTML(paste(files(), collapse = '<br/>' ))
    }
  })
  # Condition for showing other controls
  output$showPanel <- reactive({
    length(files()) > 0
  })
  outputOptions(output, "showPanel", suspendWhenHidden = FALSE)
  # Condition for showing plot button
  output$showPlot <- reactive({
    (input$C != "") & (input$H != "") & (input$N != "") & (input$O != "")
  })
  outputOptions(output, "showPlot", suspendWhenHidden = FALSE)
}

# Run the app ----
shinyApp(ui = ui, server = server)