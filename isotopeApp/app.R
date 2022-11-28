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
library(future)
library(future.apply)

# Define UI ----
ui <- navbarPage(
  title = HTML("HRMS"),
  
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  tabPanel(
    "Data Visualizer",
    sidebarLayout(
      sidebarPanel(
        conditionalPanel(
          condition = "output.fileexist",
          uiOutput("chooseFile"),
          uiOutput("addPlot")
        ),
      ),
      
      mainPanel(
        shinyDirButton('datadir', label = "Choose Data Directory", title = "Please select data path", multiple = FALSE),
        br(),
        htmlOutput("filelist")
      )
    )
  ),
  
  tabPanel(
    "Isotope Abundance Calculator",
    
    div(
      actionButton("1","t1"),
      actionButton("2","t2")
    )
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
    list.files(path = parseDirPath(roots, input$datadir), pattern = "^.*\\.mzML$",
                 all.files = TRUE, full.names = TRUE,
                 recursive = FALSE, ignore.case = FALSE,
                 include.dirs = FALSE, no.. = TRUE)  %>%
      mixedsort(decreasing = F)
  })
  # Print filelist
  output$filelist <- renderUI({
    HTML(paste(files(), collapse = '<br/>' ))
  })
  # Read reactive file list
  d <- reactive({
    readMSData(files(), mode = "onDisk")
  })
  # Make reactive file numbers
  filenums <- reactive({
    d() %>% fromFile() %>% unique()
  })
  # Make reactive status of file existence
  output$fileexist <- reactive({
    sum(filenums()) > 0
  })
  outputOptions(output, "fileexist", suspendWhenHidden = FALSE)
  # Generate conditional radio buttons
  output$chooseFile <- renderUI({
    radioButtons("filechosen","Choose a file to plot",choices = filenums(),selected = character(0), inline = TRUE)
  })
  # Filter centroided files depending on input
  # Generate conditional button
  output$addPlot <- renderUI({
    actionButton("addButton","Add a plot")
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)