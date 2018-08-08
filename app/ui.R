library(DT)
library(shiny)
library(dplyr)
library(shinycssloaders)
library(rcytoscapejs)
## library("reactomeVisualizer")

cellLines <- readRDS("../data/cellLines.rds")
targets <- readRDS("./targets.rds")
navbarPage("CIE",
           tabPanel("Analysis",
                    ## Using fluid bootstrap layout
                    fluidPage(
                        
                        ## Title the page
                        titlePanel("Causal Inference and Directional Enrichment Methods on Biological Networks"),
                        sidebarPanel(
                            selectInput(inputId = "databaseType",
                                        label = "Database Type",
                                        choices = c("TRED" = "TRED",
                                                    "StringDB" = "string",
                                                    "ChIP Atlas" = "ChIP",
                                                    "TRRUST" = "trrust"),
                                        selected = "ChIP"),
                            conditionalPanel (
                                condition = "input.databaseType == 'ChIP'",
                                selectInput(inputId = "distance",
                                            label = "Distance from TSS",
                                            choices = c("1kb" = "1",
                                                        "5kb" = "5",
                                                        "10kb" = "10"),
                                            selected = "1"),
                                selectInput(inputId = "cutoffType",
                                            label = "Cutoff Type",
                                            choices = c("Minimum" = "min",
                                                        "Maximum" = "max",
                                                        "Average" = "average",
                                                        "Automatic" = "auto"),
                                            selected = "average"),
                                conditionalPanel (condition = "input.cutoffType != 'auto'",
                                                  sliderInput(inputId = "cutoff",
                                                              label = "Binding Score Cutoff",
                                          min = 0, max = 1000,
                                          value = 500,
                                          round=TRUE)
                                          ),
                                
                                selectInput(inputId = "cellLines",
                                            label = "Limit ChIP Atlas results by cell line of experiment",
                                            choices=as.character(c(NA, cellLines$Cell.Line.Name)),
                                            selected= "NA"),
                                selectInput(inputId = "cellLineType",
                                            label = "Limit ChIP Atlas results by cell line tissue origin",
                                            choices=as.character(c(NA, unique(cellLines$Primary.Tissue))),
                                            selected= "NA"),
                                selectInput(inputId = "cellLineDiagnosis",
                                            label = "Limit ChIP Atlas results by the diagnosis of the individual the cell line is from",
                                            choices=as.character(c(NA, unique(cellLines$Tissue.Diagnosis))),
                                            selected= "NA")
                            ),
                            selectInput(inputId = "targetsOfInterest",
                                        label = "Show if the Protien Targets Selected Genes",
                                        choices = c(NA, targets),
                                        selected = NA,
                                        multiple = TRUE),
                            
                            fileInput(inputId = "degFiles",
                                      label = "Upload your differentially expressed gene table, should be in .tsv format",
                                      multiple = FALSE,
                                      accept = c("text/tsv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".tsv")),
                            numericInput(inputId = "p.thresh",
                                         label = "p-Value Threshold",
                                         value = 0.05,
                                         min = 0,
                                         max = 1,
                                         step = 0.01),
                            numericInput(inputId = "fc.thresh",
                                         label = "Select a value, the log of which will be the fold count threshold",
                                         value = "1.5"),
                            selectInput(inputId = "method",
                                        label = "Method of Enrichment Analysis",
                                        choices = c("Ternary", "Quaternary", "Enrichment", "Fisher"),
                                        selected = "Quaternary",
                                        multiple=FALSE),
                            checkboxInput(inputId = "useMart",
                                          label = "Include whether a protein is a transcription factor in enrichment results",
                                          value=TRUE),
                          checkboxInput(inputId = "useBHLH",
                                        label = "Include whether a protein is a BHLH in enrichment results",
                                        value = TRUE),
                          actionButton("run", "Run analysis")
                          ## actionButton("cancel", "Cancel")
            
                        ),
                        mainPanel(
                            tabsetPanel(
                              tabPanel("Graph of Regulators", 
                                       rcytoscapejsOutput("graph") %>%
                                         withSpinner(color="#3498DB", type=8))
                              ## tabPanel("Pathway Overview",
                              ##          tags$div(id = "fireworksHolder"),
                              ##          reactomeVisualizerOutput("pathwaysGraph") %>% 
                              ##           withSpinner(color="#3498DB", type=8))
                            ),
                            ## Working interface for single output
                            column(6, uiOutput("targSlider")),
                            column(6, uiOutput("pathsToShow")),
                            
                            column(12, uiOutput("tableTitle")),
                            column(6, uiOutput("downloadButton")),
                            column(6, uiOutput("pathButton")),
                            
                            tabsetPanel(id="tables",
                                tabPanel("Analysis", DT::dataTableOutput("table") %>%
                                                     withSpinner(color="#3498DB", type = 8)),
                                tabPanel("Pathway Enrichment", DT::dataTableOutput("pathwaysTable") %>%
                                                               withSpinner(color="#3498DB", type = 8)) 
                            )
                            
                            ## Dropping attempt at tab functionality for now
                            ## uiOutput("tabs") %>% withSpinner(color="#0dc5c1")
                        )
                    )),
           tabPanel("Download files",
                    fluidPage(
                        titlePanel('File download'),
                        sidebarLayout(
                            sidebarPanel(
                                selectInput("dataset",
                                            "Choose data to download:", 
                                            choices = list.files("../data/"),
                                            multiple=TRUE),
                                downloadButton('downloadCompr', 'Download')
                            ),
                            mainPanel(
                                tableOutput('tableDownload')
                            )
                        )
                    )
           )
           
)
