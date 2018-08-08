
source("../networkAssembly/R/filterChIP.R")
source("../inferenceModels/R/runCIE.R")
source("../inferenceModels/R/runCytoscape.R")
library("stringr")
library("DT")
library("parallel")
library("utils")
library("shinyjs")
## library("reactomeVisualizer")

helperFunctionTable <- function(input, ents, rels, degs, progress) {
    enrichment <- runCIE(NULL, NULL, DEGs = degs,
                         p.thresh = input$p.thresh,
                         fc.thresh = log(input$fc.thresh),
                         methods = input$method,
                         useFile=FALSE,
                         ents = ents,
                         rels = rels,
                         useMart = input$useMart,
                         martFN = "../data/mart_human_TFs.csv",
                         useBHLH = input$useBHLH,
                         BHLHFN = "../data/BHLH_TFs.txt",
                         hypTabs="1",
                         verbose=F,
                         databaseDir = "../data/",
                         expectProgressObject=TRUE,
                         progress=progress)
    ## Code which will display the enrichment given that only one output has been
    ## provided, which is the case due to minor changes in the interface.
    enrichment

    ## Code with unknown bug that would allow for multiple outputs as
    ## I intended when writing the function.
    ## if((class(enrichment) == "data.frame") &&
    ##     (class(degs) == "data.frame")) {
    ##   resultTitles <- paste(input$method, input$degFiles$name, sep = ".")
    ##   tables <- enrichment
    ## }
    ## else if((class(enrichment) == "list") &&
    ##    (class(enrichment[[1]]) == "data.frame") &&
    ##    (class(degs == "list"))) {
    ##   resultTitles <- paste(input$method, names(enrichment), sep = ".")
    ##   tables <- lapply(enrichment, function(x) {x} )
    ## }
    ## else if((class(enrichment) == "list") &&
    ##         (class(enrichment[[1]]) == "data.frame") &&
    ##         (class(degs == "data.frame"))) {
    ##   resultTitles <- paste(names(enrichment), input$degFiles$name, sep = ".")
    ##   tables <- lapply(enrichment, function(x) {x} )
    ## }
    ## else if((class(enrichment) == "list") &&
    ##    (class(enrichment[[1]]) == "list")) {
    ##     resultTitles <- lapply(1:length(enrichment), function(x) {
    ##         lapply(1:length(enrichment[[1]]), function(y) {
    ##             paste(names(enrichment)[x], names(enrichment[[x]])[y],
    ##                   sep = "\\.") } ) } )
    ##     tables <- lapply(enrichment, function(x) {
    ##         lapply(enrichment[[x]], function(y) {
    ##             y } ) } )
    ## }
    ## names(tables) <- resultTitles
    ## tabTitles <- str_replace_all(as.character(resultTitles), ".", " ")
    ## print("Return")
    ## return(list(tables = tables, tabTitles = tabTitles))
}


server <- function(input, output, session) {
    ## Code for analysis
    ## Code to have the ents and rels tables be calculated across reactive calls
    ## with mcparallel
    rValEntsRels <- reactiveValues()
    rValEntsRels$process <- NULL
    rValEntsRels$msg <- NULL

    startEnr <- NULL

    ## Partially used code to store enrichment in a similar manner
    ## rValEnrichment <- reactiveValues()
    ## rValEnrichment$process <- NULL
    ## rValEnrichment$msg <- NULL

    observeEvent({input$run}, {
        shinyjs::disable("run")
        shinyjs::disable("method")
        shinyjs::disable("fc.thresh")
        shinyjs::disable("p.thresh")
        shinyjs::disable("degFiles")
        shinyjs::disable("targetsOfInterest")
        shinyjs::disable("cellLineDiagnosis")
        shinyjs::disable("cellLineType")
        shinyjs::disable("cellLines")
        shinyjs::disable("cutoffType")
        shinyjs::disable("distance")
        shinyjs::disable("databaseType")
        shinyjs::disable("useBHLH")
        shinyjs::disable("useMart")
    })
    observeEvent({enrichment()}, {
        shinyjs::enable("run")
        shinyjs::enable("method")
        shinyjs::enable("fc.thresh")
        shinyjs::enable("p.thresh")
        shinyjs::enable("degFiles")
        shinyjs::enable("targetsOfInterest")
        shinyjs::enable("cellLineDiagnosis")
        shinyjs::enable("cellLineType")
        shinyjs::enable("cellLines")
        shinyjs::enable("cutoffType")
        shinyjs::enable("distance")
        shinyjs::enable("databaseType")
        shinyjs::enable("useBHLH")
        shinyjs::enable("useMart")
   
    })

    ## Get degs
    degs <- eventReactive({input$run}, {
        degs <- read.table(input$degFiles$datapath, header=T, sep="\t")
        degs
    })

    progress <- reactiveValues()
    progress$result <- NULL

    ## Get ents rels
    observeEvent({input$run}, {
        req(input$databaseType)
        progress$result <- shiny::Progress$new()
        progress$result$set(message="Preparing network for enrichment analysis",
                             value=0)
        if(input$cutoffType == "auto") {
            cutoff=NA
        }
        else{
            cutoff=input$cutoff
        }
        if(input$cellLines == "NA") {
            cellLines = NA
        }
        else {
            cellLines <- input$cellLines
        }
        
        if(input$cellLineType == "NA") {
            cellLineType = NA
        }
        else {
            cellLineType <- input$cellLineType
        }
        
        if(input$cellLineDiagnosis== "NA") {
            cellLineDiagnosis = NA
        }
        else {
            cellLineDiagnosis <- input$cellLineDiagnosis
        }
        if(input$databaseType == "ChIP") {
            if(!is.null(rValEntsRels$process)) {
                return()
            }
            else {
                rValEntsRels$result <- NULL
                rValEntsRels$process <- mcparallel({
                    isolate(filterChIPAtlas(distance = as.numeric(input$distance),
                                            cutoff = cutoff,
                                            cutoffType = input$cutoffType,
                                            cellLines = cellLines,
                                            cellLineType = cellLineType,
                                            cellLineDiagnosis = cellLineDiagnosis,
                                            writeToFile=FALSE,
                                            databaseDir = "../data/"))
                })
                rValEntsRels$msg <- paste("Process ",rValEntsRels$process$pid,
                                          " started")
                print("Started ents rels")
            }
            
        }
        else {
            entsFile <- paste("../data/", input$databaseType, ".ents", sep = "")
            relsFile <- paste("../data/", input$databaseType, ".rels", sep = "")
            ents <- read.table(entsFile, sep = "\t", header = T,
                               stringsAsFactors = F)
            rels <- read.table(relsFile, sep = "\t", header = T,
                               stringsAsFactors = F)
            rValEntsRels$result = list(ents=ents, rels=rels)
            rValEntsRels$process = 0
            rValEntsRels$msg = "Loaded from file"
        }
    })

    # Collect ents rels (if processed, not read)
    observeEvent(rValEntsRels$process, {
      if(is.null(rValEntsRels$result[[1]])) {
          entsRels <- suppressWarnings(isolate(mccollect(rValEntsRels$process$pid)))
          entsRels <- entsRels[[1]]
          names(entsRels) <- c("ents", "rels")
          rValEntsRels$result <- entsRels
      }
    })

    ##Get enrichment
    enrichment <- reactive({
        req(rValEntsRels$result)
        value <- progress$result$getValue()
        value <- value + ((progress$result$getMax() - value) / 10)
        progress$result$set(message="Preparing DGE files to run Enrichment",
                             value=value)
        entsRels <- rValEntsRels$result
        degs <- degs()
        on.exit(progress$result$close())
        isolate(helperFunctionTable(input, entsRels$ents, entsRels$rels, degs, progress$result))
     })

    ## Render pathway enrichment button
    output$pathButton <- renderUI({
      req(enrichment())
      req(input$table_rows_selected)
      actionButton("pathEnr", "Run pathway enrichment")
    })
    ## Add slider for pathway enrichment
    output$pathsToShow <- renderUI({
      req(input$pathEnr)
      req(enrichment())
      req(input$table_rows_selected)
      sliderInput(inputId = "pathsToDisplay",
                  label = "Choose How Many Paths to Display",
                  min = 0, max = 100,
                  value = 10,
                  round=TRUE)
    })
    ## Render pathway enrichment
    pEnr <- reactive({
        req(input$pathsToDisplay)
        print("Calculating pathway enrichment")
        prot <- enrichment()$name[input$table_rows_selected]
        pEnr <- pathwayEnrichment(prot, input$pathsToDisplay)      
    })
    
    observeEvent({is.null(pEnr())}, {
        print("Switching tabs")
        updateTabsetPanel(session, "tables",
                          selected = "Pathway Enrichment")
    })

    output$pathwaysTable <- DT::renderDataTable({
        req(pEnr())
        print("Showing table")
        DT::datatable(pEnr(), selection=list(mode='single', selected=1))
    })

    ## Render enrichment
    output$table <- DT::renderDataTable({
        req(enrichment())
        enrichment <- enrichment()
        index <- grep("pval|pvalue|p.value|p-value|p-val|p.val",
                      colnames(enrichment))
        index <- unlist(index)
        index2 <- grep("adj", colnames(enrichment))
                                  index <- index[!(index %in% index2)]
        if(length(index) == 1) {
            table <- enrichment %>%
                dplyr::select(c(name,
                                total_targets = total.reachable,
                                significant_targets = significant.reachable,
                                index))
        }
        else {
            table <- enrichment %>%
                dplyr::select(c(name,
                                total_targets = total.reachable,
                                significant_targets = significant.reachable,
                                index[1], index[2]))
        }
        rownames(table) <- enrichment$uid
        DT::datatable(table, selection=list(mode='multiple',
                                            selected=1:5))  
    })

    ## Title analysis
    output$tableTitle <- renderUI({
        h3(paste("Enrichment Results for", input$method, "analysis with data base",
                 input$databaseType))
    })

    ## Value for slider for rcytoscape js
    maxSlider <- reactive({
        req(enrichment())
        req(input$table_rows_selected)
        ids <- input$table_rows_selected
        numTarg <- enrichment()$significant.reachable[ids]
        if(max(numTarg) > 100) {
            100
        }
        else {
            max(numTarg)
        }
    })
    ## Slider for rcytoscape js
    output$targSlider <- renderUI({
        req(enrichment())
        req(maxSlider())
       
        sliderInput(inputId = "numTargets",
                    label = "Targets to Display",
                    min = 1, max = maxSlider(),
                    value = 10,
                    round=TRUE)
    })
    ##render r cytoscape js
    output$graph <- renderRcytoscapejs({
        req(input$numTargets)
        req(enrichment())
        req(rValEntsRels$result)
        degs <- degs()
        entsRels <- rValEntsRels$result
        enrichment <- enrichment()
        req(input$table_rows_selected)
        createCytoGraph(enrichment, entsRels$ents, entsRels$rels, degs,
                        ids=input$table_rows_selected,
                        numTargets=input$numTargets)
    })
    
    ## output$pathwaysGraph <- renderReactomeVisualizer({
    ##     print("Started output visualization")
    ##     req(output$pathwaysTable)
    ##     print(paste("Visualizing row: ", input$pathwaysTable_rows_selected,
    ##                 sep = ""))
    ##     pathId <- pEnr()$id[input$pathwaysTable_rows_selected]
    ##     reactomeVisualizer(pathId, "fireworksHolder")})

    ## Download full table
    output$downloadButton <- renderUI({
        req(enrichment())
        downloadButton('download', label="Download Full Table")
    })
    ## Server of files
    output$download <- downloadHandler(
        filename = function() {
            paste(input$method, input$databaseType,
                  input$degFiles, ".tsv", sep="")
        },
        content = function(filename) {
            write.table(enrichment(), filename, sep="\t",
                        row.names=FALSE, quote=FALSE)
        }
    )
    ## Download for running script locally
    ## Server of files
    output$downloadCompr <- downloadHandler(
        filename = function() {
            names <- names <- sapply(input$dataset, function(x) {strsplit(x, split="\\.")[[1]][1]})
            name <- paste(names, collapse="-")
            paste(name, ".zip", sep="")
        },
        content = function(filename) {
            zip(filename, paste("../data/", input$dataset, sep=""))
        },
        contentType = "application/zip"
    )
    
    output$description <- renderUI({
        tags$h3("Running CIE on Your Own Machine")
        tags$p(paste("Running CIE on your own ", sep=""))
    })
    
    
    
    
    ## Code which should allow a varying number of tabs to display the results
    ## of the pipeline but does not due to an unknown bug
    ## output$tabs <- renderUI({
    ##     req(input$degFiles)
    ##     tablesTabs <- helperFunctionTable(input)
    ##     if(length(tablesTabs$tabTitles) > 1) {
    ##       tables <- do.call(renderTable, tablesTabs$tables)
    ##       tablesOut <- do.call(tableOutput, tables)
    ##       tabs <- lapply(1:length(tabTitles), function(x) {
    ##         tabPanel(tablesTabs$tabTitles[x], tablesOut)
    ##       })
    ##       do.call(tabsetPanel, tabs)
    ##     }
    ##     else {
    ##       table <- renderTable(tablesTabs$tables, striped = TRUE)
    ##       tableOut <- tableOutput(table)
    ##       tabPanel(tablesTabs$tabTitles, table)
    ##     }
    ## })
}
