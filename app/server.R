
source("../networkAssembly/R/filterChIP.R")
source("../inferenceModels/R/runCIE.R")
source("../inferenceModels/R/runCytoscape.R")
library("stringr")
library("DT")
library("subprocess")
library("parallel")

helperFunctionTable <- function(input, ents, rels, degs) {
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
                         databaseDir = "../data/")
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


server <- function(input, output) {
    rValEntsRels <- reactiveValues()
    rValEntsRels$process <- NULL
    rValEntsRels$msg <- NULL
    rValEnrichment <- reactiveValues()
    rValEnrichment$process <- NULL
    rValEnrichment$msg <- NULL
    
    degs <- eventReactive(input$run, {
      degs <- read.table(input$degFiles$datapath, header=T, sep="\t")
      degs
    })
    observeEvent({input$run}, {
        req(input$databaseType)
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
    observeEvent(rValEntsRels$process, {
      if(is.null(rValEntsRels$result[[1]])) {
        entsRels <- suppressWarnings(isolate(mccollect(rValEntsRels$process$pid)))
        entsRels <- entsRels[[1]]
        names(entsRels) <- c("ents", "rels")
        rValEntsRels$result <- entsRels
      }
    })
    enrichment <- reactive({
      req(rValEntsRels$result)
      entsRels <- rValEntsRels$result
      degs <- degs()
      isolate(helperFunctionTable(input, entsRels$ents, entsRels$rels, degs))
    })
    ## observeEvent({rValEntsRels$process}, {
    ##     degs <- degs()
    ##     if(is.null(rValEntsRels$result[[1]])) {
    ##         entsRels <- suppressWarnings(isolate(mccollect(rValEntsRels$process$pid)))
    ##         entsRels <- entsRels[[1]]
    ##         names(entsRels) <- c("ents", "rels")
    ##         rValEntsRels$result <- entsRels
    ##     }
    ##     else{
    ##         entsRels <- rValEntsRels$result
    ##     }
    ##     rValEnrichment$result <- NULL
    ##     rValEnrichment$process <- mcparallel({ 
    ##         isolate(helperFunctionTable(input, entsRels$ents, entsRels$rels, degs))
    ##     })
    ##     rValEnrichment$msg <- paste("Process",
    ##                                 rValEnrichment$process$pid,
    ##                                 "started")
    ##     print("Started enrichment process")
    ## })
    ## observeEvent({rValEnrichment$process}, {
    ##     rValEnrichment$result <- suppressWarnings(isolate(mccollect(rValEnrichment$process$pid)))
    ##     rValEnrichment$result <- rValEnrichment$result[[1]]
    ##     rValEnrichment$process <- NULL
    ##     rValEnrichment$msg <- paste("Process for enrichment completed!")
    ##     print(head(rValEnrichment$result))
    ##     print("Process for enrichment completed!")
    ## })
    
    ## eventReactive({input$cancel},  {
    ##     if(!is.null(rValEntsRels$process$pid)) {
    ##         tools::pskill(rValEntsRels$process$pid)
    ##         rValEnrichment$msg <-paste("Process",
    ##                                    rValEntsRels$process$pid,
    ##                                    "killed")
    ##         print("EntsRels process killed")
    ##         rValEnrichment$process=NULL
    ##         rValEnrichment$result=NULL
    ##     }
    ##     if(!is.null(rValEnrichment$process$pid)) {
    ##         tools::pskill(rValEnrichment$process$pid)
    ##         rValEnrichment$msg <-paste("Process",
    ##                                    rValEnrichment$process$pid,
    ##                                    "killed")
    ##         print("Enrichment process killed")
    ##         rValEnrichment$process=NULL
    ##         rValEnrichment$result=NULL
    ##     }
    ## })
    output$pathButton <- renderUI({
      req(enrichment())
      req(input$table_rows_selected)
      actionButton("pathEnr", "Run pathway enrichment")
    })
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
    output$pathwaysTable <- DT::renderDataTable({
      req(input$pathsToDisplay)
      prot <- rValEnrichment$result$name[input$table_rows_selected]
      pEnr <- pathwayEnrichment(prot, input$pathsToDisplay)
    })
    output$table <- DT::renderDataTable({
        req(enrichment())
        table <- withProgress(message="Calculating Enrichment\n",
                              detail = "\nThis may take some time",
                              expr = {
                                  enrichment <- rValEnrichment$result
                                  index <- grep("pval|pvalue|p.value|p-value|p-val|p.val",
                                                colnames(enrichment))
                                  index <- unlist(index)
                                  index2 <- grep("adj", colnames(enrichment))
                                  index <- index[index %in% index2]
                                  if(length(index) == 1) {
                                      table <- enrichment %>%
                                          dplyr::select(c(name,
                                                          total_targets = total.reachable,
                                                          significant_targets = significant.reachable,
                                                          index))
                                  }
                                  else {
                                      index2 <- sapply(c("adj"), function(x) {
                                          grep(x, colnames(enrichment))})
                                      index2 <- unlist(index2)
                                      table <- enrichment %>%
                                          dplyr::select(c(name,
                                                          total_targets = total.reachable,
                                                          significant_targets = significant.reachable,
                                                          index2[1], index2[2]))
                                  }
                                  rownames(table) <- enrichment$uid
                                  DT::datatable(table, selection=list(mode='multiple',
                                                                      selected=1:5))
                                  })
        table
    })
    output$tableTitle <- renderUI({
        h3(paste("Enrichment Results for", input$method, "analysis with data base",
                 input$databaseType))
    })
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
    output$targSlider <- renderUI({
        req(enrichment())
        req(maxSlider())
       
        rValEnrichment$result <- enrichment()
        
        sliderInput(inputId = "numTargets",
                    label = "Targets to Display",
                    min = 1, max = maxSlider(),
                    value = 10,
                    round=TRUE)
    })
    output$graph <- renderRcytoscapejs({
        req(input$numTargets)
        req(enrichment())
        req(rValEntsRels$result)
        rValEnrichment$result <- enrichment()
        degs <- degs()
        entsRels <- rValEntsRels$result
        enrichment <- rValEnrichment$result
        req(input$table_rows_selected)
        createCytoGraph(enrichment, entsRels$ents, entsRels$rels, degs,
                        ids=input$table_rows_selected,
                        numTargets=input$numTargets)
    })
    output$downloadButton <- renderUI({
        req(enrichment())
        rValEnrichment$result <- enrichment()
        downloadButton('download', label="Download Full Table")
    })
    output$download <- downloadHandler(
        filename = function() {
            paste(input$method, input$databaseType,
                  input$degFiles, ".tsv", sep="")
        },
        content = function(file) {
            write.table(enrichment(), file, sep="\t",
                        row.names=FALSE, quote=FALSE)
        }
    )
    
    
    
    
    
    
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
