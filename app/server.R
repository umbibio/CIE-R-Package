
source("../networkAssembly/R/filterChIP.R")
source("../inferenceModels/R/runCIE.R")
source("../inferenceModels/R/runCytoscape.R")
library("stringr")
library("DT")

helperFunctionTable <- function(input, ents, rels) {
    if(length(input$degFiles$datapath) > 1) {
        degs <- lapply(input$degFiles$datapath, function(x) {
            read.table(x, header=T, sep="\t") } )
    }
    else {
        degs <- read.table(input$degFiles$datapath, header=T, sep="\t")
    }
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
                         verbose=F)
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
    entsRels <- reactive({
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
            ChIP <- filterChIPAtlas(distance = as.numeric(input$distance),
                                    cutoff = cutoff,
                                    cutoffType = input$cutoffType,
                                    cellLines = cellLines,
                                    cellLineType = cellLineType,
                                    cellLineDiagnosis = cellLineDiagnosis,
                                    writeToFile=FALSE)
            ents <- ChIP$filteredChIP.ents
            rels <- ChIP$filteredChIP.rels
        }
        else {
            entsFile <- paste("../data/", input$databaseType, ".ents", sep = "")
            relsFile <- paste("../data/", input$databaseType, ".rels", sep = "")
            ents <- read.table(entsFile, sep = "\t", header = T,
                               stringsAsFactors = F)
            rels <- read.table(relsFile, sep = "\t", header = T,
                           stringsAsFactors = F)
        }
        list(ents=ents, rels=rels)
    })
    ## output$targetSelector <- renderUI({
    ##  req(entsRels())
    ##  ents <- entsRels()$ents
    ##  targets <- ents[ents$type == "mRNA", colnames(ents) %in% "name"]
    ##  selectInput(inputId = "targetsOfInterest",
    ##              label = "Show if the Protien Targets Selected Genes",
    ##              choices = c(NA, targets),
    ##              selected = NA,
    ##              multiple = TRUE)
    ##})
    enrichment <- eventReactive(input$run, {
        req(input$degFiles)
        helperFunctionTable(input, entsRels()$ents, entsRels()$rels)
    })
    output$table <- DT::renderDataTable({
        table <- withProgress(message="Calculating Enrichment\n",
                     detail = "\nThis may take some time",
                     expr = {
                         req(enrichment())
                         enrichment <- enrichment()
                         index <- grep("pval|pvalue|p.value|p-value|p-val|p.val",
                                       colnames(enrichment))
                         index <- unlist(index)
                         if(length(index) == 1) {
                             table <- enrichment %>%
                                 dplyr::select(c(name, total_targets = total.reachable,
                                                 significant_targets = significant.reachable,
                                                 index))
                         }
                         else {
                             index2 <- sapply(c("adj"), function(x) {grep(x, colnames(enrichment))})
                             index2 <- unlist(index2)
                             table <- enrichment %>%
                                 dplyr::select(c(name, total_targets = total.reachable,
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
    output$graph <- renderRcytoscapejs({
        req(enrichment())
        req(input$table_rows_selected)
        if(length(input$degFiles$datapath) > 1) {
            degs <- lapply(input$degFiles$datapath, function(x) {
                read.table(x, header=T, sep="\t") } )
        }
        else {
            degs <- read.table(input$degFiles$datapath, header=T, sep="\t")
        }
        createCytoGraph(enrichment(), entsRels()$ents, entsRels()$rels, degs,
                        ids=input$table_rows_selected,
                        numTargets=input$numTargets)
        
    })
    output$downloadButton <- renderUI({
        req(enrichment())
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
