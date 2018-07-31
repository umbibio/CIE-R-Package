
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
    
    degs <- reactive({
      req(input$degFiles)
      if(length(input$degFiles$datapath) > 1) {
        degs <- lapply(input$degFiles$datapath, function(x) {
          read.table(x, header=T, sep="\t") } )
      }
      else {
        degs <- read.table(input$degFiles$datapath, header=T, sep="\t")
      }
    })
    observeEvent(input$run, {
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
        degs <- degs()
        if(is.null(rValEntsRels$result[[1]])) {
            entsRels <- suppressWarnings(isolate(mccollect(rValEntsRels$process$pid)))
            entsRels <- entsRels[[1]]
            names(entsRels) <- c("ents", "rels")
            rValEntsRels$result <- entsRels
        }
        else{
            entsRels <- rValEntsRels$result
        }
        rValEnrichment$result <- NULL
        rValEnrichment$process <- mcparallel({ 
            isolate(helperFunctionTable(input, entsRels$ents, entsRels$rels, degs))
        })
        rValEnrichment$msg <- paste("Process",
                                    rValEnrichment$process$pid,
                                    "started")
        print("Started enrichment process")
    })
    observeEvent(rValEnrichment$process, {
        rValEnrichment$result <- suppressWarnings(isolate(mccollect(rValEnrichment$process$pid)))
        rValEnrichment$result <- rValEnrichment$result[[1]]
        rValEnrichment$process <- NULL
        rValEnrichment$msg <- paste("Process for enrichment completed!")
        print("Process for enrichment completed!")
    })
    
    ## observeEvent(input$cancel,  {
    ##   print("test")
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
    output$pathways <- DT::renderDataTable({
      req(rValEnrichment$result)
      req(input$table_rows_selected)
      prot <- rValEnrichment$result$name[input$table_rows_selected]
      print(prot)
      pEnr <- pathwayEnrichment(prot)
      
      DT::datatable(pEnr)
    })
    output$table <- DT::renderDataTable({
        req(rValEnrichment$result)
        table <- withProgress(message="Calculating Enrichment\n",
                              detail = "\nThis may take some time",
                              expr = {
                                  enrichment <- rValEnrichment$result
                                  index <- grep("pval|pvalue|p.value|p-value|p-val|p.val",
                                                colnames(enrichment))
                                  index <- unlist(index)
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
        req(rValEnrichment$result)
        req(rValEntsRels$result)
        req(input$table_rows_selected)
        degs <- degs()
        entsRels <- rValEntsRels$result
        ents <- entsRels$ents
        rels <- entsRels$rels
        enrichment <- rValEnrichment$result
        ids <- input$table_rows_selected
        ## Code written by Dr. Kourosh Zarringhalam
        pval.ind = grep('qval|q.val|q-val|q-val|P-value|P.value|pvalue|pval|Pval',
                        colnames(degs), ignore.case = T)
        fc.ind = grep('fc|FC|fold|FoldChange', colnames(degs), ignore.case = T)
        id.ind = grep('id|entr|Entrez', colnames(degs), ignore.case = T)
        
        if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
            stop('Please make sure the expression files column names are labled as entrez, fc, pvalue')
        }
            
        colnames(degs)[pval.ind] <- 'pvalue'
        colnames(degs)[fc.ind] <- 'foldchange'
        colnames(degs)[id.ind] <- 'id'
        
        sigDEG <- degs %>% filter(abs(foldchange) >= input$fc.thresh & 
                                    pvalue <= input$p.thresh) %>%
            transmute(id = id, val = ifelse(foldchange > 0, 1, -1), pval = pvalue) %>%
                distinct(id, .keep_all = T)
        
        ## End Dr. Zarringhalam cotde
        sigProt <- enrichment$uid[ids]
        sigEnts <- ents[ents$id %in% sigDEG$id, ]
        sigRels <- rels[rels$srcuid %in% sigProt, ]
        numTarg <- sapply(sigProt, function(x) {
            tarRels <- sigRels %>% dplyr::filter(srcuid == x)
            targs <- sigRels$trguid
            targs <- targs[targs %in% sigEnts$uid]
            length(targs)
        })
        if(max(numTarg) > 100) {
          100
        }
        else {
          max(numTarg)
        }
    })
    output$targSlider <- renderUI({
      req(rValEnrichment$result)
      req(maxSlider())
      sliderInput(inputId = "numTargets",
                  label = "Targets to Display",
                  min = 1, max = maxSlider(),
                  value = 10,
                  round=TRUE)
    })
    output$graph <- renderRcytoscapejs({
        req(input$numTargets)
        req(rValEnrichment$result)
        req(rValEntsRels$result)
        degs <- degs()
        entsRels <- rValEntsRels$result
        enrichment <- rValEnrichment$result
        req(input$table_rows_selected)
        createCytoGraph(enrichment, entsRels$ents, entsRels$rels, degs,
                        ids=input$table_rows_selected,
                        numTargets=input$numTargets)
        
    })
    output$downloadButton <- renderUI({
        req(rValEnrichment$result)
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
