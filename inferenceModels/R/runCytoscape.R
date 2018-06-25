require(rcytoscapejs)
require(dplyr)
require(RColorBrewer)

#' Create Graphs of significant interactions
#'
#' @description Creates graphs of protein-gene interactions based on enrichment analyses
#' by this pipeline
#'
#' @usage createCytoGraph(enrichment, ents, rels, DEGs)
#'
#' @param enrichment The enrichment tables from the analysis pipeline.  Can be a list
#' of lists (multiple methods and conditions) or a list (multiple conditions, one
#' method), or a single data frame
#'
#' @param ents The ents table used in the enrichment, must be a single data frame
#'
#' @param rels The rels table used in the enrichment, must be a single data frame
#'
#' @param DEGs The differentially expressed gene tables, can be a single data frame
#' or list of data frames matching that which the pipeline was run on.
#'
#' @return Opens a browser window with the graph, green = proteins, purple = mRNA
#'
#' @export
#'
#' @examples
#'
#' ChIP1ap <- filterChIPAtlas(1, NA, "auto", NA, "prostate", NA, NA, FALSE)
#' files <- list.files("./", ".txt")
#' degs <- lapply(files, function(x) {
#'    read.table(x, header=T, sep="\t") } )
#' names(degs) <- files
#' methods <- c("Ternary", "Quaternary", "Enrichment", "Fisher")
#' enrichment <- runInferenceModels(NULL, NULL, DEGs=degs,
#'                                 method = methods,
#'                                 ents=ChIP1ap$filteredChIP.ents,
#'                                 rels=ChIP1ap$filteredChIP.rels,
#'                                 useFile=F, useMart=TRUE, useBHLH=TRUE,
#'                                 martFN="../CIE/data/mart_human_TFs.csv",
#'                                 BHLHFN="../CIE/data/BHLH_TFs.txt")
#' createCytoGraph(enrichment, ChIP1ap$filteredChIP.ents, ChIP1ap$filteredChIP.rels, degs)
    
createCytoGraph <- function(enrichment, ents, rels, DEGs) {
    if(class(ents) != "data.frame" &&
       class(rels) != "data.frame") {
        stop("Please make sure that the provied ents and rels tables are data.frames")
    }
    if(class(enrichment) == "data.frame" &&
       class(DEGs) == "data.frame") {
        createCytoGraphHelper(enrichment, ents, rels, DEGs)
    }
    else if(class(enrichment) == "list" &&
            class(enrichment[[1]]) == "list" &&
            class(DEGs) == "list") {
        lapply(1:length(enrichment), function(x) {
            lapply(1:length(DEGs), function(y) {
                createCytoGraphHelper(enrichment[[x]][[y]],
                                      ents, rels, DEGs[[y]],
                                      names(enrichment)[x],
                                      names(enrichment[[x]])[y]) } ) } )
    }
    else if(class(enrichment) == "list" &&
            class(enrichment[[1]]) == "data.frame" &&
            class(DEGs) == "list") {
        lapply(1:length(enrichment), function(x) {
            createCytoGraphHelper(enrichment[[x]], ents, rels, DEGs[[x]],
                                  condition=names(enrichment)[x]) }  )
    }
}
createCytoGraphHelper <- function(enrichment, ents, rels, DEG, method=NA,
                                  condition=NA) {
    pValIndicesEnr <- grep("pval|pvalue|p-val|p-value|p.val|p.value",
                           colnames(enrichment), ignore.case=TRUE)
    pValIndicesEnr <- unlist(pValIndicesEnr)
    sigProt <- enrichment[sapply(enrichment[,pValIndicesEnr], function(x) {
        x < 0.05} ), ]
    if(nrow(sigProt) == 0) {
        returnString <- paste("No significant proteins found.",
                              "For Method: ", method,
                              " and condition: ", condition)
        print(returnString)
        return
    }
    else {
        returnString <- paste("Significant proteins found.",
                              "For Method: ", method,
                              " and condition: ", condition,
                              "Graphing...")
        print(returnString)
    
        sigProt <- sigProt$uid
    
        sigRels <- rels %>% dplyr::filter(srcuid %in% sigProt)
        
        sigEnts <- ents[sigRels$trguid,]
        sigEnts <- sigEnts[complete.cases(sigEnts),]
        
        pValIndicesDEG <- grep("pval|pvalue|p-val|p-value|p.val|p.value",
                               colnames(DEGs), ignore.case=TRUE)
        pValIndicesDEG <- unlist(pValIndicesDEG)
        
        
        pValIndicesDEG <- grep("pval|pvalue|p-val|p-value|p.val|p.value",
                               colnames(DEGs), ignore.case=TRUE)
        pValIndicesDEG <- unlist(pValIndicesDEG)
        
        id.ind = grep('id|entr|Entrez', colnames(DEGs), ignore.case = T)
        
        sigEnts <- sigEnts %>% dplyr::filter(id %in% DEGs[[id.ind]])
        sigEnts <- rbind(ents[sigProt,], sigEnts)
        sigEnts <- sigEnts[complete.cases(sigEnts),]
        
        sigRels <- sigRels %>% dplyr::filter(trguid %in% sigEnts$uid)
        
        
        colorPal <- brewer.pal(3, "Dark2")
        colors <- as.character(sigEnts$type)
        colors[colors=="Protein"] <- colorPal[1]
        colors[colors=="mRNA"] <- colorPal[3]
        
        nodeD <- data.frame(id=as.character(sigEnts$uid),
                        name=sigEnts$name,
                        type=sigEnts$type,
                        stringsAsFactors=FALSE,
                        color=colors)
        edgeD <- data.frame(id = sigRels$uid,
                        source=as.character(sigRels$srcuid),
                        target=as.character(sigRels$trguid),
                        stringsAsFactors=FALSE)
    
    
        network <- createCytoscapeJsNetwork(nodeD, edgeD, edgeColor=colorPal[2])
        rcytoscapejs(network$nodes, network$edges, showPanzoom=FALSE)
        }
}
