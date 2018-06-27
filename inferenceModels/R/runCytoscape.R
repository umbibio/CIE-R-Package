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
#' @param numProt The number of protiens, ranked by p.value from enrichment, to include
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
    
createCytoGraph <- function(enrichment, ents, rels, DEGs, p.thresh=0.05,
                            fc.thresh = log(1.5), numProt=5) {
    print("Starting Graph Generation")
    if(class(ents) != "data.frame" &&
       class(rels) != "data.frame") {
        stop("Please make sure that the provied ents and rels tables are data.frames")
    }
    if(class(enrichment) == "data.frame" &&
       class(DEGs) == "data.frame") {
        print("Data type detected")
        createCytoGraphHelper(enrichment, ents, rels, DEGs, p.thresh=p.thresh,
                              fc.thresh = fc.thresh, numProt=numProt)
    }
    else if(class(enrichment) == "list" &&
            class(enrichment[[1]]) == "list" &&
            class(DEGs) == "list") {
        print("Data type detected")
        lapply(1:length(enrichment), function(x) {
            lapply(1:length(DEGs), function(y) {
                createCytoGraphHelper(enrichment = enrichment[[x]][[y]],
                                      ents = ents, rels = rels, DEG = DEGs[[y]],
                                      p.thresh = p.thresh, fc.thresh = fc.thresh,
                                      method = names(enrichment)[x],
                                      condition = names(enrichment[[x]])[y],
                                      numProt = numProt) } ) } )
    }
    else if(class(enrichment) == "list" &&
            class(enrichment[[1]]) == "data.frame" &&
            class(DEGs) == "list") {
        print("Data type detected")
        lapply(1:length(enrichment), function(x) {
            createCytoGraphHelper(enrichment = enrichment[[x]],
                                  ents = ents, rels = rels, p.thresh = p.thresh,
                                  fc.thresh=fc.thresh, DEG= DEGs[[x]],
                                  condition=names(enrichment)[x], numProt=numProt) }  )
    }
}
createCytoGraphHelper <- function(enrichment, ents, rels, DEG,
                                  p.thresh, fc.thresh, method=NA,
                                  condition=NA, numProt) {

    print("Starting Analysis")
    sigProt <- enrichment[1:numProt,]
    if(nrow(sigProt) == 0) {
        returnString <- paste("No significant proteins found.",
                              "For Method: ", method,
                              " and condition: ", condition)
        print(returnString)
        return
    }  else {
        returnString <- paste("Significant proteins found.",
                              "For Method: ", method,
                              " and condition: ", condition,
                              "Graphing...")
        print(returnString)
    
        sigProt <- sigProt$uid
        
        ## Code written by Dr. Kourosh Zarringhalam
        pval.ind = grep('qval|q.val|q-val|q-val|P-value|P.value|pvalue|pval|Pval',
                        colnames(DEG), ignore.case = T)
        fc.ind = grep('fc|FC|fold|FoldChange', colnames(DEG), ignore.case = T)
        id.ind = grep('id|entr|Entrez', colnames(DEG), ignore.case = T)
  
        if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
            stop('Please make sure the expression files column names are labled as entrez, fc, pvalue')
        }
  
        colnames(DEG)[pval.ind] <- 'pvalue'
        colnames(DEG)[fc.ind] <- 'foldchange'
        colnames(DEG)[id.ind] <- 'id'

        sigDEG <- DEG %>% filter(abs(foldchange) >= fc.thresh & pvalue <= p.thresh) %>%
            transmute(id = id, val = ifelse(foldchange > 0, 1, -1), pval = pvalue) %>%
            distinct(id, .keep_all = T)

        ## End Dr. Zarringhalam code
        
        sigEnts <- ents[ents$id %in% sigDEG$id,]
        sigEntsTempUIDs <- lapply(sigProt, function(x) {
                    targets <- rels %>% dplyr::filter(srcuid == x) 
                    sigEntsTemp <- sigEnts[sigEnts$id %in% targets$trguid, ]
                    if(nrow(sigEntsTemp) > 10) {
                        sigDEGsub <- sigDEG[sigDEG$id %in% sigEntsTemp$id,] %>%
                            dplyr::arrange(pval)
                        sigDEGsub <- sigDEGsub[1:10,]
                        sigEntsTemp <- sigEntsTemp[sigEntsTemp$id %in% sigDEGsub$id,]
                        uids <- as.numeric(sigEntsTemp$uid)
                        uids
                    } else { sigEntsTemp$uid } } )
        sigEntsTempUIDs <- unique(unlist(sigEntsTempUIDs))
        sigEnts <- ents[sigEntsTempUIDs,]
                
        
        if(nrow(sigEnts) == 0) {
            warning("No significantly differentially expressed genes found as targets of proteins specified, using top 10 by p value")
            sigEntsTempUIDs <- lapply(sigProt, function(x) {
                uids <- rels %>% dplyr::filter(srcuid == x) %>%
                    dplyr::select(trguid)
                tempEnts <- ents[as.numeric(uids$trguid),]
                tempEnts <- tempEnts[tempEnts$id %in% DEG$id, ]
                uids <- as.numeric(tempEnts$uids)
                if(length(uids) > 10) {
                    pVal <- sapply(tempEnts$id, function(x) {
                        index <- which(DEG$id == x)
                        DEG$pVal[index]
                    })
                    print(head(pVal))
                    uids <- data.frame(uids, pVal)
                    uids <- uids %>% dplyr::arrange(pVal)
                    uids$uids[1:10]
                } else { uids } } )
            sigEntsTempUIDs <- unique(unlist(sigEntsTempUIDs))
            sigEnts <- ents[sigEntsTempUIDs,]
        } 
        
        sigRels <- rels[(rels$srcuid %in% sigProt), ]
        sigRels <- sigRels[(sigRels$trguid %in% sigEnts$uid),]

        sigEnts <- rbind(ents[sigProt,], sigEnts)
        
        mRNAfc <- sapply(sigEnts$id, function(x) {
            index <- which(sigDEG$id == x)
            if(length(index) == 0) {
                NA
            } else {
                sigDEG$val[index]
            }
        })
     
        colorPal <- brewer.pal(11, "Spectral")
        type <- as.character(sigEnts$type)
        
        colors <- sapply(1:nrow(colorMap), function(x) {
            if(type[x] == "mRNA") {
                if(is.na(mRNAfc[x])) { "#7b7c7c" }
                else if(mRNAfc[x] == 1) { colorPal[2] }
                else if(mRNAfc[x] == 0) { "#FFFFFF" }
                else if(mRNAfc[x] == -1) { colorPal[10] }
            } else { colorPal[11] } } )

        nodeD <- data.frame(id=as.character(sigEnts$uid),
                        name=sigEnts$name,
                        type=sigEnts$type,
                        color=colors,
                        stringsAsFactors=FALSE)
        
        edgeD <- data.frame(id = sigRels$uid,
                        source=as.character(sigRels$srcuid),
                        target=as.character(sigRels$trguid),
                        stringsAsFactors=FALSE)

        network <- createCytoscapeJsNetwork(nodeD, edgeD, edgeColor=colorPal[4])
        rcytoscapejs(network$nodes, network$edges, showPanzoom=FALSE)
    }
}
