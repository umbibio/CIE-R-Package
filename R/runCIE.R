#' Run Causual Inference Enrichment
#'
#' @description Runs inference models selected on network selected, with differentially
#' expressed genes provided, and with methods requested.
#' 
#' @usage runCIE(databaseType = c("TRED", "string", "ChIP", "TRRUST"), filter = FALSE,
#'                               DGEs, p.thresh = 0.05, fc.thresh=log(1.5),
#'                               methods,
#'                               filteredDataName=NA, ents=NA, rels=NA, useFile=TRUE,
#'                               useMart=FALSE, useBHLH=FALSE, martFN=NA, BHLHFN=NA,
#'                               targetsOfInterest=NA
#'                               hypTabs= c("1", "2"), verbose=T, databaseDir = NA,
#'                               expectProgressObject=FALSE,
#'                               progress=NA, numCores=NA)
#' @param databaseType Currently we only support selecting one database at a time.  Using
#' this option and having the database automatically loaded for you only currently works
#' if your working directory is the folder where this script exists in the package.  You
#' may also override the path it uses to look for the databse with the parameter
#' databaseDir.  There  are ways to overide this and get the databse from the environment
#' (see parameters ents and rels). If not needed, enter NULL.
#'
#' @param filter If using the repositories that come with this function, and you have made
#' filter ChIP Atlas save the filtered data to file, use this option as TRUE.  If your data
#' is not in the data folder with the format ChIPfilter.rels|ents, please specify a full
#' path with parameter filteredDataName.  If you are not reading from file, enter NULL.
#'
#' @param DGEs Either a list of dataframes of differentally expressed gene data, or a
#' single data frame.  Must have columns consisting of p or q value, fold change, and
#' entrez id.
#'
#' @param p.thresh 0.05 by default, used as a cuttoff for filtering the DGEs by their
#' q or p value.
#'
#' @param fc.thres log(1.5) by default, used as cutoff for filtering the DGEs by their
#' fold change
#'
#' @param ents Optional argument, specifies a .ents table that you have in your environment
#' in the event that you did not wish to write the .rels or .ents tables to file
#'
#' @param rels Optional argument, specifies a .rels table that you have in your environment
#' in the event that you did not wish to write the .rels or .ents tables to file.
#' **NOTE** You cannot mix and match your sources of .rels and .ents tables at this time.
#' If supplying as data frames from environment, none can be read from file and visa versa
#'
#' @param useFile Boolean value, indicating the source of the database.  If using the
#' els and rels paramters, set this to FALSE.  Otherwise, the default is TRUE
#'
#' @param useMart Boolean value, indicating whether you are adding additional information
#' to your enrichment about whether or not an entry is TF.  If you are, you must supply a
#' valid martFN.  Default value for this is false.
#'
#' @param useBHLH Boolean value, indicating whether you are adding additional information
#' to your enrichment about whether or not this enry is a BHLH.  If you are, you must
#' supply a valid BHLHFN.  Default value for this is false
#'
#' @param martFN The full path to the file where the Mart Human TF data can be found
#'
#' @param BHLHFN The full path to the file where the BHLH information can be found
#'
#' @param targetsOfIntrest A vector of strings of targets that you want to
#' show whether or not a protien targets in the analysis.
#'
#' @param hypTabs One of 1 or 2, indicating which hypTab function you wish to use
#'
#' @param verbose Make output verbose, true by default
#'
#' @param databaseDir The path where the .rels and .ents files for a given Database type
#' can be found.
#'
#' @param expectProgressObject Boolean, a switch to use a shiny progress bar instead of
#' text progress bar.
#'
#' @param progress Optional parameter, specifying a shiny progress bar to be updated.
#'
#' @param numCores Optional parameter, specifies number of cores
#' 
#' @return Either a single data frame or list of data frames, depending on the number of
#' methods and DGEs provided, which contain the enrichment results (if a list, they are
#' named with the methods and the original names of the DGE list)
#' @export
#'
#' @examples
#' 
#' ## This example assumes that you have differential gene expression
#' ## files ending in "edgeR.txt" in the current working directory and
#' ## that the database directory is "../data".  The filtering
#' ## parameters for ChIP in this instance are tailored for the
#' ## experiment, namely restricting results to interactions observed
#' ## in prostatic tissues.  For information about filtering the ChIP
#' ## Atlas database, please visit our Wiki or review the
#' ## documentation for filterChIPAtlas.
#' ## Please note that if you ever
#' ## wish to run enrichment on multiple conditions you may do so in
#' ## the following manner and have your results be accessible by
#' ## [out var name]$[method, if more than one]$[name of dge table
#' ## named list]
#' ## ex. enrichment$cell_line_cxcl12_tgfb_evidence_edgeR.txt
#' ## ex. enrichment[[1]]
#' ## Or as shown below:
#' 
#' ChIP1ap <- filterChIPAtlas(1, NA, "auto", NA, "prostate", NA, NA, FALSE)
#' files <- list.files("./", "edgeR.txt")
#' dges <- lapply(files, function(x) {
#'    read.table(x, header=T, sep="\t") } )
#' names(dges) <- files
#' methods <- c("Enrichment", "Fisher")
#' enrichment <- runCIE(NULL, NULL, DGEs=dges,
#'                      method = methods,
#'                      ents=ChIP1ap$filteredChIP.ents,
#'                      rels=ChIP1ap$filteredChIP.rels,
#'                      useFile=F, useMart=TRUE, useBHLH=TRUE,
#'                      martFN="../CIE/data/mart_human_TFs.csv",
#'                      BHLHFN="../CIE/data/BHLH_TFs.txt")
#' View(enrichment$Fisher$cell_line_cxcl12_tgfb_evidence_edgeR.txt)
#' 

runCIE <- function(databaseType = c("TRED", "string", "ChIP", "TRRUST"),
                   filter = FALSE,
                   DGEs, p.thresh = 0.05, fc.thresh=log(1.5),
                   methods,
                   filteredDataName=NA, ents=NA, rels=NA, useFile=TRUE,
                   useMart=FALSE, useBHLH=FALSE, martFN=NA, BHLHFN=NA,
                   targetsOfInterest=NA,
                   hypTabs= c("1", "2"), verbose=T, databaseDir = NA,
                   expectProgressObject=FALSE,
                   progress=NA, numCores=NA) {
    hypTabs = match.arg(hypTabs)
    if(useFile) {
        if(!filter & is.na(filteredDataName)) {
            if(is.na(databaseDir)) {
                relsFN <- paste("../data/", databaseType, ".rels", sep="")
                entsFN <- paste("../data/", databaseType, ".ents", sep="")
            }
            else {
                relsFN <- paste(databaseDir, databaseType, ".rels", sep="")
                entsFN <- paste(databaseDir, databaseType, ".ents", sep="")
            }
        }
        else if(filter & databaseType != "ChIP") {
            stop("Filtering is currenlty only suppored for ChIP Atlas data")
        }
        else if(filter & is.na(filteredDataName)) {
            relsFN <- paste(databaseType, "filter.rels", sep="")
            entsFN <- paste(databaseType, "filter.ents", sep="")
        }
        else {
            relsFN <- paste(filteredDataName, ".rels", sep="")
            entsFN <- paste(filteredDataName, ".ents", sep="")
        }
        rels <- read.table(relsFN, sep="\t", header=T, stringsAsFactors=F)
        ents <- read.table(entsFN, sep="\t", header=T, stringsAsFactors=F)
    }
    entCols <- c("uid", "name", "id", "type")
    result <- sapply(entCols, function(x) {grep(x, colnames(ents))})
    result <- unlist(result)
    if(length(result) == 0) {
        stop("Please make sure that your ents table contains the correct columns. They should be uid, name, id, and type. For clarification, please see our Wiki")
    }
    relsCols <- c("uid","srcuid", "trguid", "type", "pmids", "nls")
    result2 <- sapply(relsCols, function(x) {grep(x, colnames(rels))})
    result2 <- unlist(result)
    if(length(result) == 0) {
        stop("Please make sure that your ents table contains the correct columns. They should be uid, srcuid, trguid, type, pmids, and nls. For clarification, please see our Wiki")
    }    

    if(useMart) {
        if(!is.na(martFN)){
            hs.TFs <- read.csv(martFN)
            ents <- ents %>% mutate(isTF = name %in% hs.TFs$HGNC.symbol)
        }
        else{
            stop("Please provide a file where the Mart Human TFs can be found")
        }
    }
    if(useBHLH) {
        if(!is.na(BHLHFN)) {
            BHLH.TFs <- read.table(BHLHFN, sep="\t", header=T, stringsAsFactors=F)
            ents <- ents %>% mutate(isBHLH = name %in% BHLH.TFs$Approved.Symbol)
        }
        else {
            stop("Please provide a file where the BHLH TFs can be found")
        }
    }
    if(!is.na(targetsOfInterest[1])) {
        if(length(targetsOfInterest) >= 1) {
            if(is.na(databaseDir)) {
                relsTarg <- readRDS("../CIE/data/relsForTarg.rds")
                entsTarg <- readRDS("../CIE/data/entsForTarg.rds")
            }
            else {
                relsTarg <- readRDS(paste(databaseDir, "relsForTarg.rds",
                                          sep=""))
                entsTarg <- readRDS(paste(databaseDir, "entsForTarg.rds",
                                          sep=""))
            }
            sigRels <- relsTarg %>%
                dplyr::filter(entsTarg$name[trguid] %in% targetsOfInterest)

            colTitle <- paste(targetsOfInterest, collapse = "_")
            colTitle <- paste("targets", colTitle, sep = "_")
            ents <- ents %>% dplyr::mutate(colTitle = ents$uid %in% sigRels$srcuid)
            colnames(ents)[which(colnames(ents) == "colTitle")] = colTitle
        }
    }
    if(class(DGEs) == "list") {
        DGEs.E <- lapply(DGEs, function(x) {
            processDGEs(x, ents, rels, p.thresh, fc.thresh, expectProgressObject) } )
        names(DGEs.E) <- names(DGEs)
    }
    else {
        DGEs.E <- processDGEs(DGEs, ents, rels, p.thresh, fc.thresh,
                              expectProgressObject)
    }
    print("Running Enrichment")
    if(length(methods) > 1 & class(DGEs.E) == "list") {
        enrichment <- lapply(methods, function(x) {
            lapply(DGEs.E, function(y) {
                runEnrichment(ents, rels, y, verbose, hypTabs, x,
                              expectProgressObject,
                              progress, numCores) } ) } )
        names(enrichment) <- methods
    }
    else if(length(methods) > 1 & class(DGEs.E) == "data.frame") {
        enrichment <- lapply(methods, function(x) {
            runEnrichment(ents, rels, DGEs.E, verbose, hypTabs, x,
                          expectProgressObject,
                          progress, numCores)
        } )
        names(enrichment) <- methods
    }
    else if(length(methods) == 1 & class(DGEs.E) == "list") {
        enrichment <- lapply(DGEs.E, function(x) {
            runEnrichment(ents, rels, x, verbose, hypTabs, methods, expectProgressObject,
                          progress, numCores) } )
    }
    else {
        enrichment <- runEnrichment(ents, rels, DGEs.E, verbose, hypTabs,
                                    methods, expectProgressObject,
                                    progress, numCores)
    }
    print("Complete!")
    return(enrichment)
}

runEnrichment <- function(ents, rels, DGEtable, verbose, hypTabs, method,
                          expectProgressObject, progress, numCores) {
    if(hypTabs == 1) {
        enrichment <- generateHypTabs(ents, rels, DGEtable, verbose=verbose,
                                      method = method,
                                      expectProgressObject=expectProgressObject,
                                      progress=progress,
                                      numCores=numCores)
        index <- grep("pval|pvalue|p.value|p-value|p-val|p.val", colnames(enrichment))
        index <- unlist(index)
        if(length(index) == 0) {
            stop("Enrichment analysis not returning expected results")
        }
        else if(length(index) == 1) {
            enrichment <- enrichment %>% arrange(.[[index]])
        }
        else {
            index2 <-grep("up", colnames(enrichment))
            index3 <- grep("adj", colnames(enrichment))
            indexFinal <- index[!(index %in% index2)]
            indexFinal <- indexFinal[!(indexFinal %in% index3)]
            enrichment <- enrichment %>% arrange(.[[indexFinal]])
        }
    }
    else {
        enrichment <- generateHypTabs2(ents,rels, DGEtable, verbose=verbose,
                                       method=method, numCores=numCores)
        index <- grep("pval|pvalue|p.value|p-value|p-val|p.val", colnames(enrichment))
        index <- unlist(index)
        if(length(index) == 0) {
            stop("Enrichment analysis not returning expected results")
        }
        else if(length(index) == 1) {
            enrichment <- enrichment %>% arrange(.[[index]])
        }
        else {
            index2 <-grep("up", colnames(enrichment))
            index3 <- grep("adj", colnames(enrichment))
            indexFinal <- index[!(index %in% index2)]
            indexFinal <- indexFinal[!(indexFinal %in% index3)]
            enrichment <- enrichment %>% arrange(.[[indexFinal]])
        }
    }

}
#' Pathway Enrichment
#' @description Gets enriched pathways from the Reactome database using their analysis
#' service
#' 
#' @usage pathwayEnrichment(sigProtiens, numPathways=10)
#' 
#' @param sigProtiens A vector of protein names, intended to be used with the names column
#' of the enrichment table returned by runCIE().  It is also suggested that you filter
#' the results of enrichment by p-value so that you do not have insignificant protiens
#' in your analysis.
#' 
#' @param numPathways The number of pathways to pull from Reactome, defaults to `10
#' 
#' @return A single data frame or list of dataframes  consisting of the pathway,
#' it's p-value and fdr from analysis
#' using the Reactome analysis service, and the list of proteins given in sigProtiens found
#' in the pathway.
#' 
#' @export
#'
#' @examples
#'
#' ## This example assumes that you have differential gene expression
#' ## files ending in "edgeR.txt" in the current working directory and
#' ## that the database directory is "../data".  The filtering
#' ## parameters for ChIP in this instance are tailored for the
#' ## experiment, namely restricting results to interactions observed
#' ## in prostatic tissues.  For information about filtering the ChIP
#' ## Atlas database, please visit our Wiki or review the
#' ## documentation for filterChIPAtlas.
#' 
#' ChIP1ap <- filterChIPAtlas(1, NA, "auto", NA, "prostate", NA, NA, FALSE, databaseDir="../data")
#' files <- list.files("./", "egeR.txt")
#' dges <- lapply(files, function(x) {
#'    read.table(x, header=T, sep="\t") } )
#' names(dges) <- files
#' methods <-  "Fisher"
#' enrichment <- runCIE(NULL, NULL, DGEs=dges,
#'                      method = methods,
#'                      ents=ChIP1ap$filteredChIP.ents,
#'                      rels=ChIP1ap$filteredChIP.rels,
#'                      useFile=F, useMart=TRUE, useBHLH=TRUE,
#'                      martFN="../data/mart_human_TFs.csv",
#'                      BHLHFN="../data/BHLH_TFs.txt",
#'                      databaseDir="../data")
#'
#' sigProt <- enrichment[[1]] %>%
#'   dplyr::filter(pval < 0.01) %>%
#'   dplyr::select(name)
#' sigProt <- sigProt$name
#' pathEnr <-pathwayEnrichment(sigProt)
#' View(pathEnr)


pathwayEnrichment <- function(sigProtiens, numPathways=10) {
    if(length(sigProtiens[[1]]) == 1) {
        pathEnr <- pathwayEnrichmentHelper(sigProtiens, numPathways)
        pathEnr
    }
    else if(length(sigProtiens[[1]]) > 1 &&
       length(sigProtiens[[1]][[1]]) == 1) {
        pathEnr <- lapply(sigProtiens, function(x) {
            pathwayEnrichmentHelper(x, numPathways)
        } )
        names(pathEnr) <- names(sigProtiens)
        pathEnr
    }
    else if(length(sigProtiens[[1]]) > 1 &&
            length(sigProtiens[[1]][[1]]) > 1 &&
            length(sigProtiens[[1]][[1]][[1]])) {
        pathEnr <- lapply(sigProtiens, function(x) {
            enrList <- lapply(x, function(y) {
                pathwayEnrichmentHelper(x, numPathways) } )
            names(enrList) <- names(x)
            enrList } )
        pathEnr
    }
    
}
pathwayEnrichmentHelper <- function(sigProtiens, numPathways) {
    ## write(sigProtiens, "proteins.txt")
    analysis <- fromJSON(system(paste("curl -H 'Content-Type: text/plain' --data-binary ",
                                      paste(sigProtiens, collapse=","),
                                      " --url https://reactome.org/",
                                      "AnalysisService/identifiers/projection/\\",
                                      "?pageSize\\=",
                                      numPathways, "\\&page\\=1", sep=""), intern=TRUE,
                                ignore.stderr=TRUE))
    ## system("rm proteins.txt")
    numPaths <- length(analysis$pathways)
    if(numPaths == 0) {
        print("None of the significant regulators you provided were found in Reactome by that Id")
        return(NA)
    }
    tableOut <- data.frame(id = sapply(1:numPaths,
                                       function(x) {analysis$pathways[[x]]$dbId}),
                           name = sapply(1:numPaths,
                                         function(x) {analysis$pathways[[x]]$name}),
                           pValue = sapply(1:numPaths,
                                           function(x) {
                                               analysis$pathways[[x]]$entities$pValue}),
                           fdr = sapply(1:numPaths,
                                        function(x) {analysis$pathways[[x]]$entities$fdr}))
    
    pathwayIds <- unlist(lapply(1:length(analysis$pathways),
                                function (x) { analysis$pathways[[x]][[2]] }))
    protiensFound <- sapply(pathwayIds, function(x) {
        reactions <- fromJSON(system(paste("curl https://reactome.org/AnalysisService/",
                                           "token/", analysis$summary$token,
                                           "/summary/", x,
                                           sep = ""), intern=TRUE,
                                     ignore.stderr=TRUE))
        ids <- unlist(lapply(reactions$identifiers, function(x) {x[1]}))
        sigProtiens[sigProtiens %in% ids] } )
    tableOut <- tableOut %>%
        dplyr::mutate(protiensFound = sapply(1:length(protiensFound),
                                             function(x) {
                                                 paste(protiensFound[[x]], collapse="; ")}))
    tableOut
    
}

## The following protion of code was written by Dr. Kourosh Zarringhalam, presented here
## with minor edits. The bulk of these edits was changing the response to input that
## could not be processed from quitting R to stopping function execution with a message.
processDGEs <- function(DGEs, ents, rels, p.thresh = 0.05, fc.thresh = log(1.5), expectProgressObject){
  ents.mRNA = ents[which(ents$type == 'mRNA'),]
  evidence = DGEs
  pval.ind = grep('qval|q.val|q-val|q-val|P-value|P.value|pvalue|pval|Pval',
                  colnames(evidence), ignore.case = T)
  fc.ind = grep('fc|FC|fold|FoldChange', colnames(evidence), ignore.case = T)
  id.ind = grep('id|entr|Entrez', colnames(evidence), ignore.case = T)
  
  if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
      stop('Please make sure the expression files column names are labled as entrez, fc, pvalue')
  }
  
  colnames(evidence)[pval.ind] <- 'pvalue'
  colnames(evidence)[fc.ind] <- 'foldchange'
  colnames(evidence)[id.ind] <- 'id'

  evidence  <- evidence %>% dplyr::filter(!is.na(id))
  
  evidence <- evidence %>% dplyr::filter(abs(foldchange) >= fc.thresh & pvalue <= p.thresh) %>%
    transmute(id = id, val = ifelse(foldchange > 0, 1, -1)) %>% distinct(id, .keep_all = T)
  
  n.e1 <- nrow(evidence)
  evidence <- evidence %>% filter(id %in% ents.mRNA$id)
  n.e2 = nrow(evidence)
  if(n.e2 == 0) {
      stop(paste("Please double-check to make sure the Entrez ids you have ",
                     "provided are mRNAs in the human genome (or genome you have",
                     " configured CIE to run with)", sep = ""))
  }
  print(paste((n.e1-n.e2), "evidence removed!"))
  map <- sapply(evidence$id, function(x) {
      which(ents.mRNA$id %in% x)})
  
  ##Change id back to uid
  evidence <- evidence %>%
      dplyr::mutate(uid = ents.mRNA$uid[map]) %>%
      dplyr::select(uid, val)
  
  evidence <- rbind(evidence,
                    data.frame(uid = ents.mRNA$uid[!(ents.mRNA$uid %in% evidence$uid)], val = 0))
  
  return(evidence)
}
    
    
runCRE <- function(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, 
                   direction = c('up', 'down'),
                   method = c("Ternary","Quaternary","Enrichment","Fisher")){
  
  direction = match.arg(direction)
  method = match.arg(method)
  
  if(method == 'Quaternary'){
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nrp + nzp
    nMinus = npm + nmm + nrm + nzm
    nZero  = npz + nmz + nrz + nzz
    
    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm + nrp + nrm - (npm + nmp)
    if(direction == 'up'){
      pval   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                            q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    }else if(direction == 'down'){
      ## Assume down-regulated
      qPlus  = nmp + nmm + nmz
      qMinus = npp + npm + npz
      score  = nmp + npm + nrp + nrm - (npp + nmm)
      pval   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                              q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    }
  }else if (method == 'Ternary'){
    qR     = 0
    qZero  = nzp + nzm + nzz
    nPlus  = npp + nmp + nzp
    nMinus = npm + nmm + nzm
    nZero  = npz + nmz + nzz
    
    ## Assume up-regulated
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    score  = npp + nmm - (npm + nmp)
    if(direction == 'up'){
      pval   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                         q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
    }else if(direction == 'down'){
      ## Assume down-regulated
      qPlus  = nmp + nmm + nmz
      qMinus = npp + npm + npz
      score  = nmp + npm - (npp + nmm)
      pval   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,q_r = qR, 
                         n_p = nPlus, n_m = nMinus, n_z = nZero)
    }
    
  } else if (method == 'Enrichment'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz
    
    qPlus  = 0
    qMinus = 0
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    
    nPlus  = nrp + nzp
    nMinus = nrm + nzm
    nZero  = nrz + nzz
    
    score  = nrp + nrm
    
    pval   = QP_Pvalue(score = score, q_p = qPlus, q_m = qMinus, q_z = qZero,
                       q_r = qR, n_p = nPlus, n_m = nMinus, n_z = nZero)
  } else if (method == 'Fisher'){
    nrp    = npp + nmp + nrp
    nrm    = npm + nmm + nrm
    nrz    = npz + nmz + nrz
    
    M = matrix(0, nrow = 2, ncol = 2)
    M[1,1] = nrp + nrm
    M[1,2] = nrz
    M[2,1] = nzp + nzm
    M[2,2] = nzz
    
    pval = fisher.test(M, alternative = 'greater')$p.val
  }
  
  return(pval)
}


generateHypTabs <- function(ents, rels, evidence, verbose=TRUE,
                            method = c("Ternary","Quaternary","Enrichment","Fisher"),
                            expectProgressObject, progress, numCores){

    if(is.na(numCores)) {
        numCores <- detectCores()
    }
    method = match.arg(method)
    ents.mRNA = ents[which(ents$type == 'mRNA'), ]
    
    D <- left_join(rels, evidence, by = c('trguid' = 'uid'))
    cluster <- create_cluster(cores=numCores, quiet=TRUE)
    intoGroups <- suppressWarnings(D %>% partition(srcuid, cluster=cluster))
    
    intoGroups %>%
        cluster_assign_value("evidence", evidence) %>%
        cluster_assign_value("D", D) %>%
        cluster_assign_value("rels", rels)
    if(!expectProgressObject) {
        pb <- txtProgressBar(min=0, max=10,
                             initial="0", char="=", width=NA,
                             style=3, file="")
        setTxtProgressBar(pb, 0)
    }
    else {
        value <- progress$getValue()
        value <- value + ((progress$getMax() - value) / 10)
        progress$set(message="Running Enrichment",
                             value=value)
    }
    D <- intoGroups %>% 
        summarise(
            npp = sum(val == 1 & type == 'increase', na.rm=T),
            npm = sum(val == -1 & type == 'increase', na.rm=T),
            npz = sum(val == 0 & type == 'increase', na.rm=T),
            nmp = sum(val == 1 & type == 'decrease', na.rm=T),
            nmm = sum(val == -1 & type == 'decrease', na.rm=T),
            nmz = sum(val == 0 & type == 'decrease', na.rm=T),
            nrp = sum(val == 1 & type == 'conflict', na.rm=T),
            nrm = sum(val == -1 & type == 'conflict', na.rm=T),
            nrz = sum(val == 0 & type == 'conflict', na.rm=T),
            nzp = sum(evidence$val[match(unique(D$trguid[!(D$trguid %in% trguid)]), evidence$uid)] == 1, na.rm = T),
            nzm = sum(evidence$val[match(unique(D$trguid[!(D$trguid %in% trguid)]), evidence$uid)] == -1, na.rm = T),
            nzz = sum(evidence$val[match(unique(D$trguid[!(D$trguid %in% trguid)]), evidence$uid)] == 0, na.rm = T),
            correct.pred = 
                sum((val == 1 & type == 'increase') | (val == -1 & type == 'decrease') 
                    | (val != 0 & type == 'conflict'), na.rm=T),
            incorrect.pred = 
                sum((val == -1 & type == 'increase') | (val == 1 & type == 'decrease') , na.rm=T),
            total.reachable = n(),
            significant.reachable = sum(val != 0, na.rm=T),
            total.ambiguous = length(type == 'conflict'),
            significant.ambiguous = sum(val != 0 & type == 'conflict', na.rm=T),
            unreachable = length(unique(D$trguid[!(D$trguid %in% trguid)])),
            total.genes = length(unique(rels$trguid)),
            total.sig.genes = sum(evidence$val != 0, na.rm = T)) %>%
        collect()
    rm(intoGroups)
    invisible(gc())
    if(!expectProgressObject) {
        setTxtProgressBar(pb, 5)
    }
    else {
        value <- progress$getValue()
        value <- value + 4*((progress$getMax() / 10))        
        progress$set(message="Calculating p-values",
                     value=value,
                     detail="This can take several minutes")
    }
    
    if(method %in% c("Enrichment","Fisher")){
      cluster2 <- parallel::makeCluster(numCores, type="PSOCK")
      registerDoParallel(cluster2)
      cluster2 %>% cluster_assign_value("runCRE", runCRE) %>%
          cluster_assign_value("D", D) %>%
          cluster_assign_value("method", method) %>%
          cluster_assign_value("QP_Pvalue", QP_Pvalue) %>%
          cluster_assign_value("fisher.test", fisher.test)
      
      D <- D %>% mutate(pval = foreach(i = 1:nrow(D), .combine = c) %dopar% {
          runCRE(D$npp[i], D$npm[i], D$npz[i], D$nmp[i], D$nmm[i], D$nmz[i],
                 D$nrp[i], D$nrm[i], D$nrz[i], D$nzp[i], D$nzm[i],
                 D$nzz[i], method = method) } )
      parallel::stopCluster(cluster2)
      registerDoParallel()
      invisible(gc())
      library(stats)
      D <- D %>% mutate(adj.pval = stats::p.adjust(pval, method = 'fdr'))

      if(!expectProgressObject) {
          setTxtProgressBar(pb, 10)
      }
      else {
          value <- progress$getMax()
          progress$set(message="Complete!",
                             value=value)
      }
    
    }else{
      cluster2 <- parallel::makeCluster(numCores, type="PSOCK")
      registerDoParallel(cluster2)
      cluster2 %>% cluster_assign_value("runCRE", runCRE) %>%
          cluster_assign_value("D", D) %>%
          cluster_assign_value("method", method) %>%
          cluster_assign_value("QP_Pvalue", QP_Pvalue)
      
      D <- D %>% mutate(pval.up = foreach(i = 1:nrow(D), .combine = c) %dopar% {
          runCRE(D$npp[i], D$npm[i], D$npz[i], D$nmp[i], D$nmm[i], D$nmz[i],
                 D$nrp[i], D$nrm[i], D$nrz[i], D$nzp[i], D$nzm[i],
                 D$nzz[i], direction = 'up', method = method) } )
      D <- D %>% mutate(pval.down = foreach(i = 1:nrow(D), .combine = c) %dopar% {
          runCRE(D$npp[i], D$npm[i], D$npz[i], D$nmp[i], D$nmm[i], D$nmz[i],
                 D$nrp[i], D$nrm[i], D$nrz[i], D$nzp[i], D$nzm[i],
                 D$nzz[i], direction = 'down', method = method) } )
      parallel::stopCluster(cluster2)
      registerDoParallel()
      invisible(gc())
      D <- D %>% mutate(adj.pval.up = p.adjust(pval.up, method = 'fdr'),
                        adj.pval.down = p.adjust(pval.down, method = 'fdr'))
      
      
      if(!expectProgressObject) {
          setTxtProgressBar(pb, 10)
      }
      else {
          value <- progress$getMax()
          progress$set(message="Complete!",
                             value=value)
      }
  }

  if(!expectProgressObject) {
      close(pb)
  }
  D <- inner_join(ents, D, by = c('uid' = 'srcuid'))
  
  return(D)
  
}
generateHypTabs2 <- function(ents, rels, evidence, verbose=TRUE,
                            method = c("Ternary","Quaternary","Enrichment","Fisher"))
{
  method = match.arg(method)
  
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  
  ## For each hypothesis, identify the children and non-children and thier evidence values
  u.hyps = unique(rels$srcuid)
  child.uid = lapply(u.hyps, function(x) rels$trguid[which(rels$srcuid == x)])
  child.sgn = lapply(u.hyps, function(x) ifelse(rels$type[which(rels$srcuid == x)] == 'increase',
                                                1, ifelse(rels$type[which(rels$srcuid == x)] == 'decrease', -1, 0)))
  
  child.val = lapply(child.uid, function(x) getGeneVals(x, evidence))
  
  non.child.uid = lapply(child.uid, function(x) unique(ents.mRNA$uid[which(!(ents.mRNA$uid %in% x))]))
  non.child.val = lapply(non.child.uid, function(x) getGeneVals(x, evidence))
  
  ## Get the data slices corresponding to each hypothesis
  child.id = lapply(child.uid, function(x) as.numeric(ents.mRNA$id[match(x,ents.mRNA$uid)])) ## to get the id
  
  if (verbose == TRUE)
    cat("\n Computing pvalues")
  results = data.frame(matrix(0, nrow  = 2 * length(u.hyps), ncol = 12), stringsAsFactors = F)
  colnames(results) = c('uid', 'name', 'regulation', 'correct.pred', 'incorrect.pred', 'score',
                        'total.reachable', 'significant.reachable', 'total.ambiguous', 'significant.ambiguous',
                        'unknow', 'pval')
  
  if (verbose == TRUE)
    cat("\n Total number of hypothesis to consider:", length(u.hyps))
  
  for(p.s in 1:length(u.hyps)){
    #cat('.')
    results[(2*(p.s-1)+1), 1] = u.hyps[p.s]
    results[(2*p.s), 1]       = u.hyps[p.s]
    results[(2*(p.s-1)+1), 2] = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*p.s), 2]       = ents$name[which(ents$uid == u.hyps[p.s])]
    results[(2*(p.s-1)+1), 3] = 'up'
    results[(2*p.s), 3]       = 'down'
    
    npp = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 1))
    npm = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == -1))
    npz = length(which(child.sgn[[p.s]] == 1 &  child.val[[p.s]] == 0))
    
    nmp = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 1))
    nmm = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == -1))
    nmz = length(which(child.sgn[[p.s]] == -1 &  child.val[[p.s]] == 0))
    
    nrp = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 1))
    nrm = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == -1))
    nrz = length(which(child.sgn[[p.s]] == 0 &  child.val[[p.s]] == 0))
    
    nzp = length(which(non.child.val[[p.s]] == 1))
    nzm = length(which(non.child.val[[p.s]] == -1))
    nzz = length(which(non.child.val[[p.s]] == 0))
    
    pval = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, method = method)
    
    qPlus  = npp + npm + npz
    qMinus = nmp + nmm + nmz
    qR     = nrp + nrm + nrz
    qZero  = nzp + nzm + nzz
    
    results[(2*(p.s-1)+1), 4]  = npp + nmm
    results[(2*(p.s-1)+1), 5]  = npm + nmp
    results[(2*(p.s-1)+1), 6]  = npp + nmm - (npm + nmp)
    results[(2*(p.s-1)+1), 7]  = qPlus + qMinus + qR
    results[(2*(p.s-1)+1), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*(p.s-1)+1), 9]  = qR
    results[(2*(p.s-1)+1), 10] = nrp + nrm
    results[(2*(p.s-1)+1), 11] = qZero
    results[(2*(p.s-1)+1), 12] = pval$pval.up
    
    results[(2*p.s), 4]  = nmp + npm
    results[(2*p.s), 5]  = npp + nmm
    results[(2*p.s), 6]  = nmp + npm - (npp + nmm)
    results[(2*p.s), 7]  = qPlus + qMinus + qR
    results[(2*p.s), 8]  = npp + npm + nmp + nmm + nrp + nrm
    results[(2*p.s), 9]  = qR
    results[(2*p.s), 10] = nrp + nrm
    results[(2*p.s), 11] = qZero
    results[(2*p.s), 12] = pval$pval.down
  }
  
  
  nhyps  <- nrow(results)/2
  scores <- results$score
  odd  <- seq.int(1L,by=2L,len=nhyps)
  even <- seq.int(2L,by=2L,len=nhyps)
  ind  <- ifelse(scores[even] > scores[odd], even, odd)
  results <- results[ind,]
  
  fdr = fdrtool(results$pval, statistic = 'pvalue', plot=F, verbose=F)
  results = cbind(results, fdr$qval)
  colnames(results)[ncol(results)] = 'FDR'
  
  return(results)
  
}
