require(QuaternaryProd)
require(fdrtool)
require(dplyr)

#' Run Causual Inference Enrichment
#'
#' @description Runs inference models selected on network selected, with differentially
#' expressed genes provided, and with methods requested.
#' 
#' @usage runCIE(databaseType = c("TRED", "string", "ChIP"), filter = FALSE,
#'                               DEGs, p.thresh = 0.05, fc.thresh=log(1.5),
#'                               methods,
#'                               filteredDataName=NA, ents=NA, rels=NA, useFile=TRUE,
#'                               useMart=FALSE, useBHLH=FALSE, martFN=NA, BHLHFN=NA,
#'                               hypTabs= c("1", "2"), verbose=T, databaseDir = NA)
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
#' @param DEGs Either a list of dataframes of differentally expressed gene data, or a
#' single data frame.  Must have columns consisting of p or q value, fold change, and
#' entrez id.
#'
#' @param p.thresh 0.05 by default, used as a cuttoff for filtering the DEGs by their
#' q or p value.
#'
#' @param fc.thres log(1.5) by default, used as cutoff for filtering the DEGs by their
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
#' @param hypTabs One of 1 or 2, indicating which hypTab function you wish to use
#'
#' @param verbose Make output verbose, true by default
#'
#' @param databaseDir The path where the .rels and .ents files for a given Database type
#' can be found.
#' 
#' @return Either a single data frame or list of data frames, depending on the number of
#' methods and DEGs provided, which contain the enrichment results (if a list, they are
#' named with the methods and the original names of the DEG list)
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
#' View(enrichment$Quaternary$cell_line_cxcl12_tgfb_evidence_edgeR.txt)
#' 

runCIE <- function(databaseType = c("TRED", "string", "ChIP"),
                               filter = FALSE,
                               DEGs, p.thresh = 0.05, fc.thresh=log(1.5),
                               methods,
                               filteredDataName=NA, ents=NA, rels=NA, useFile=TRUE,
                               useMart=FALSE, useBHLH=FALSE, martFN=NA, BHLHFN=NA,
                               hypTabs= c("1", "2"), verbose=T, databaseDir = NA) {
    hypTabs = match.arg(hypTabs)
    if(useFile) {
        if(!filter & is.na(filteredDataName)) {
            if(is.na(databaseDir)) {
                relsFN <- paste("../../data/", databaseType, ".rels", sep="")
                entsFN <- paste("../../data/", databaseType, ".ents", sep="")
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
    if(class(DEGs) == "list") {
        DEGs.E <- lapply(DEGs, function(x) {
            processDEGs(x, ents, rels, p.thresh, fc.thresh) } )
        names(DEGs.E) <- names(DEGs)
    }
    else {
        DEGs.E <- processDEGs(DEGs, ents, rels, p.thresh, fc.thresh)
    }
    print("Running Enrichment")
    if(length(methods) > 1 & class(DEGs.E) == "list") {
        enrichment <- lapply(methods, function(x) {
            lapply(DEGs.E, function(y) {
                runEnrichment(ents, rels, y, verbose, hypTabs, x) } ) } )
        names(enrichment) <- methods
    }
    else if(length(methods) > 1 & class(DEGs.E) == "data.frame") {
        enrichment <- lapply(methods, function(x) {
            runEnrichment(ents, rels, DEGs.E, verbose, hypTabs, x)
        } )
        names(enrichment) <- methods
    }
    else if(length(methods) == 1 & class(DEGs.E) == "list") {
        enrichment <- lapply(DEGs.E, function(x) {
            runEnrichment(ents, rels, x, verbose, hypTabs, methods) } )
    }
    else {
        enrichment <- runEnrichment(ents, rels, DEGs.E, verbose, hypTabs, method)
    }
    print("Complete!")
    return(enrichment)
}
runEnrichment <- function(ents, rels, DEGtable, verbose, hypTabs, method) {
    if(hypTabs == 1) {
        enrichment <- generateHypTabs(ents, rels, DEGtable, verbose=verbose,
                                      method = method)
        index <- grep("pval|pvalue|p.value|p-value|p-val|p.val", colnames(enrichment))
        index <- unlist(index)
        if(length(index) == 0) {
            stop("Enrichment analysis not returning expected results")
        }
        else if(length(index) == 1) {
            enrichment <- enrichment %>% arrange(.[[index]])
        }
        else {
            index2 <- sapply(c("adj", "up"), function(x) {grep(x, colnames(enrichment))})
            index2 <- unlist(index2)
            indexFinal <- index[!(index %in% index2)]
            enrichment <- enrichment %>% arrange(.[[indexFinal]])
        }
    }
    else {
        enrichment <- generateHypTabs2(ents,rels, DEGtable, verbose=verbose,
                                       method=method)
        index <- grep("pval|pvalue|p.value|p-value|p-val|p.val", colnames(enrichment))
        index <- unlist(index)
        if(length(index) == 0) {
            stop("Enrichment analysis not returning expected results")
        }
        else if(length(index) == 1) {
            enrichment <- enrichment %>% arrange(.[[index]])
        }
        else {
            index2 <- sapply(c("adj", "up"), function(x) {grep(x, colnames(enrichment))})
            index2 <- unlist(index2)
            indexFinal <- index[!(index %in% index2)]
            enrichment <- enrichment %>% arrange(.[[indexFinal]])
        }
    }

}

## The following protion of code was written by Dr. Kourosh Zarringhalam, presented here
## with minor edits. The bulk of these edits was changing the response to input that
## could not be processed from quitting R to stopping function execution with a message.
processDEGs <- function(DEGs, ents, rels, p.thresh = 0.05, fc.thresh = log(1.5)){
  ents.mRNA = ents[which(ents$type == 'mRNA'),]
  evidence = DEGs
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
  
  evidence <- evidence %>% filter(abs(foldchange) >= fc.thresh & pvalue <= p.thresh) %>%
    transmute(id = id, val = ifelse(foldchange > 0, 1, -1)) %>% distinct(id, .keep_all = T)
  
  n.e1 <- nrow(evidence)
  evidence <- evidence %>% filter(id %in% ents.mRNA$id)
  n.e2 = nrow(evidence)
  print(paste((n.e1-n.e2), "evidence removed!"))
  
  evidence <- rbind(evidence,
                    data.frame(id = ents.mRNA$id[!(ents.mRNA$id %in% evidence$id)], val = 0))
  
  ##Change id back to uid
  evidence <- left_join(evidence, ents.mRNA, by = 'id') %>%
    dplyr::select(uid, val)
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
                            method = c("Ternary","Quaternary","Enrichment","Fisher")){
  
  method = match.arg(method)
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]

  D <- left_join(rels, evidence, by = c('trguid' = 'uid'))
  D <- D %>% group_by(srcuid) %>% 
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
      total.sig.genes = sum(evidence$val != 0, na.rm = T))
  
  if(method %in% c("Enrichment","Fisher")){
    D <- D %>% rowwise() %>% 
      mutate(pval = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz, 
                           direction = 'up', method = method))
    
    D <- D %>% mutate(adj.pval = p.adjust(pval, method = 'fdr'))
    
  }else{
    D <- D %>% rowwise() %>% 
      mutate(
        pval.up = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz,
                         direction = 'up', method = method),
        pval.down = runCRE(npp, npm, npz, nmp, nmm, nmz, nrp, nrm, nrz, nzp, nzm, nzz,
                           direction = 'down', method = method))
    D <- D %>% mutate(adj.pval.up = p.adjust(pval.up,method = 'fdr'),
                        adj.pval.down = p.adjust(pval.down,method = 'fdr'))
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
