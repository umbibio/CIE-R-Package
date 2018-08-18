
#' parseTRED
#' @description Process TRED data into a format usable by this enrichment pipeline. 
#' 
#' @usage parseTRED(in.file)
#' 
#' @param in.file A csv file containing the TRED database.
#' 
#' @return A list of two dataframes, ents and rels.  The ents dataframe is all the entries from
#' the selected database. The rels table describes the relations between the entries
#' 
#' @export
#'
#' @examples
#'
#' entsRels <- parseTRED("TRED.csv")
#' ents <- entsRels$ents
#' rels <- entsRels$rels
#' 

parseTRED <- function(in.file){
  TRED <- read.csv(in.file)
  TRED.TF.ents <- TRED %>% dplyr::select(regulator_symbol, regulator_id) %>% distinct() %>%
    mutate(type = 'Protein')
  colnames(TRED.TF.ents) <- c('name', 'id', 'type')
  
  TRED.genes.ents <- TRED %>% dplyr::select(target_symbol, target_id) %>% distinct() %>%
    mutate(type = 'mRNA')
  colnames(TRED.genes.ents) <- c('name', 'id', 'type')
  
  TRED.ents <- rbind(TRED.TF.ents, TRED.genes.ents) 
  TRED.ents <- TRED.ents %>% mutate(uid = 1:nrow(TRED.ents))
  TRED.ents <- TRED.ents[, c('uid', 'name', 'id', 'type')]
  
  TRED.ents.protein <- TRED.ents %>% filter(type == "Protein")
  TRED.ents.mRNA <- TRED.ents %>% filter(type == "mRNA")
 
  reg.uid  <- TRED.ents.protein$uid[match(TRED$regulator_id, TRED.ents.protein$id)]
  targ.uid <- TRED.ents.mRNA$uid[match(TRED$target_id, TRED.ents.mRNA$id)]
  
  TRED <- TRED %>% mutate(srcuid = reg.uid, trguid = targ.uid)
  TRED.rels <- TRED %>% transmute(uid = 1:nrow(TRED), srcuid = srcuid, trguid = trguid, 
                                  type = 'conflict', pmids = 'NA', nls = 'NA')
  
  L <- list(ents = TRED.ents, rels = TRED.rels)
  
  return(L)
  
}

#' parseStringKB
#' @description Process STRINGdb data into a format usable by this enrichment pipeline. Assumes
#' that if a protein targets a protein it also targets the gene.
#' 
#' @usage parseTRED(in.file)
#' 
#' @param ents.file A tab separated file containing the entries of STRINGdb
#'
#' @param rels.file A tab separated file describing the interactions between the entries of
#' STRINGdb
#'
#' @param verbose A boolean value, indicating whether the output should include verbose messaging.
#' Defaults to FALSE.
#' 
#' @return A list of two dataframes, ents and rels.  The ents dataframe is all the entries from
#' the selected database. The rels table describes the relations between the entries
#' 
#' @export
#'
#' @examples
#'
#' entsRels <- parseStringKB("string.ents", "string.rels")
#' ents <- entsRels$ents
#' rels <- entsRels$rels
#' 

parseStringKB <- function(ents.file, rels.file, verbose = F){
  ## Generating the one-level network form the knowledge base
  ## Note: There are ' and # characters in KB that mess up the tables. 
  ## The following will clean them up.
  ents <- read.table(ents.file, header = TRUE, stringsAsFactors = FALSE, strip.white=TRUE, sep = '\t',
                     quote = NULL, comment.char = '')
  colnames(ents) = c('uid', 'name', 'id', 'type')
  ents <- na.omit(ents) # remove incomplete entries
  
  rels <- read.table(rels.file, header = TRUE, stringsAsFactors = FALSE, strip.white=TRUE, sep = '\t', 
                     quote = NULL, comment.char = '')
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'nls')
  
  cleanup <- function(x) {
    x <- gsub('\"', '', x)
    x <- gsub('\'', 'p', x)
    x <- gsub('#', '_', x)
    x <- gsub(' ', '', x)
    return(x) 
  } 
  
  ents <- as.data.frame(lapply(ents,cleanup), stringsAsFactors = F)
  rownames(ents) <- 1:nrow(ents)
  rels <- as.data.frame(lapply(rels,cleanup), stringsAsFactors = F)
  rownames(rels) = 1:nrow(rels)
  
  ## Convert uids to integers
  uid.orig = ents$uid
  uid.new = seq(1, length(ents$uid))
  id.map = data.frame(uid.orig = uid.orig, uid.new = uid.new, stringsAsFactors = F)
  ents$uid = uid.new
  rels$srcuid = uid.new[match(rels$srcuid, uid.orig)]
  rels$trguid = uid.new[match(rels$trguid, uid.orig)]
  
  
  ## Protein, Compound or mRNA entries
  ents = unique(ents[which(ents$type %in% c('mRNA', 'Protein', 'Compound')),])
  
  ## (unique) Protein/Compound entries
  ents.pc = unique(ents[which(ents$type %in% c('Protein', 'Compound')),])
  ## (unique) mRNA entries
  ents.mRNA = unique(ents[which(ents$type == 'mRNA'),])
  
  ## src has to be protein or compound
  rels = rels[which(rels$srcuid %in% ents.pc$uid),]
  rownames(rels) = 1:nrow(rels)
  ## target has to be mRNA
  rels = rels[which(rels$trguid %in% ents.mRNA$uid),]
  rownames(rels) = 1:nrow(rels)
  
  ## remove -1's if any
  ents = ents[(ents$id != -1 | is.na(ents$uid)),]
  rels = rels[which(rels$trguid %in% ents$uid & rels$srcuid %in% ents$uid), ] ##update rels
  ents = ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  rownames(ents) = 1:nrow(ents)
  
  ## Remove duplicates
  if (anyDuplicated(ents)) ents <- ents[!duplicated(ents),]
  if (anyDuplicated(rels)) rels <- rels[!duplicated(rels),]
  
  ## Check that relation types are valid
  #valid <- (rels$type %in% c("increase","decrease","conflict"))
  #if (!all(valid)) rels <- rels[valid,]
  rels$type[which(!(rels$type %in% c("increase","decrease")))] = 'conflict'
  
  if (anyDuplicated(rels[,c("srcuid","trguid")])) {
    ## Identify duplicated rows
    dup.rels = rels[duplicated(rels[,c(2,3)]), ]
    ## Take one example from each
    dup.rels.uniq = dup.rels[!duplicated(dup.rels[,c(2,3)]), ]
    
    for(i in 1:nrow(dup.rels.uniq)){
      ind = which(rels$srcuid == dup.rels.uniq$srcuid[i] & rels$trguid == dup.rels.uniq$trguid[i])
      if(all(c("increase","decrease") %in% unique(rels$type[ind]))){
        rels$type[ind] = 'conflict'
      }else if("increase" %in% unique(rels$type[ind])){
        rels$type[ind] = 'increase'
      }else if("decrease" %in% unique(rels$type[ind])){
        rels$type[ind] = 'decrease'
      }else{
        rels$type[ind] = 'conflict'
      }
    }
    rels$type[which(!(rels$type %in% c('increase', 'decrease', 'conflict')))] = 'conflict'
  }
  
  ## remove duplicated rels
  if (any(duplicated(rels[,c("srcuid","trguid", "type")]))) 
    rels <- rels[!duplicated(rels[,c("srcuid","trguid", "type")]),]
  
  ## Check that sources in rels are protein/compound and targets are mRNA
  # Note: this also takes care of removing lines in rels with uids not matched in ents
  ind.pc <- which(ents$type %in% c("Protein","Compound"))
  inds <- which(!(rels$srcuid %in% ents$uid[ind.pc]))
  if (length(inds)) rels <- rels[-inds,]
  ind.mRNA <- which(ents$type=="mRNA")
  inds <- which(!(rels$trguid %in% ents$uid[ind.mRNA]))
  if (length(inds)) rels <- rels[-inds,]
  
  
  if (verbose == TRUE) {
    cat("\n Processed network dimensions:")
    cat("\n ents:", dim(ents)[1])
    cat("\n rels:", dim(rels)[1])
  }
  
  L = list(ents = ents, rels = rels)
}


#' parseChIP
#' @description Process unfiltered ChIP Atlas data, pulled into a
#' folder with a script like the one that can be found in
#' ../inst/shellScripts/.  Can either return a list of data
#' frames suitable for filtering with filterChIPAtlas (regardless
#' of species) or .rels and .ents tables suitable for running
#' enrichment or other purposes.  
#' 
#' @usage parseChIP(tsv_dir, processToFilterLater=FALSE) 
#' 
#' @param tsv_dir The directory containing the data pulled from ChIP
#' Atlas.
#' @param processToFilterLater A switch controlling whether the output
#' will be .ents and .rels tables (FALSE) or a list of data frames
#' suitable for later filtering by filterChIPAtlas (TRUE).  Defaults
#' to FALSE
#' 
#' @return A list of two data frames (ents and rels) or a list of
#' dataframes suitable for later processing.
#' 
#' @export
#'
#' @examples
#' ## Both these examples assume you pulled the ChIP Atlas data
#' ## and it can be found at the location "../data/chip-atlas-10kb/"
#' 
#' ## Example One
#' chip-atlas-10kb <- parseChIP("../data/chip-atlas-10kb/",
#'                              processToFilterLater=TRUE)
#' saveRDS(chip-atlas-10kb, "../data/chip-atlas-10kb.rds")
#'
#' entsRels <- filterChIPAtlas(10, 500, "average",
#'                             databaseDir="../data")
#' View(entsRels$filteredChIP.ents)
#' View(entsRels$filteredChIP.rels)
#'
#' ## Example Two
#'
#' entsRels <- parseChIP("../data/chip-atlas-10kb/")
#' View(entsRels$ents)
#' View(entsRels$rels)
#' 

parseChIP <- function(tsv_dir, processToFilterLater=FALSE){
    if(!processToFilterLater) {
        tsv_files <- list.files(tsv_dir)
        ChIP <- data.frame(Target_genes = NA, TF = NA,
                           stringsAsFactors = F)
        for(f in tsv_files){
            tmp <- read_tsv(paste(tsv_dir,f, sep = ''))
            if(ncol(tmp) < 2) next
            tmp <- tmp %>% dplyr::select(Target_genes) %>%
                mutate(TF = f)
            ChIP <- rbind(ChIP, tmp)
        }
  
        mRNAs <- unique(na.omit(ChIP$Target_genes))
        TFs <- unique(na.omit(ChIP$TF))
  
        ChIP.ents.TF <- data.frame(name = TFs,
                                   type = Rle('Protein',
                                              length(TFs)),
                                   stringsAsFactors = F)
        ChIP.ents.mRNA <- data.frame(name = mRNAs,
                                     type = Rle('mRNA', length(mRNAs)),
                                     stringsAsFactors = F)
        
        ## converting ids
        XX <- select(org.Hs.eg.db, keys=ChIP.ents.TF$name,
                     columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
        
        YY <- select(org.Hs.eg.db, keys=ChIP.ents.mRNA$name,
                     columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  
  #YY.counts <- YY %>% group_by(SYMBOL) %>% summarise(n = n()) 
  #YY <- left_join(YY, YY.counts)
  #YY %>% filter(n > 1)
  
        YY <- YY %>% distinct(SYMBOL, .keep_all = T)
        ChIP.ents.TF <- left_join(ChIP.ents.TF, XX,
                                  by = c('name' = 'SYMBOL'))
        ChIP.ents.mRNA <- left_join(ChIP.ents.mRNA,
                                    YY, by = c('name' = 'SYMBOL'))
        
        ChIP.ents <- rbind(ChIP.ents.TF, ChIP.ents.mRNA)
        ChIP.ents <- ChIP.ents %>% mutate(uid = 1:nrow(ChIP.ents)) %>% 
            transmute(uid = uid, name = name, id = ENTREZID, type = type)
        
        ChIP.ents.TF <- ChIP.ents %>% filter(type == 'Protein')
        ChIP.ents.mRNA <- ChIP.ents %>% filter(type == 'mRNA')
        
        ChIP.rels <- ChIP %>% na.omit() %>% group_by(TF) %>% 
            mutate(srcuid = ChIP.ents.TF$uid[match(TF, ChIP.ents.TF$name)]) %>% 
            mutate(trguid = ChIP.ents.mRNA$uid[match(Target_genes, ChIP.ents.mRNA$name)]) %>%
            ungroup()
        ChIP.rels <- ChIP.rels %>%
            transmute(uid = 1:nrow(ChIP.rels),
                      srcuid = srcuid, trguid = trguid,
                      type = 'conflict', pmids = NA, nls = NA)
  
        L = list(ents = ChIP.ents, rels = ChIP.rels)
        L
    }
    else {
        tsvFiles <- list.files(tsv_dir)
        ChIPlist <- list(1:length(tsvFiles))
        index <- 1;
        for(f in tsvFiles) {
            ChIPlist[[index]] <- read_tsv(paste(tsv_dir, f, sep =""))
            names(ChIPlist)[index] <- strsplit(f, "\\.")[[1]][1]
            ChIPlist[[index]] <- cbind(names(ChIPlist)[index], ChIPlist[[index]])
            colnames(ChIPlist[[index]])[1] <- "TF"
            index <- index + 1
        }
        ChIPlist
    }
}
