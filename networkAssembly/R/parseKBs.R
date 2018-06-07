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

parseChIP <- function(tsv_dir, cutoff){
  tsv_files <- list.files(tsv_dir)
  ChiP <- data.frame(Target_genes = NA, TF = NA, stringsAsFactors = F)
  for(f in tsv_files){
    tmp <- read_tsv(paste(tsv_dir,f, sep = ''))
    if(ncol(tmp) < 2) next
    tmp <- tmp %>% filter(.[[2]] > cutoff) %>% dplyr::select(Target_genes) %>% mutate(TF = f)
    ChiP <- rbind(ChiP, tmp)
  }
  
  mRNAs <- unique(na.omit(ChiP$Target_genes))
  TFs <- unique(na.omit(ChiP$TF))
  
  ChiP.ents.TF <- data.frame(name = TFs, type = Rle('Protein', length(TFs)),stringsAsFactors = F)
  ChiP.ents.mRNA <- data.frame(name = mRNAs, type = Rle('mRNA', length(mRNAs)),stringsAsFactors = F)
  
  
  ## converting ids
  XX <- select(org.Hs.eg.db, keys=ChiP.ents.TF$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  
  YY <- select(org.Hs.eg.db, keys=ChiP.ents.mRNA$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  
  #YY.counts <- YY %>% group_by(SYMBOL) %>% summarise(n = n()) 
  #YY <- left_join(YY, YY.counts)
  #YY %>% filter(n > 1)
  
  YY <- YY %>% distinct(SYMBOL, .keep_all = T)
  ChiP.ents.TF <- left_join(ChiP.ents.TF, XX, by = c('name' = 'SYMBOL'))
  ChiP.ents.mRNA <- left_join(ChiP.ents.mRNA, YY, by = c('name' = 'SYMBOL'))
  
  ChiP.ents <- rbind(ChiP.ents.TF, ChiP.ents.mRNA)
  ChiP.ents <- ChiP.ents %>% mutate(uid = 1:nrow(ChiP.ents)) %>% 
    transmute(uid = uid, name = name, id = ENTREZID, type = type)
  
  ChiP.ents.TF <- ChiP.ents %>% filter(type == 'Protein')
  ChiP.ents.mRNA <- ChiP.ents %>% filter(type == 'mRNA')
  
  ChiP.rels <- ChiP %>% na.omit() %>% group_by(TF) %>% 
    mutate(srcuid = ChiP.ents.TF$uid[match(TF, ChiP.ents.TF$name)]) %>% 
    mutate(trguid = ChiP.ents.mRNA$uid[match(Target_genes, ChiP.ents.mRNA$name)]) %>%
    ungroup()
  ChiP.rels <- ChiP.rels %>%
    transmute(uid = 1:nrow(ChiP.rels), srcuid = srcuid, trguid = trguid, type = 'conflict', pmids = NA, nls = NA)
  
  L = list(ents = ChiP.ents, rels = ChiP.rels)
  return(L)
}