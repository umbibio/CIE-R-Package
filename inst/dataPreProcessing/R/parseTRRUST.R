
#' parseTRRUST
#' @description Process TRRUST data into a format usable by this enrichment pipeline. 
#' 
#' @usage parseTRED(in.file)
#' 
#' @param in.file A tab seperated file containing the TRRUST database.
#' 
#' @return A list of two dataframes, ents and rels.  The ents dataframe is all the entries from
#' the selected database. The rels table describes the relations between the entries
#' 
#' @export
#'
#' @examples
#'
#' entsRels <- parseTRRUST("trrust_rawdata.human.tsv")
#' ents <- entsRels$ents
#' rels <- entsRels$rels
#' 


parseTRRUST <- function(in.file){
  
  ## TRRUST Database
  TRRUST <- read_tsv(in.file)
  colnames(TRRUST) <- c("TF", "Target_genes", "type", "other")
  
  mRNAs <- unique(na.omit(TRRUST$Target_genes))
  TFs <- unique(na.omit(TRRUST$TF))
  
  TRRUST.ents.TF <- data.frame(name = TFs, type = Rle('Protein', length(TFs)),stringsAsFactors = F)
  TRRUST.ents.mRNA <- data.frame(name = mRNAs, type = Rle('mRNA', length(mRNAs)),stringsAsFactors = F)
  
  ## converting ids
  XX <- select(org.Hs.eg.db, keys=TRRUST.ents.TF$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  YY <- select(org.Hs.eg.db, keys=TRRUST.ents.mRNA$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  YY <- YY %>% distinct(SYMBOL, .keep_all = T)
  TRRUST.ents.TF <- left_join(TRRUST.ents.TF, XX, by = c('name' = 'SYMBOL'))
  TRRUST.ents.mRNA <- left_join(TRRUST.ents.mRNA, YY, by = c('name' = 'SYMBOL'))
 
  TRRUST.ents <- rbind(TRRUST.ents.TF, TRRUST.ents.mRNA)
  TRRUST.ents <- TRRUST.ents %>% mutate(uid = 1:nrow(TRRUST.ents)) %>%
    transmute(uid = uid, name = name, id = ENTREZID, type = type)
  
  
  TRRUST.ents.TF <- TRRUST.ents %>% filter(type == 'Protein')
  TRRUST.ents.mRNA <- TRRUST.ents %>% filter(type == 'mRNA')
  
  TRRUST.rels <- TRRUST %>% na.omit() %>% group_by(TF) %>%
    mutate(srcuid = TRRUST.ents.TF$uid[match(TF, TRRUST.ents.TF$name)]) %>%
    mutate(trguid = TRRUST.ents.mRNA$uid[match(Target_genes, TRRUST.ents.mRNA$name)]) %>%
    ungroup()
  TRRUST.rels <- TRRUST.rels[, c("srcuid","TF","Target_genes","trguid", "type", "other")]
  TRRUST.rels <- TRRUST.rels %>% transmute(uid = 1:nrow(TRRUST.rels),srcuid = srcuid, trguid = trguid, type = type, pmids = NA, nls = NA)
  
  TRRUST.rels$type <- gsub("Unknown", "conflict", TRRUST.rels$type)
  TRRUST.rels$type <- gsub("Activation", "increase", TRRUST.rels$type)
  TRRUST.rels$type <- gsub("Repression", "decrease", TRRUST.rels$type)
  
  L <- list(ents = TRRUST.ents, rels = TRRUST.rels)
  return(L)
  
  }




