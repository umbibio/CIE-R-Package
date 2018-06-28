library(dplyr)
library(org.Hs.eg.db)

require(readr)
parseTRUST <- function(in.file){
  
  ## TRUST Database
  TRUST <- read_tsv(in.file)
  colnames(TRUST) <- c("TF", "Target_genes", "type", "other")
  
  mRNAs <- unique(na.omit(TRUST$Target_genes))
  TFs <- unique(na.omit(TRUST$TF))
  
  TRUST.ents.TF <- data.frame(name = TFs, type = Rle('Protein', length(TFs)),stringsAsFactors = F)
  TRUST.ents.mRNA <- data.frame(name = mRNAs, type = Rle('mRNA', length(mRNAs)),stringsAsFactors = F)
  
  ## converting ids
  XX <- select(org.Hs.eg.db, keys=TRUST.ents.TF$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  YY <- select(org.Hs.eg.db, keys=TRUST.ents.mRNA$name, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
  YY <- YY %>% distinct(SYMBOL, .keep_all = T)
  TRUST.ents.TF <- left_join(TRUST.ents.TF, XX, by = c('name' = 'SYMBOL'))
  TRUST.ents.mRNA <- left_join(TRUST.ents.mRNA, YY, by = c('name' = 'SYMBOL'))
 
  TRUST.ents <- rbind(TRUST.ents.TF, TRUST.ents.mRNA)
  TRUST.ents <- TRUST.ents %>% mutate(uid = 1:nrow(TRUST.ents)) %>%
    transmute(uid = uid, name = name, id = ENTREZID, type = type)
  
  
  TRUST.ents.TF <- TRUST.ents %>% filter(type == 'Protein')
  TRUST.ents.mRNA <- TRUST.ents %>% filter(type == 'mRNA')
  
  TRUST.rels <- TRUST %>% na.omit() %>% group_by(TF) %>%
    mutate(srcuid = TRUST.ents.TF$uid[match(TF, TRUST.ents.TF$name)]) %>%
    mutate(trguid = TRUST.ents.mRNA$uid[match(Target_genes, TRUST.ents.mRNA$name)]) %>%
    ungroup()
  TRUST.rels <- TRUST.rels[, c("srcuid","TF","Target_genes","trguid", "type", "other")]
  TRUST.rels <- TRUST.rels %>% transmute(uid = 1:nrow(TRUST.rels),srcuid = srcuid, trguid = trguid, type = type, pmids = NA, nls = NA)
  
  TRUST.rels$type <- gsub("Unknown", "conflict", TRUST.rels$type)
  TRUST.rels$type <- gsub("Activation", "increase", TRUST.rels$type)
  TRUST.rels$type <- gsub("Repression", "decrease", TRUST.rels$type)
  
  L <- list(ents = TRUST.ents, rels = TRUST.rels)
  return(L)
  
  }




