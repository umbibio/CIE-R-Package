processDEGs <- function(DEGs, ents, rels, p.thresh = 0.05, fc.thresh = log(1.5)){
  ents.mRNA = ents[which(ents$type == 'mRNA'),]
  evidence = DEGs
  pval.ind = grep('qval|q.val|q-val|q-val|P-value|P.value|pvalue|pval|Pval', colnames(evidence), ignore.case = T)
  fc.ind = grep('fc|FC|fold|FoldChange', colnames(evidence), ignore.case = T)
  id.ind = grep('id|entr|Entrez', colnames(evidence), ignore.case = T)
  
  if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
    print('Please make sure the expression files column names are labled as entrez, fc, pvalue')
    quit(save = "no", status = 1, runLast = FALSE)
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
