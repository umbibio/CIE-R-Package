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
