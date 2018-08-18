require(readr)
require(dplyr)
datasets.file <- '../../GAP/datasets.csv'


datasets <- read.csv(datasets.file,header = F)
colnames(datasets) <- c('id', 'src', 'method', 'celltype')

DEGs.dir <- '../../GAP/DEGs/'
DEG.files <- list.files(DEGs.dir)

DEG.dat <- data.frame(id = character(), src = character(), target = character(), 
                      fc = double(), pval = double(), qval = double(), 
                      method = character(), celltype = character(),
                      stringsAsFactors = F)
for(f in DEG.files){
  tmp <- read_tsv(paste(DEGs.dir, f, sep = ''),col_names = F)
  if(nrow(tmp) == 0)
    next;
  colnames(tmp) <- c('target', 'caseMean','controlMean', 'tStat', 'fc', 'pval', 'qval', 
                    paste('other', 1:(ncol(tmp) - 7), sep = ''))
  id <- datasets$id[match(gsub('.txt', '', f), datasets$id)]
  src <- datasets$src[match(gsub('.txt', '', f), datasets$id)]
  method <- datasets$method[match(gsub('.txt', '', f), datasets$id)]
  celltype <- datasets$celltype[match(gsub('.txt', '', f), datasets$id)]
  tmp <- tmp %>% dplyr::select(target, fc, pval, qval) %>% 
    transmute(id = id, src = src, target = target, fc = fc, pval = pval, 
              qval = qval, method = method, celltype = celltype)
  DEG.dat <- rbind(DEG.dat, tmp)
}

DEG.dat <- DEG.dat %>% transmute(id = as.character(id), src = as.character(src), 
                                 target = as.character(target), fc = as.double(fc), 
                                 pval = as.double(pval), qval = as.double(qval),
                                 method = as.character(method), celltype = as.character(celltype))
## Fet Entrez IDs.
require(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
srcs <- unique(DEG.dat$src)
trgs <- unique(DEG.dat$target)

XX <- select(org.Hs.eg.db, key = trgs, columns = c('SYMBOL', 'ENTREZID'), keytype = 'SYMBOL')

#counts <- XX %>% group_by(SYMBOL) %>% summarise(n = n()) 
#XX1 <- left_join(XX, counts, by = 'SYMBOL') %>% filter(n > 1)

XX <- XX %>% distinct(SYMBOL, .keep_all = T)


YY <- select(org.Hs.eg.db, key = srcs, columns = c('SYMBOL', 'ENTREZID'), keytype = 'SYMBOL')

#counts <- XX %>% group_by(SYMBOL) %>% summarise(n = n()) 
#XX1 <- left_join(XX, counts, by = 'SYMBOL') %>% filter(n > 1)

DEG.dat <- left_join(DEG.dat, XX, by = c('target' = 'SYMBOL')) %>% rename('ENTREZID' = 'trg_id')
DEG.dat <- left_join(DEG.dat, YY, by = c('src' = 'SYMBOL')) %>% rename('ENTREZID' = 'src_id')

hs.TSs <- read.csv('../../KB/mart_human_TFs.csv')
bHLH.TFs <- read.table('../../KB/BHLH_TFs.txt', sep = '\t', header = T, stringsAsFactors = F)

DEG.dat <- DEG.dat %>% mutate(isTF = src %in% hs.TSs$HGNC.symbol)
DEG.dat <- DEG.dat %>% mutate(isBHLH = src %in% bHLH.TFs$Approved.Symbol)
DEG.dat <- DEG.dat %>% mutate(isDEG = ifelse(qval <= 0.01 & abs(fc) >= 1.5, T, F))


hasDEG <- DEG.dat %>% group_by(src) %>% summarise(n = sum(isDEG)) %>% filter(n > 1)

DEG.dat.hasDEG <- DEG.dat %>% filter(src %in% hasDEG$src & isDEG) 

## generate single perturbation DEG files for quaternaryprod usage
DEG.dat.hasDEG.isTF <- DEG.dat.hasDEG %>% filter(isTF)
#out.dir <- '../../SinglePertubation/'
#for(gene in unique(DEG.dat.hasDEG.isTF$src)){
#  file.name = paste(out.dir, gene, '.txt', sep = '')
#  tmp <- DEG.dat.hasDEG %>% filter(src == gene) %>% transmute(Entrez = trg_id, pvalue = pval, FoldChange = fc)
#  write.table(tmp, file.name, col.names = T, row.names = F, quote = F, sep = '\t')
#}

source("../CIE/inferenceModels/R/runCIE.R")
source("../CIE/inferenceModels/R/runCytoscape.R")
source("../CIE/networkAssembly/R/filterChIP.R")

ChIP5ave500 <- filterChIPAtlas(5, 500, "average", writeToFile = FALSE, databaseDir = './data/')
ents <- ChIP5ave500$filteredChIP.ents
rels <- ChIP5ave500$filteredChIP.rels

ents.tf <- ents %>% dplyr::filter(type == 'Protein')
ents.mRNA <- ents %>% dplyr::filter(type == 'mRNA')

rels <- left_join(rels, ents.tf, by = c('srcuid' = 'uid'))
rels <- left_join(rels, ents.mRNA, by = c('trguid' = 'uid'))
ChIP.targs <- rels %>% dplyr::transmute(src = name.x, chip.trg = name.y)
GPA <- DEG.dat.hasDEG.isTF %>% dplyr::transmute(src = src, gpa.trg = target)
ChIP.GPA <- inner_join(ChIP.targs, GPA, by = 'src') %>% arrange(src) %>% group_by(src) %>%
  summarise(n1 = length(unique(chip.trg)), n2 = length(unique(gpa.trg)), n3 = sum(unique(chip.trg) %in% unique(gpa.trg)))
