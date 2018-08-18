require(dplyr)
require(readr)
require(org.Hs.eg.db)
source('./parseKBs.R')

## Processing string KB
ents.file     <- './origKB/string.ents'
rels.file     <- './origKB/string.rels'
L <- parseStringKB(ents.file, rels.file, verbose = T)
ents <- L$ents
rels <- L$rels

write.table(ents, './KB/string.ents', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(rels, './KB/string.rels', col.names = T, row.names = F, sep = '\t', quote = F)

## Processing TRED KB
in.file     <- './origKB/TRED.csv'
L <- parseTRED(in.file)
ents <- L$ents
rels <- L$rels

write.table(ents, './KB/TRED.ents', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(rels, './KB/TRED.rels', col.names = T, row.names = F, sep = '\t', quote = F)


## Processing ChIP-seq data
tsv_dir <- '../../CHIP/TSVs/TSVs/'
cutoff <- 50

L <- parseChIP(tsv_dir, cutoff)
ents <- L$ents
rels <- L$rels

write.table(ents, '../../KB/ChIP.ents', col.names = T, row.names = F, sep = '\t', quote = F)
write.table(rels, '../../KB/ChIP.rels', col.names = T, row.names = F, sep = '\t', quote = F)

