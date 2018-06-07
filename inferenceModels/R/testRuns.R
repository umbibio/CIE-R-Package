## Required libraires
require(QuaternaryProd)
require(fdrtool)
require(dplyr)
source('./enrichments.R')
source('./processDEGs.R')

## Using TRED
ents.file.TRED     <- '../KB/TRED.ents'
rels.file.TRED     <- '../KB/TRED.rels'

ents.file.STRING     <- '../KB/string.ents'
rels.file.STRING     <- '../KB/string.rels'



ents.TRED <- read.table(ents.file.TRED, sep = '\t', header = T, stringsAsFactors = F)
rels.TRED <- read.table(rels.file.TRED, sep = '\t', header = T, stringsAsFactors = F)

ents.STRING <- read.table(ents.file.STRING, sep = '\t', header = T, stringsAsFactors = F)
rels.STRING <- read.table(rels.file.STRING, sep = '\t', header = T, stringsAsFactors = F)


### Get TF Info
hs.TSs <- read.csv('../KB/mart_human_TFs.csv')
bHLH.TFs <- read.table('../../KB/BHLH_TFs.txt', sep = '\t', header = T, stringsAsFactors = F)
ents.KB <- ents.KB %>% mutate(isTF = name %in% hs.TSs$HGNC.symbol)
ents.KB <- ents.KB %>% mutate(isBHLH = name %in% bHLH.TFs$Approved.Symbol)
###
## Get differentailly expressed genes
## conditions:
## c1 : Vehicle
## c2 : CXCL12 treated
## c3 : TGFb treated

Vehicle.vs.CXCL12 <- read.table('../data/cell_line_vehicle_cxcl12_evidence_edgeR.txt', header = T, sep = '\t')
Vehicle.vs.TGFb   <- read.table('../data/cell_line_vehicle_tgfb_evidence_edgeR.txt', header = T, sep = '\t')
CXCL12.vs.TGFb    <- read.table('../data/cell_line_cxcl12_tgfb_evidence_edgeR.txt', header = T, sep = '\t')

VC.TRED.e <- processDEGs(Vehicle.vs.CXCL12, ents.TRED, rels.TRED)
VT.TRED.e <- processDEGs(Vehicle.vs.TGFb, ents.TRED, rels.TRED)
CT.TRED.e <- processDEGs(CXCL12.vs.TGFb, ents.TRED, rels.TRED)


VC.STRING.e <- processDEGs(Vehicle.vs.CXCL12, ents.STRING, rels.STRING)
VT.STRING.e <- processDEGs(Vehicle.vs.TGFb, ents.STRING, rels.STRING)
CT.STRING.e <- processDEGs(CXCL12.vs.TGFb, ents.STRING, rels.STRING)


## Running enrichment analysis
CT.STRING.quaternary <- generateHypTabs(ents.STRING, rels.STRING, CT.STRING.e, verbose=T, method = "Quaternary")

CT.STRING.quaternary <- CT.STRING.quaternary %>% arrange(pval.down)


## Running enrichment analysis
VC.TRED.fisher <- generateHypTabs(ents.TRED, rels.TRED, VC.TRED.e, verbose=T, 
                                  method = "Fisher")
VC.TRED.fisher <- VC.TRED.fisher %>% arrange(pvalue)

