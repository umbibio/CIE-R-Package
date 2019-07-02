library(CIE)

dge  <- read.table("/path/to/sample/inputs/sample_input_myc.txt",
                   header=T, sep="\t")

tissChIP  <- filterChIPAtlas(NA, NA, NA, cellLineType="all", tissueCorrect=TRUE,
                             databaseDir="path/to/tissueCorrectedChIP/")

tissChIP.myc  <- runCIE(DGEs=dge,
                        methods="Ternary",
                        ents=tissChIP$ChIPall.ents,
                        rels=tissChIP$ChIPall.rels,
                        useFile=FALSE)
write.table(tissChIP.myc, "documentation/tissChIPmyc.tsv", quote=F,
            sep="\t", row.names=FALSE)

string.myc  <- runCIE(DGEs=dge,
                      methods="Quaternary",
                      databaseType="string",
                      databaseDir="path/to/ents/and/rels/files/")
write.table(string.myc, "documentation/stringMyc.tsv", quote=F,
            sep="\t", row.names=FALSE)

trrust.myc  <- runCIE(DGEs=dge,
                      methods="Fisher",
                      databaseType="TRRUST",
                      databaseDir="path/to/ents/and/rels/files/")

write.table(trrust.myc, "documentation/trrustMYC.tsv", quote=F,
            sep="\t", row.names=FALSE)

library(rcytoscapejs)

string  <- lapply(c("ents", "rels"), function(x) {
    read.table(paste0("path/to/ents/rels/string.", x),
               header=T, sep="\t")
})
names(string)  <- c("ents", "rels")

trrust  <- lapply(c("ents", "rels"), function(x) {
    read.table(paste0("path/to/ents/rels/TRRUST.", x),
               header=T, sep="\t")
})
names(trrust)  <- c("ents", "rels")

createCytoGraph(enrichment=tissChIP.myc, ents=tissChIP$ChIPall.ents,
                rels=tissChIP$ChIPall.rels, DGEs=dge)

createCytoGraph(enrichment=string.myc, ents=string$ents,
                rels=string$rels, DGEs=dge)

createCytoGraph(enrichment=trrust.myc, ents=trrust$ents,
                rels=trrust$rels, DGEs=dge)

pTcMYC  <- pathwayEnrichment(tissChIP.myc$name[1:20],)
pStMYC  <- pathwayEnrichment(string.myc$name[1:20],)
pTrMYC  <- pathwayEnrichment(trrust.myc$name[1:20],)

write.table(pTcMYC, "documentation/pathTcMYC.tsv", quote=F,
            sep="\t", row.names=FALSE)

write.table(pStMYC, "documentation/pathStMYC.tsv", quote=F,
            sep="\t", row.names=FALSE)

write.table(pTrMYC, "documentation/pathTrMYC.tsv", quote=F,
            sep="\t", row.names=FALSE)
