#' filterChIPAtlas
#'
#' @description Given selection parameters, retrieve a portion of the ChIP Atlas database in the
#' format of an .ents and a .rels file, which are useful for network and graph generation.  This
#' is also the format required by the runCIE and createCytoGraph functions in this library.
#' 
#' @usage readChIPAtlas(distance, cutoff, cutoffType, cellLines = NA, cellLineType=NA, cellLineDiagnosis = NA, outFileName = NA, writeToFile=TRUE, databaseDir=NA)
#'
#' @param distance A number, either 1, 5, or 10 which indicates what distance (in kb)
#' from the transcription start site should be considered.
#'
#' @param cutoff A number to be used in evaluating which targets have a sufficient binding score.
#' Binding score filtering depends on the next parameter, cutoffType
#'
#' @param cutoffType One of "average", "min", "max", or "auto".  If "average" is selected, only
#' targets which have an average binding score above the cutoff will be included.  If "min" or
#' "max" is selected, the binding scores from individual experiments will be evaluated to see
#' if they are above or below that threshold, respectively.  If none are, the target will not be
#' included.  If "auto" is selected, the 75th quantile of experiment-specific bindings scores for
#' each transcription factor will be determined.  Only targets with a bindings score above the
#' 75th quantile will be included.
#'
#' @param cellLines A vector of strings which are names of cell lines.  Only interactions which
#' have been documented in that cell line will be filtered by cutoff and then returned.
#' It is set to NA by default, and is an optional parameter.
#'
#' @param cellLineType A single string which corresponds to a name of an organ.  Only interactions
#' which have been documented in cell lines originating from that tissue will be filtered by
#' cutoff and then returned. It is set to NA by default, and is an optional parameter.
#'
#' @param cellLineDiagnosis A single string which corresponds to the diagnosis of the person the
#' cell line is from.  Only interactions which have been documented in cell lines originating from
#' a person with that diagnosis will be filtered by cutoff and then returned. It is set to NA by
#' default, and is an optional parameter.
#' 
#' @param outFileName A string of the desired output file name, optional.  It should be
#' the part of the file name before the extension (ex. "filename" becomes filename.rels).
#'
#' @param writeToFile A boolean value, determining whether the output will be returned to
#' the environment or written to file (FALSE and TRUE, respectively.  Defaults to TRUE
#'
#' @param databaseDir The path to the folder containing cellLines.rds and chip-atlas-*kb.rds.  The
#' function will search the current working directory if a database directory  is not provided.
#'
#' @param annotate Use annotations, defaults to false
#'
#' @param expectProgressObject optional parameter specifying if a Shiny progress
#' bar will be passed
#'
#' @param progressObject optional parameter providing the Shiny progress bar
#' 
#' @return Two dataframes, differentiated by the suffix to their names both in a returned object
#' and  when written to file. One, ending in .ents, is the entries from the database which passed
#' the filters specified.  The other, ending in .rels describes the relations between the entries.
#' If returned to the environment, the ents table is accessible by [outVarName]$filteredChIP.ents
#' and the rels table is accessible by [outVarName]$filteredChIP.rels
#' 
#' @export
#'
#' @examples
#' For 5kb distance from the TSS, an average score of 500, and the cell line SEM:
#' readChIPAtlas(5, 500, "average", "SEM")
#' For 5kb distance from the TSS, and automatic cutoff
#' readChIPAtlas(5, NA, "automatic")
#' readChIPAtlas(5, NA, "auto")
#' For 1kb distance from the TSS and select cell line by primary tissue type  blood
#' and automatic cutoff
#' readChIPAtlas(1, NA, "auto", cellLineType="blood")
#' Select by 1kb distance, autocutoff, and only cell lines from blood cells in people
#' diagnosed with leukemia 
#' ChIP1autoBldLeuk <- filterChIPAtlas(1, NA, "auto", NA, "blood", "leukemia")
#' ents <- ChIP1autoBldLeuk$filteredChIP.ents
#' rels <- ChIP1autoBldLeuk$filteredChIP.rels
#' 
filterChIPAtlas <- function(distance, cutoff, cutoffType, cellLines = NA,
                            cellLineType=NA, cellLineDiagnosis = NA,
                            outFileName = NA, writeToFile=FALSE, databaseDir=NA,
                            annotate=FALSE, expectProgressObject=FALSE,
                            progressObject=NA) {
    if(annotate) {
        relsFN  <- paste0("chip-atlas-", distance, "kb-anno.rels.rds")
    }
    else {
        relsFN  <- paste0("chip-atlas-", distance, "kb-collapse.rels.rds")
    }
    entsFN  <- paste0("chip-atlas-", distance, "kb-collapse.ents")
    if(is.na(databaseDir) & (!file.exists(entsFN) | file.exists(relsFN))) {
         stop(paste("Please provide a directory where the cellLines.rds file and",
                   " chip-atlas-*kb.rds files can be found or place the files in your",
                   " current working directory"))
    }
    if(!is.na(databaseDir)) {
        relsFN  <- paste0(databaseDir, relsFN)
        entsFN  <- paste0(databaseDir, entsFN)
    }
    if(!file.exists(relsFN) | !file.exists(entsFN)) {
        stop(paste("The directory you specified either does not exist or does not contain",
                   "the required files"))
    }
    rels  <- readRDS(relsFN)
    ents  <- read.table(entsFN, header=T, sep="\t", stringsAsFactors=F)
    if(!is.na(cellLineType[1])) {
        if(!is.na(cellLineDiagnosis)) {
            cellLinesTemp <- findCellLines(cellLineType, cellLineDiagnosis, databaseDir)
        }
        else {
            cellLinesTemp <- findCellLines(cellLineType, databaseDir=databaseDir)
        }
        if(!is.na(cellLines[1])) {
            cellLines <- c(cellLinesTemp, cellLines)
            cellLines <- unique(cellLines)
        }
        else {
            cellLines <- cellLinesTemp
        }
    }
    if(cutoffType == "average") {
        relsTemp  <- rels %>% dplyr::filter(cellLines=="Average", score > cutoff)
        rels  <- rels %>% dplyr::filter(srcuid %in% relsTemp$srcuid,
                                        trguid %in% relsTemp$trguid,
                                        cellLines != "Average")
    }
    else if(cutoffType == "min") {
        rels  <- rels %>% dplyr::filter(cellLines!="Average", score > cutoff)
    }
    else if(cutoffType == "max") {
        rels  <- rels %>% dplyr::filter(cellLines!="Average", score < cutoff)
    }
    else if(cutoffType == "automatic" || cutoffType == "auto") {
        prot <- unique(rels$srcuid)
        relsTemp  <- rels %>% dplyr::filter(cellLines != "Average") %>%
            group_by(srcuid)
        cutoffs  <- unlist(relsTemp %>%
                           summarize(cutoff = quantile(score)[4]) %>%
                           dplyr::select(cutoff))
        relsOut  <- rels %>% dplyr::filter(cellLines != "Average")
        if(!expectProgressObject) {
            print(paste0("Fitlering rels table by fourth quantile of scores",
                         " for each TF"))
            rels  <- do.call(rbind, pblapply(1:length(prot), function(x) {
                relsOut %>% dplyr::filter(score > cutoffs[x], srcuid==prot[x])    
            }))
        }
        else {
            rels  <- do.call(rbind, pblapply(1:length(prot), function(x) {
                relsOut %>% dplyr::filter(score > cutoffs[x], srcuid==prot[x])    
                progressObject$set(message="Filtering entries",
                             value=((x/length(prot))/2))
            }))
        }    
    }
    else {
        print("Please provide a valid cutoff type")
    }
    if(!is.na(cellLines[1])) {
        cellLines.temp  <- cellLines
        rels <- rels %>% dplyr::filter(cellLines %in% cellLines.temp)
    }
    print("Processing to output format")
    ents.prot  <- ents %>% dplyr::filter(type=="Protein") %>%
        mutate(id = NULL, type=NULL, srcName= name,
               name=NULL)
    ents.mRNA  <- ents %>% dplyr::filter(type=="mRNA") %>%
        mutate(id = NULL, type=NULL, trgName =name,
               name =NULL)
    rels  <- rels %>% left_join(ents.prot, by=c("srcuid"="uid"))
    rels  <- rels %>% left_join(ents.mRNA, by=c("trguid"="uid"))
    entsNew.mRNA  <- data.frame(name=unique(rels$trgName),
                                stringsAsFactors=FALSE)
    entsNew.prot  <- data.frame(name=unique(rels$srcName),
                                stringsAsFactors=FALSE)
    entsToJoinP  <- ents %>% dplyr::filter(type=="Protein") %>% mutate(uid=NULL)
    entsToJoinM  <- ents %>% dplyr::filter(type=="mRNA") %>% mutate(uid=NULL)
    entsNew.prot  <-  entsNew.prot %>%
        left_join(entsToJoinP, by=c("name"="name"))
    entsNew.mRNA  <- entsNew.mRNA %>%
        left_join(entsToJoinM, by=c("name"="name"))
    ents  <- rbind(entsNew.prot, entsNew.mRNA)
    ents  <- cbind(uid=1:nrow(ents), ents)
    rels  <-  rels %>% dplyr::mutate(trgName=NULL, srcName=NULL)
    prot  <- unique(rels$srcuid)
    if(expectProgressObject) {
        if(cutoffType == "automatic" || cutoffType == "auto") {
            rels  <- do.call(rbind, lapply(1:length(prot), function(x) {
                targs  <- rels %>% dplyr::filter(srcuid==prot[x]) %>% group_by(trguid)
                type  <- targs %>%
                    summarize(type = ifelse((sum(type=="increase")==(0.75*length(type))),
                                        "increase",
                                     ifelse(sum(type=="decrease")==(0.75*length(type)),
                                            "decrease", "conflict")),
                              cellLinesTotal  = paste(cellLines, collapse=", "),
                              scoreTotal  = paste(score, collapse=", "))
                progressObject$set(message="Processing databse for enrichment",
                                   value=((x/length(prot))/2))
                cbind(scruid=x, as.data.frame(type))
            }))
        }
        else {
            rels  <- do.call(rbind, lapply(1:length(prot), function(x) {
                targs  <- rels %>% dplyr::filter(srcuid==prot[x]) %>% group_by(trguid)
                type  <- targs %>%
                    summarize(type = ifelse((sum(type=="increase")==(0.75*length(type))),
                                        "increase",
                                     ifelse(sum(type=="decrease")==(0.75*length(type)),
                                            "decrease", "conflict")),
                              cellLinesTotal  = paste(cellLines, collapse=", "),
                              scoreTotal  = paste(score, collapse=", "))
                progressObject$set(message="Processing databse for enrichment",
                                   value=(x/length(prot)))
                cbind(scruid=x, as.data.frame(type))
            }))
        }
    }
    else {
        rels  <- do.call(rbind, pblapply(prot, function(x) {
            targs  <- rels %>% dplyr::filter(srcuid==x) %>% group_by(trguid)
            type  <- targs %>%
                summarize(type = ifelse((sum(type=="increase")==(0.5*length(type))),
                                        "increase",
                                 ifelse(sum(type=="decrease")==(0.5*length(type)),
                                        "decrease", "conflict")),
                          cellLinesTotal  = paste(cellLines, collapse=", "),
                          scoreTotal  = paste(score, collapse=", "))
            cbind(scruid=x, as.data.frame(type))
        }))
    }
    rels  <- cbind(uid=1:nrow(rels), rels)
    if(writeToFile) {
        if(is.na(outFileName)) {
            write.table(ents, "ChIPfiltered.ents", row.names=FALSE, sep="\t",
                        quote=F)
            write.table(rels, "ChIPfiltered.rels", row.names=FALSE, sep="\t",
                        quote=F)
        }
        else {
            write.table(ents, paste0(outFileName, ".ents"), row.names=FALSE,
                        sep="\t", quote=F)
            write.table(rels, paste0(outFileName, ".rels"), row.names=FALSE,
                        sep="\t", quote=F)
        }
    }
    else {
        list("ents"=ents, "rels"=rels)
    }
    
}

findCellLines <- function(cellLinePTtype, cellLineDiagnosis = NA, databaseDir) {
    rdsFN <- "cellLines.rds"
    if(is.na(databaseDir) & !file.exists(rdsFN)) {
        stop(paste("Please provide a directory where the cellLines.rds file and",
                   " chip-atlas-*kb.rds files can be found or place the files in your",
                   " current working directory"))
    }
    else if(!is.na(databaseDir)) {
        rdsFN <- paste(databaseDir, rdsFN, sep = "")
    }
    if(!file.exists(rdsFN)) {
        stop(paste("The directory you specified either does not exist or does not contain",
                   "the required files"))
    }
    else {
        cellTypeDatabase <- readRDS(rdsFN)
    }
    if(!is.na(cellLineDiagnosis)) {
        cellLines <- cellTypeDatabase  %>%
            dplyr::filter(grepl(cellLinePTtype, Primary.Tissue, ignore.case=T),
                          grepl(cellLineDiagnosis, Tissue.Diagnosis, ignore.case=T)) %>%
            dplyr::select(Cell.Line.Name)
        gc()
        cellLines <- unlist(cellLines)
        as.character(cellLines)
    }
    else {
        cellLines <- cellTypeDatabase %>%
            dplyr::filter(grepl(cellLinePTtype, Primary.Tissue, ignore.case=T)) %>%
            dplyr::select(Cell.Line.Name)
        gc()
        cellLines <- unlist(cellLines)
        as.character(cellLines)
    }
}
