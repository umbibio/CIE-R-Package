require(dplyr)
require(readr)
require(org.Hs.eg.db)
require(xlsx)
## This needs to be improved
## source("https://bioconductor.org/biocLite.R")
## biocLite("GenomicRanges")


#' filterChIPAtlas
#'
#' @description Given selection parameters, retrieve a portion of the ChIP Atlas database in the
#' format of an .ents and a .rels file, which are useful for network and graph generation.
#'
#' @usage readChIPAtlas(distance, cutoff, cutoffType, cellLines = NA)
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
#' have been documented in that cell line will be filtered by cutoff and then written to file.
#' It is set to NA by default, and is an optional parameter.
#'
#' @param cellLineType A single string which corresponds to the organ you wish to see cell lines
#' from
#'
#' @param cellLineDiagnosis A single string which corresponds to the diagnosis of the person the
#' cell line is from
#' @param outFileName A string of the desired output file name, optional.  It should be
#' the part of the file name before the extension (ex. "filename" becomes filename.rels)
#'
#' @param writeToFile A boolean value, determining whether the output will be returned to
#' the environment or written to file (FALSE and TRUE, respectively.  Default to TRUE
#' 
#' @return Writes two files to disk, one with the file type .ents, which is the recorded entries
#' from the database that passed the filters specified, and one with the file type .rels, which
#' describes the relations between them.
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
#' filterChIPAtlas(1, NA, "auto", NA, "blood", "leukemia")
#' 
filterChIPAtlas <- function(distance, cutoff, cutoffType, cellLines = NA,
                            cellLineType=NA, cellLineDiagnosis = NA,
                            outFileName = NA, writeToFile=TRUE) {
    rdsFN <- paste("chip-atlas-", distance, "kb.rds", sep = "")
    if(!file.exists(rdsFN)) {
        folderName <- paste("chip-atlas-", distance, "kb/", sep = "")
        tsvFiles <- list.files(folderName)
        ChIPlist <- list(1:length(tsvFiles))
        index <- 1;
        for(f in tsvFiles) {
            ChIPlist[[index]] <- read_tsv(paste(folderName, f, sep =""))
            names(ChIPlist)[index] <- strsplit(f, "\\.")[[1]][1]
            ChIPlist[[index]] <- cbind(names(ChIPlist)[index], ChIPlist[[index]])
            colnames(ChIPlist[[index]])[1] <- "TF"
            index <- index + 1
        }
        saveRDS(ChIPlist, rdsFN)
    }
    else {
        ChIPlist <- readRDS(rdsFN)
    }
    if(!is.na(cellLineType)) {
        if(!is.na(cellLineDiagnosis)) {
            cellLinesTemp <- findCellLines(cellLineType, cellLineDiagnosis)
        }
        else {
            cellLinesTemp <- findCellLines(cellLineType)
        }
        if(!is.na(cellLines)) {
            cellLines <- c(cellLinesTemp, cellLines)
            cellLines <- unique(cellLines)
        }
        else {
            cellLines <- cellLinesTemp
        }
    }
    if(!is.na(cellLines)) {
        ChIPlist <- filterByCellLine(ChIPlist, cellLines)
    }
    if(writeToFile==TRUE) {
        if(!is.na(outFileName)) {
            if(cutoffType == "average") {
                getByAverageBS(ChIPlist, cutoff, outFileName)
            }
            else if(cutoffType == "min") {
                getByMinBS(ChIPlist, cutoff, outFileName)
            }
            else if(cutoffType == "max") {
                getByMaxBS(ChIPlist, cutoff, outFileName)
            }
            else if(cutoffType == "automatic" || cutoffType == "auto") {
                getByAutoBS(ChIPlist, outFileName)
            }
            else {
                print("Please provide a valid cutoff type")
            }
        }
        else {
            if(cutoffType == "average") {
                getByAverageBS(ChIPlist, cutoff)
            }
            else if(cutoffType == "min") {
                getByMinBS(ChIPlist, cutoff)
            }
            else if(cutoffType == "max") {
                getByMaxBS(ChIPlist, cutoff)
            }
            else if(cutoffType == "automatic" || cutoffType == "auto") {
                getByAutoBS(ChIPlist)
            }
            else {
                print("Please provide a valid cutoff type")
            }
        }
    }
    else {
        if(cutoffType == "average") {
            getByAverageBS(ChIPlist, cutoff, NA, FALSE)
        }
        else if(cutoffType == "min") {
            getByMinBS(ChIPlist, cutoff, NA, FALSE)
        }
        else if(cutoffType == "max") {
            getByMaxBS(ChIPlist, cutoff, NA, FALSE)
        }
        else if(cutoffType == "automatic" || cutoffType == "auto") {
            getByAutoBS(ChIPlist, NA, FALSE)
        }
        else {
            print("Please provide a valid cutoff type")
        }   
    }
}
## Helper function which writes the .rels and .ents files
writeRelsEnts <- function(ChIPdata, outFileName, writeToFile) {
    mRNAs <- unique(ChIPdata$Target_genes)
    TFs <- unique(ChIPdata$TF)
    
    ChIP.ents.TF <- data.frame(name = TFs, type = Rle('Protein', length(TFs)),
                               stringsAsFactors =F)
    ChIP.ents.mRNA <- data.frame(name = mRNAs, type = Rle('mRNA', length(mRNAs)),
                                 stringsAsFactors = F)
    ChIP.ents.TF$name <- as.character(ChIP.ents.TF$name)
    
    XX <- select(org.Hs.eg.db, ChIP.ents.TF$name,
                 columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
    YY <- select(org.Hs.eg.db, keys=as.character(ChIP.ents.mRNA$name),
                 columns=c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")

    YY <- YY %>% distinct(SYMBOL, .keep_all = T)
    ChIP.ents.TF <- left_join(ChIP.ents.TF, XX, by = c('name' = 'SYMBOL'))
    ChIP.ents.mRNA <- left_join(ChIP.ents.mRNA, YY, by = c('name' = 'SYMBOL'))

    ChIP.ents <- rbind(ChIP.ents.TF, ChIP.ents.mRNA)
    ChIP.ents <- ChIP.ents %>% mutate(uid = 1:nrow(ChIP.ents)) %>%
        transmute(uid = uid, name = name, id = ENTREZID, type = type)

    ChIP.ents.TF <- ChIP.ents %>% filter(type == 'Protein')
    ChIP.ents.mRNA <- ChIP.ents %>% filter(type == 'mRNA')

    ChIP.rels <- ChIPdata %>% na.omit() %>% group_by(TF) %>%
        mutate(srcuid = ChIP.ents.TF$uid[match(TF, ChIP.ents.TF$name)]) %>%
        mutate(trguid = ChIP.ents.mRNA$uid[match(Target_genes, ChIP.ents.mRNA$name)]) %>%
        ungroup()
    ChIP.rels <- ChIP.rels %>%
        transmute(uid = 1:nrow(ChIP.rels), srcuid = srcuid, trguid = trguid,
                  type = 'conflict', pmids = NA, nls = NA)
    if(writeToFile) {
        if(!is.na(outFileName)) {
            write.table(ChIP.ents, paste(outFileName, ".ents", sep=""),
                        col.names=T,row.names = F, sep = '\t', quote = F)
            write.table(ChIP.rels, paste(outFileName, ".rels", sep=""),
                        col.names=T, row.names = F, sep = '\t', quote = F)
        }
        write.table(ChIP.ents, 'ChIPfilter.ents', col.names=T,
                    row.names = F, sep = '\t', quote = F)
        write.table(ChIP.rels, 'ChIPfilter.rels', col.names=T,
                    row.names = F, sep = '\t', quote = F)
    }
    else {
        returnData <- list(ChIP.ents, ChIP.rels)
        names(returnData) <- c("filteredChIP.ents", "filteredChIP.rels")
        return(returnData)
    }
}

## Helper function which allows a vectorized call on the list of objects for selection by
## average binding score
getByAverageBShelper <- function(ChIPlistItem, cutoff) {
    returnData <- data.frame(Target_genes = NA, TF = NA, stringsAsFactors = F)
    tf <- ChIPlistItem[1,1]
    rawData <- ChIPlistItem[,2:ncol(ChIPlistItem)]
    returnData <- rawData %>% filter(rawData[,2] > cutoff) %>%
        dplyr::select(Target_genes) %>%
        mutate(TF = tf)
    returnData <- returnData[complete.cases(returnData),]
    returnData
}

## Function to get by average binding score
getByAverageBS <- function(ChIPlist, cutoff, outFileName = NA, writeToFile=TRUE) {
    ChIP <- do.call(rbind, lapply(ChIPlist, function(x) {getByAverageBShelper(x, cutoff)}))
    if(writeToFile) {
        writeRelsEnts(ChIP, outFileName, writeToFile)
    }
    else {
        result <- writeRelsEnts(ChIP, outFileName, writeToFile)
        return (result)
    }
}


getByMinBShelper <- function(ChIPlistItem, cutoff) {
    returnData <- data.frame(Target_genes = NA, TF = NA, stringsAsFactors = F)
    tf <- ChIPlistItem[1,1]
    rawData <- ChIPlistItem[,2:ncol(ChIPlistItem)]
    if(ncol(rawData) == 3) {
        data <- rawData %>% filter(.[[3]] > cutoff)
    }
    else {
        data <- rawData[apply(rawData[,3:ncol(rawData)], 1,
                              function(x) {any(x > cutoff)}),]
    }
    returnData <- data %>%
        dplyr::select(Target_genes) %>%
        mutate(TF = tf)
    returnData <- returnData[complete.cases(returnData),]
    returnData

}

## Function to get by average binding score,
getByMinBS <- function(ChIPlist, cutoff, outFileName = NA, writeToFile=TRUE) {
    ## Data frame from repository with specified average
    ChIP <- do.call(rbind, lapply(ChIPlist, function(x) {getByMinBShelper(x, cutoff)}))

    if(writeToFile) {
        writeRelsEnts(ChIP, outFileName, writeToFile)
    }
    else {
        result <- writeRelsEnts(ChIP, outFileName, writeToFile)
        return (result)
    }
}

## Helper function to get by max binding score
getByMaxBShelper <- function(ChIPlistItem, cutoff) {
    returnData <- data.frame(Target_genes = NA, TF = NA, stringsAsFactors = F)
    tf <- ChIPlistItem[1,1]
    rawData <- ChIPlistItem[,2:ncol(ChIPlistItem)]
    if(ncol(rawData) == 3) {
        data <- rawData %>% filter(.[[3]] < cutoff)
    }
    else {
        data <- rawData[apply(rawData[,3:ncol(rawData)], 1,
                              function(x) {any(x < cutoff)}),]
    }
    returnData <- data %>%
        dplyr::select(Target_genes) %>%
        mutate(TF = tf)
    returnData <- returnData[complete.cases(returnData),]
    returnData
}

## Function to get by average binding score,
getByMaxBS <- function(ChIPlist, cutoff, outFileName = NA, writeToFile=TRUE) {
    ChIP <- do.call(rbind, lapply(ChIPlist, function(x) {getByMinBShelper(x, cutoff)}))
    if(writeToFile) {
        writeRelsEnts(ChIP, outFileName, writeToFile)
    }
    else {
        result <- writeRelsEnts(ChIP, outFileName, writeToFile)
        return (result)
    }
}

## Helper function to get by automatic threshold
getByAutoBShelper <- function(ChIPlistItem) {
    returnData <- data.frame(Target_genes = NA, TF = NA, stringsAsFactors = F)
    tf <- ChIPlistItem[1,1]
    rawData <- ChIPlistItem[,2:ncol(ChIPlistItem)]
    allScores <- unlist(rawData[,3:ncol(rawData)])
    cutoff <- quantile(allScores)[4]
    if(ncol(rawData) == 3) {
        data <- rawData %>% filter(.[[3]] > cutoff)
    }
    else {
        data <- rawData[apply(rawData[,3:ncol(rawData)], 1,
                              function(x) {any(x > cutoff)}),]
    }
    returnData <- data %>%
        dplyr::select(Target_genes) %>%
        mutate(TF = tf)
    returnData <- returnData[complete.cases(returnData),]
    returnData
}

## Function to get by an automatic threshold, calculated for each transcription factor
getByAutoBS <- function(ChIPlist, outFileName = NA, writeToFile=TRUE) {
    ChIP <- do.call(rbind, lapply(ChIPlist, getByAutoBShelper))
    if(writeToFile) {
        writeRelsEnts(ChIP, outFileName, writeToFile)
    }
    else {
        result <- writeRelsEnts(ChIP, outFileName, writeToFile)
        return(result)
    }
}


## Function which helps with selecting by cell line
filterByCellLineHelper <- function(ChIPlistItem, cellLines) {
    if(length(cellLines) > 1) {
        colIndices <- sapply(cellLines,
                             function(x) {grep(x, ignore.case = T, colnames(ChIPlistItem))})
        colIndices <- unlist(colIndices)
        if(length(colIndices) == 0) {
            NA
        }
        else {
            colIndices <- c(1, 2, 3, colIndices)
            returnData <- ChIPlistItem[,colIndices]
            returnData
        }
    }
    else {
        colIndex <- grep(cellLines, ignore.case = T, colnames(ChIPlistItem))
        if(length(colIndex) == 0) {
            NA
        }
        else {
            colIndex <- c(1, 2, 3, colIndex)
            returnData <- ChIPlistItem[,colIndex]
            returnData
        }
    }
}

## Function which returns a modified list of data frames with entries which do not have the
## cell line removed
filterByCellLine <- function(ChIPlist, cellLines) {
    returnList <- lapply(ChIPlist, function(x) {
        filterByCellLineHelper(x, cellLines) })
    returnList <- returnList[!is.na(returnList)]
    returnList
}


findCellLines <- function(cellLinePTtype, cellLineDiagnosis = NA) {
    rdsFN <- "cellLines.rds"
    if(!file.exists(rdsFN)) {
        cellTypeDatabase <- read.xlsx("cellLines.xlsx", sheetIndex=1, stringsAsFactors=F)
        saveRDS(cellTypeDatabase, rdsFN)
        gc()
    }
    else {
        cellTypeDatabase <- readRDS("cellLines.rds")
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
