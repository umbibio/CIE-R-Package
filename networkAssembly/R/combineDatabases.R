require(dplyr)

#' combineDatabases
#' @description Given either a vector of file paths or a list of data frames, combine the
#' ents and rels tables into two master tables.
#' @usage combineDatabases(databaseEnts, databaseRels)
#'
#' @param databaseEnts Either a character vector of file paths or a list of data frames
#' containing the ents tables
#'
#' @param databaseRels Eithre a character vector of file paths or a list of data frames
#' containing the rels tables  Note: the two paramters must be the same length and
#' containing matching databases.
#'
#' @return A list containing the master rels and ents tables.
#'
#' @export
#' @examples
#' Read from file
#' combineDatabases(c("TRED.ents", "string.ents"), c("TRED.rels", "string.rels"))
#' Use list of data frames
#' ChIP1ap <- filterChIPAtlas(1, NA, "auto", NA, "prostate", NA, NA, FALSE)
#' ChIP1a500 <- filterChIPAtlas(1, 500, "average", NA, NA, NA, NA, FALSE)
#'
#' ents <- list(ChIP1ap$filteredChIP.ents, ChIP1a500$filteredChIP.ents)
#' rels <- list(ChIP1ap$filteredChIP.rels, ChIP1a500$filteredChIP.rels)
#' combineDatabases(ents, rels)


combineDatabases <- function(databaseEnts, databaseRels) {
    if(class(databaseEnts) == "character" & class(databaseRels) == "character") {
        databaseEntsList <- lapply(databaseEnts, function(x) {
            read.table(x, sep="\t", header = T, stringsAsFactors=F) } )
        databaseRelsList <- lapply(databaseRels, function(x) {
            read.table(x, sep="\t", header = T, stringsAsFactors=F) } )
    } else if(class(databaseEnts) == "list" && class(databaseRels) == "list" &&
         class(databaseEnts[[1]]) == "data.frame" &&
         class(databaseRels[[1]]) == "data.frame") {
        databaseEntsList <- databaseEnts
        databaseRelsList <- databaseRels
    } else {
        stop("This method only supports character vectors which are lists of file paths to rels and ents files or a list of data frames consisting of rels and ents tables.  The input you provided is not supported")
    }
    if(length(databaseEntsList) != length(databaseRelsList)) {
        stop("Please make sure that the your lists of databases match")
    }

    databaseEntsList <- lapply(databaseEntsList, function(x) {
        x$type <- as.character(x$type); x } )
    TFs <- lapply(databaseEntsList, function(x) {
        x %>% dplyr::filter(grepl("protein", type, ignore.case=T)) %>%
            dplyr::select(name) } )

    mRNAs <- lapply(databaseEntsList, function(x) {
        x %>% dplyr::filter(grepl("mrna", type, ignore.case=T))%>%
            dplyr::select(name) } )
                        
    for(i in 1:length(databaseRelsList)) {
        databaseRelsList[[i]]$srcuid <- databaseEntsList[[i]]$name[databaseRelsList[[i]]$srcuid]
        databaseRelsList[[i]]$trguid <- databaseEntsList[[i]]$name[databaseRelsList[[i]]$trguid]
    }
    
    masterEnts <- do.call(rbind, lapply(databaseEntsList, function(x) {x}))
    masterRels <- do.call(rbind, lapply(databaseRelsList, function(x) {x}))

    masterEnts <- unique(masterEnts)
    masterRels <- unique(masterRels)

    masterEnts <- masterEnts %>% dplyr::arrange(desc(type)) %>%
        dplyr::mutate(uid= 1:nrow(masterEnts))

    masterRels <-  masterRels %>%
        dplyr::mutate(srcuid = masterEnts$uid[match(srcuid, masterEnts$name)]) %>%
        dplyr::mutate(trguid = masterEnts$uid[match(trguid, masterEnts$name)]) %>%
        dplyr::arrange(srcuid) %>%
        dplyr::mutate(uid = 1:nrow(masterRels))

    result <- list(masterEnts, masterRels)
    names(result) <- c("combinedEnts", "combinedRels")
    result
}
    
