# A script describing the steps involved in making the data object(s). Output of
# this script should be files on disk ready to be pushed to S3.

## This script is meant to be run in an R session with a working directory
## that is set to the inst/extdata/scripts directory of this package.
##
## The user needs to download the *.xml version of the latest MSigDB gene set
## definitions, which we parse and convert into a GeneSetDb object.
##
## This file is downloaded from:
##   http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.2/msigdb_v5.2.xml
## Set the `fn` argument to the path that the xml file was saved to
library(multiGSEA)
library(data.table)
library(XML)
library(GSEABase)
library(org.Hs.eg.db)

## Downoad the latest msigdb file:
## http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?\
## filePath=/resources/msigdb/6.1/msigdb_v6.1.xml
fn <- '~/Downloads/msigdb_v6.1.xml'

## Some utility functions that were extracted from GSEABase. ==================
##
## I am modifying it here, because some 'factories' are not recognized
## (for category == 'archived' for instance) and I don't want the whole thing
## to fail
myXMLNodeToGeneSetFactory <- function(file, membersId = "MEMBERS_SYMBOLIZED") {
  .mkSplit <- function(x) {
    if (is.null(x))
      character(0)
    else unlist(strsplit(x, GSEABase:::.BROAD_SEPARATOR))
  }
  url <- NULL
  if (length(file) == 1) {
    isUri <- grep("^(http|ftp|file)://", file)
    if (length(isUri) == 1 && isUri == 1)
      url <- file
    else if (file.exists(file)) {
      full <- path.expand(file)
      if (length(grep("^/", full)) == 1)
        url <- paste("file:/", full, sep = "")
      else if (.Platform$OS.type == "windows")
        url <- paste("file:///", full, sep = "")
      else url = full
    }
  }
  symbolId <- SymbolIdentifier()
  function(node) {
    attrs <- as.list(xmlAttrs(node))
    args <- list(symbolId, setName = attrs[["STANDARD_NAME"]],
                 setIdentifier = attrs[["SYSTEMATIC_NAME"]], geneIds = unique(.mkSplit(attrs[[membersId]])),
                 organism = attrs[["ORGANISM"]], urls = c(getBroadSet = url,
                                                          .mkSplit(attrs[["EXTERNAL_DETAILS_URL"]])), collectionType = {
                                                            categories <- .mkSplit(attrs[["CATEGORY_CODE"]])
                                                            subcategories <- .mkSplit(attrs[["SUB_CATEGORY_CODE"]])
                                                            category <- subcategory <- as.character(NA)
                                                            if (length(categories) >= 1) category <- tolower(categories[[1]])
                                                            if (length(subcategories) >= 1) subcategory <- subcategories[[1]] else if (length(categories) >=
                                                                                                                                       2) subcategory <- categories[[2]]
                                                            if (length(categories) > 2 || (length(categories) >
                                                                                           1 && length(subcategories) != 0)) {
                                                              fmt <- "Broad 'CATEGORY_CODE' too long: '%s'"
                                                              txt <- paste(categories, collapse = "' '")
                                                              warning(sprintf(fmt, txt))
                                                            }
                                                            MyBroadCollection(category = mkScalar(category),
                                                                              subCategory = mkScalar(subcategory))
                                                          }, contributor = attrs[["CONTRIBUTOR"]], pubMedIds = attrs[["PMID"]],
                 shortDescription = attrs[["DESCRIPTION_BRIEF"]],
                 longDescription = attrs[["DESCRIPTION_FULL"]], TAGS = NULL,
                 MESH = NULL, CHIP = NULL, MEMBERS = NULL)
    args <- args[!sapply(args, is.null)]
    do.call(GeneSet, args)
  }
}

MyBroadCollection <- function(category = "c1", subCategory = NA, ...)  {
  # if (length(category) != 1 || !(category %in% c("c1", "c2",
  #                                                "c3", "c4", "c5", "c6", "c7", "h")))
  #   stop(sprintf("invalid BroadCollection category: '%s'",
  #                paste(category, collapse = "', '")))
  new("BroadCollection", category = mkScalar(category), subCategory = mkScalar(as.character(subCategory)))
}

.fromXML <- function(file, node, handler, ...) {
  res <- xmlTreeParse(file, useInternalNodes = TRUE, ...)
  geneSets <- getNodeSet(res, node, fun = handler)
  free(res)
  geneSets
}

getMsigSets <- function (uri, ..., membersId = c("MEMBERS_SYMBOLIZED", "MEMBERS_EZID"))
{
  membersId <- match.arg(membersId)
  factories <- sapply(uri, myXMLNodeToGeneSetFactory, membersId = membersId)

  # tryCatch({
  #   geneSets <- unlist(mapply(GSEABase:::.fromXML, uri, "//GENESET",
  #                             factories, SIMPLIFY = FALSE, USE.NAMES = FALSE))
  # }, error = function(err) {
  #   msg <- paste("'getMsigSets' failed to create gene sets:\n  ",
  #                conditionMessage(err))
  #   msg
  # })
  # # GeneSetCollection(geneSets)
  # geneSets
  unlist(mapply(GSEABase:::.fromXML, uri, "//GENESET",
                factories, SIMPLIFY = FALSE, USE.NAMES = FALSE))

}

## Processing =================================================================

## Process the downloaded file and save into intermediary format just incase
gs.list <- getMsigSets(fn, membersId="MEMBERS_EZID")
saveRDS(gs.list ,'~/Downloads/GeneSet-list-v61-MEMBERS_EZID.rds')

## Convert into GSEABase::GeneSetCollection
gsc <- GeneSetCollection(gs.list)

if (FALSE) {
  ## sanity check
  g.go <- gsc[[1]]
  g.h <- gsc[['HALLMARK_ANGIOGENESIS']]
  sapply(slotNames(g.h), function(x) slot(g.h, x), simplify=FALSE)
}

info <- lapply(gsc, function(gs) {
  data.table(collection=bcCategory(collectionType(gs)),
             name=setName(gs),
             organism=organism(gs),
             subcategory=gs@collectionType@subCategory,
             featureId=geneIds(gs))
})
info.all <- rbindlist(info)
for (col in names(info.all)) {
  info.all[, (col) := as.character(info.all[[col]])]
}

info <- subset(info.all, collection %in% c('h', paste0('c', 1:7)))
setdiff(info.all$collection, info$collection) ## "archived"

info[, symbol := mapIds(org.Hs.eg.db, featureId, 'SYMBOL', 'ENTREZID')]

gdb <- GeneSetDb(info[, list(collection, name, featureId, symbol)])

## Update the gdb@table
gdb@table <- local({
  gs <- unique(info[, list(collection, name, subcategory, organism)])
  stopifnot(sum(duplicated(gs$name)) == 0)
  gst <- merge(gdb@table, gs, by=c('collection', 'name'))
  setkeyv(gst, key(gdb@table))
  stopifnot(all.equal(gst[, names(gdb@table), with=FALSE], gdb@table))
  gst
})

## Beef up collectionMetadata --------------------------------------------------
## URL function
url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(gdb)$collection)) {
  geneSetCollectionURLfunction(gdb, col) <- url.fn
  featureIdType(gdb, col) <- EntrezIdentifier()
  gdb <- addCollectionMetadata(gdb, col, 'source', 'MSigDB_v6.1')
}

org(gdb) <- 'Homo_sapiens'

# gdb.fn <- sprintf('MSigDB.Homo_sapiens.GeneSetDb.rds', species)
gdb.fn <- '../extdata/GeneSetDb.MSigDB.Hsapiens-entrez.v61.rds'
saveRDS(gdb, gdb.fn)

# Create Ensembl version -------------------------------------------------------
library(multiGSEA)
library(biomaRt)
library(dplyr)
library(GSEABase)
library(dtplyr)

hgdb <- gdb
# hgdb <- readRDS("inst/extdata/GeneSetDb.MSigDB.Hsapiens-entrez.v61.rds")
hdf <- hgdb@db
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
exref <- getBM(
  attributes = c("entrezgene", "ensembl_gene_id", "hgnc_symbol"),
  filters = "entrezgene",
  values = unique(hdf$featureId),
  mart = mart) %>%
  transmute(entrezgene = as.character(entrezgene),
            featureId = ensembl_gene_id, symbol = hgnc_symbol)

hdf.ens <- hdf %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  dplyr::select(collection, name, entrezgene = featureId) %>%
  dplyr::inner_join(exref, by = "entrezgene") %>%
  dplyr::arrange(collection, name, symbol) %>%
  dplyr::distinct(collection, name, featureId, .keep_all = TRUE) %>%
  dplyr::select(-entrezgene)

# Use addGeneSetMetadata function instead of this stuff in the future
gdb2 <- GeneSetDb(hdf.ens)
take.cols <- c('collection', 'name', setdiff(colnames(hgdb@table), colnames(gdb2@table)))
meta <- hgdb@table[gdb2@table, take.cols, with=FALSE]
mnew <- gdb2@table[meta]
stopifnot(
  all.equal(gdb2@table[, list(collection, name)], mnew[, list(collection, name)]),
  all.equal(gdb2@table$N, mnew$N))
gdb2@table <- mnew

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(gdb2)$collection)) {
  geneSetCollectionURLfunction(gdb2, col) <- url.fn
  featureIdType(gdb2, col) <- ENSEMBLIdentifier()
  gdb2 <- addCollectionMetadata(gdb2, col, 'source', 'MSigDB_v6.1')
}

org(gdb2) <- 'Homo_sapiens'
gdb.fn <- '../extdata/GeneSetDb.MSigDB.Hsapiens-ensembl.v61.rds'
saveRDS(gdb2, gdb.fn)
