# Design
# library(Bc)
# 
# (Reads) -> bc_runQc -> bc_filterSeq -> (Reads) -> bc_runQc -> plot
# (Reads) -> bc_extract -> bc_cure -> check -> bc_curec -> count -> biological analysis
# 
# runQc() -> class BarcodeObj(meta, rawBc)
# filterFq.Barcode() -> class BarcodeObj(meta, rawBc)
# 
# BarcodeObj -> List(meta, messyBc, cleanBc, cleanProc)
# 
# extractBc.Barcode() -> class BarcodeObj(meta, rawBc, messyBc)
# 
# cureBc.Barcode() -> class BarcodeObj(meta, rawBc, messyBc, cleanBc)
# 
# countBc.Barcode() -> class BarcodeObj(meta, rawBc, messyBc, cleanBc, countBc)
# 
# 
# barcodeObj[, 1]
# barcodeObj[, "sample_name1"]
# barcodeObj[, 1:2]
# barcodeObj[, c("sample_name1", "sample_name2")]
# barcodeObj["AACT", ]
# renames(barcodeObj, c("old_sample_name", "new_sample_name"))
# subset(barcodeObj, name = c("sample_name1", "sample_name2"))
# subset(barcodeObj, barcode_white_list = c("AAAA", "TTT"), bacode_black = c("GGG", "CCC"))
# barcodeObj1 + barcodeObj2
# barcodeObj - black_list
# barcodeObj * white_list
# samplenames(bc_obj)
# bc2genBarCode(bc_obj)


# TODO:
# One function do all the job
# filterFq() -> output the filtering detail
# make function to cross talk with genBaRcode and other packages
# create function to identify the cutoff point
# evaluate the sequencing depth (saturation of the sequene)


##' foo: A package for computating the notorious bar statistic
##'
##' The foo package provides three categories of important functions:
##' foo, bar and baz.
##' 
##' @section Foo functions:
##' The foo functions ...
##'
##' @docType package
##' @name foo

#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist
#' @importFrom ShortRead readFastq
#' @importFrom Biostrings readDNAStringSet
#' @import ggplot2
#' @import Biostrings
#' @import ShortRead
NULL

. <- list()

