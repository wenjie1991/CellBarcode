# TODO:
# One function do all the job
# filterFq() -> output the filtering detail
# make function to cross talk with genBaRcode and other packages
# create function to identify the cutoff point
# evaluate the sequencing depth (saturation of the sequene)

#' DNA Barcode Analysis toolkit
#'
#' This package performs DNA Barcode (genetic lineage tracing) analysis. The
#' package can handle all kinds of DNA barcodes, as long as the barcode within a
#' single sequencing read and has a pattern which can be matched by a regular
#' expression. This package can handle barcode with flexible length, with or
#' without UMI (unique molecular identifier). This tool also can be used for
#' pre-processing of some amplicon data such as CRISPR gRNA screening, immune
#' repertoire sequencing and meta genome data.
#'
#' @name CellBarcode
#' @docType package
#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist .N :=
#' @importFrom Biostrings readDNAStringSet
#' @importFrom methods is
#' @importFrom plyr . count
#' @import ShortRead
#' @import Rcpp
#' @import ggplot2
#' @useDynLib CellBarcode
NULL

#' A dummy BarcodeObj
#'
#' Dataset contains a BarcodeObj
#'
#' @name bc_obj
#' @docType data
#' @keywords dataset
#' @usage data(bc_obj)
#'
#' @details 
#' This is a BarcodeObj derived from dummy data.
#'
#' \preformatted{
#' d1 = data.frame(
#'     seq = c(
#'        "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'        "AATTCGATCGATCGAAGAGATCGATCGATC",
#'        "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'        "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'        "AAATCGATCGATCGAAGAGATCGATCGATC",
#'        "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'        "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'        "GGATCGATCGATCGAAGAGATCGATCGATC",
#'        "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'        "GGTTCGATCGATCGACGAGATCGATCGATC",
#'        "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'        ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#' 
#' d2 = data.frame(
#'     seq = c(
#'        "ACTTCGATCGATCGAAACGATCGATCGATC",
#'        "AATTCGATCGATCGAAGAGATCGATCGATC",
#'        "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'        "AAATCGATCGATCGAAGAGATCGATCGATC",
#'        "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'        "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'        "GGATCGATCGATCGAAGAGATCGATCGATC",
#'        "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'        "GGTTCGATCGATCGACGAGATCGATCGATC",
#'        "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'        ),
#'     freq = c(
#'         30, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#' 
#' pattern = "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj = bc_extract(
#'     list(test1 = d1, test2 = d2), 
#'     pattern, sample_name=c("test1", "test2"))
#'
#' bc_obj = bc_cure_depth(bc_obj, depth=5)
#'
#' # save the dummy data
#' # save(bc_obj, file = "./data/bc_obj.RData")
#' ###
#' }
NULL
