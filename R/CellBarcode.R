##  Global help function
##########################
ifnullelse <- function(cand1, other) {
    if (is.null(cand1)) {
        return(other)
    } else {
        return(cand1)
    }
}

# Check if a object is a pure data.frame
is_pure_dataframe <- function(obj) {
    inherits(obj, "data.frame") && !inherits(obj, "tbl") && !inherits(obj, "data.table")
}

#' DNA Barcode Analysis toolkit
#'
#' The package CellBarcode performs Cellular DNA Barcode analysis. It can handle all kinds of DNA barcodes, as long as the barcode is within a single sequencing read and has a pattern that can be matched by a regular expression. \code{CellBarcode} can handle barcodes with flexible lengths, with or without UMI (unique molecular identifier). This tool also can be used for pre-processing some amplicon data such as CRISPR gRNA screening, immune repertoire sequencing, and metagenome data.
#'
#' @name CellBarcode
#' @docType _PACKAGE
#' @importFrom magrittr %>% %<>% extract extract2
#' @importFrom data.table data.table rbindlist .N :=
#' @importFrom Biostrings readDNAStringSet
#' @importFrom plyr . count
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @importFrom egg ggarrange
#' @importFrom stringr str_glue str_match str_split_fixed str_c
#' @importFrom utils combn aregexec
#' @import ShortRead
#' @import stats
#' @import methods
#' @import Rcpp
#' @import ggplot2
#' @useDynLib CellBarcode
NULL

#' A dummy BarcodeObj object
#'
#' Dataset contains a \code{BarcodeObj} with makeup barcode data.
#'
#' @name bc_obj
#' @docType data
#' @keywords dataset
#' @usage data(bc_obj)
#' @format This is a BarcodeObj object
#' @return A \code{BarcodeObj} object.
#' @source
#' This is a \code{BarcodeObj} object derived from makeup data by:
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
#' # Save the dummy data
#' # save(bc_obj, file = "./data/bc_obj.RData")
#' ###
#' }
"bc_obj"
