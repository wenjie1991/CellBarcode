# TODO:
# One function do all the job
# filterFq() -> output the filtering detail
# make function to cross talk with genBaRcode and other packages
# create function to identify the cutoff point
# [ ] evaluate the sequencing depth (saturation of the sequene)
# [y] sequence qc output: data.frame with 1. total reads, 2. mean length
# [y] bc_extract add raw_read_count, barcode_read_count to metadata
# [y] bc_cure_depth auto_identify the depth cutoff point
# [y] Add cutoff point in metadata
# [ ] bc_plot_pairwise with interactive available (plotly)
# [ ] bc_plot_pairwise with interactive available (d3)
# [y] bc_plot_pairwise ggplot2
# [ ] bc_plot_single with interactive (plotly)
# [ ] bc_plot_single with interactive (d3)
# [y] bc_plot_single with ggplot2

# [y] Draw reads per barcode between samples pairwisely between technical replicates
#       or different samples, draw reads per barcodes between more than two samples.
# [ ] Explore the barcodes within a single sample or between two sample or more
#       interactively (explore the cutoff point interactively).
# [y] Merge technical repeats by merge sample function (subsetting barcodeObj ->
#       change samle names -> combine two or more barcodeObj).
# [ ] If no barcodes in the object, give errors when processing it with
# BarcodeObj methods.

# global help function

ifnullelse <- function(cand1, other) {
    if (is.null(cand1)) {
        return(other)
    } else {
        return(cand1)
    }
}

#' DNA Barcode Analysis toolkit
#'
#' This package performs DNA Barcode (genetic lineage tracing) analysis. The
#' package can handle all kinds of DNA barcodes, as long as the barcode within a
#' single sequencing read and has a pattern which can be matched by a regular
#' expression. \code{CellBarcode} can handle barcode with flexible length, with or
#' without UMI (unique molecular identifier). This tool also can be used for
#' pre-processing of some amplicon data such as CRISPR gRNA screening, immune
#' repertoire sequencing and meta genome data.
#'
#' @name CellBarcode
#' @docType package
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
#' # save the dummy data
#' # save(bc_obj, file = "./data/bc_obj.RData")
#' ###
#' }
"bc_obj"
