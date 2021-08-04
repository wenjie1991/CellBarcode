# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information

#' Filter barcodes by count
#'
#' The function filter barcodes by applying the filter on the reads count or the
#' UMI count.
#'
#' @param barcodeObj A BarcodeObj.
#' @param depth A single or vector of numeric, specifying the threshold
#' of the threshold of minimum count for a sequence to kept.
#' @param isUpdate A bool, if TRUE, the cleanBc element in BarcodeObj will
#' be used, otherwise the messyBc element will be used. In the case that no
#' cleanBc is available, messyBc will always be used.
#' @return A BarcodeObj with cleanBc update.
#'
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth <= 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#' 
#' # Use UMI information to filter the barcode <= 5 UMI-barcode tags
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#'
#' ### 
#' @export
bc_cure_depth <- function(
    barcodeObj
    , depth = 0 
    , isUpdate = TRUE
    ){

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata)
        , depth = depth
    )

    if (is.null(barcodeObj$cleanBc) | !isUpdate) {
        cleanBc <- barcodeObj$messyBc
    } else {
        cleanBc <- barcodeObj$cleanBc
    }

    ## count the reads directly
    cleanBc <- lapply(seq_along(cleanBc),
        function(i) {
            d <- data.table(cleanBc[[i]])
            ## TODO: If the output is empty (with 0 row) ...
            d[, 
                .(count = sum(count)), 
                by = barcode_seq][count > parameter_df[i, "depth"]]
        }
    )

    ## save the result
    names(cleanBc) <- rownames(barcodeObj$metadata)
    barcodeObj$cleanBc <- cleanBc
    #   barcodeObj$cleanProc = cleanProc
    barcodeObj
}

#' Merge barcodes by editing distance
#'
#' This function performs the clustering to merge the barcodes with similar
#' sequence. This function is only appliable to the BarcodeObj with a cleanBc
#' element.
#'
#' @param barcodeObj A BarcodeObj.
#' @param distance A sigle or a vector of integer, specifying he editing
#' distance threshold if two sequence is similar enough to be merged.
#' @param dist_method A  character string specifying the distance algorithm used
#' to evaluate the barcodes similarity. It can be "hamm" for Hamming distance or
#' "leven" for Levenshtein distance.
#' @param merge_method A character string specifying the algorithm used to
#' perform the clustering merging of two barcodes. Currently only "greedy" is
#' available, the least abundent barcode to the most abundent ones, merge the
#' least abundent to the most abundent ones.
#' @param barcode_n A single or vector of integer, specifying the max sequences
#' number should be. When the most abundent barcodes are satisfied this number
#' the merging finished, and all the rest sub-abundent sequences are removed. 
#' @param dist_costs A list, the costs of the events when calculate distance
#' between to barcode sequences, applicable when Levenshtein distance is 
#' applied. The names of vector can be of "insert", "delete" and 
#' "replace". The default cost is 1.
#' @return A BarcodeObj with cleanBc element updated.
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth <= 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' 
#' # Do the clustering, merge the less abundent barcodes to the more abundent
#' # one by hamming distance <= 1 
#' bc_cure_cluster(bc_cured, distance = 1)
#' 
#' # Levenshtein distance <= 1
#' bc_cure_cluster(bc_cured, distance = 2, dist_method = "leven",
#'     dist_costs = list("insert" = 2, "replace" = 1, "delete" = 2))
#' 
#' ###
#' @export
bc_cure_cluster <- function(
    barcodeObj
    , distance = 1
    , dist_method = "hamm"
    , merge_method = "greedy"
    , barcode_n = 10000
    , dist_costs = list("replace" = 1, "insert" = 1, "delete" = 1)
    ) {
    # TODO: Add more clustering methods

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata)
        , distance = distance
        , barcode_n = barcode_n
    )

    cleanBc <- barcodeObj$cleanBc

    if (dist_method == "hamm" & merge_method == "greedy") {
        correct_out <- lapply(seq_along(cleanBc),
            function(i) {
                d <- cleanBc[[i]]
                seq_v <- d$barcode_seq
                count_v <- d$count

                # run hamming clustering
                seq_correct(
                    seq_v, 
                    count_v, 
                    parameter_df[i, "barcode_n"], 
                    parameter_df[i, "distance"],
                    1
                )
            }
        )

        ## Prepare output
        cleanBc <- lapply(correct_out,
            function(d) {
                ##  The result is default ordered
                data.frame(
                    d$seq_freq[order(d$seq_freq$count, decreasing = TRUE), ]
                )
            }
        )

        ## The correction log
        #     cleanProc <- lapply(correct_out,
        #       function(d) {
        #         d$link_table
        #       }
        #     )

        names(cleanBc) <- rownames(barcodeObj$metadata)
        # TODO: The cleanProc data structure
        #     names(cleanProc) <- rownames(barcodeObj$metadata)
    } else if (dist_method == "leven" & merge_method == "greedy") {

        insert_costs = ifelse(is.null(dist_costs$insert), 1, dist_costs$insert)
        delete_costs = ifelse(is.null(dist_costs$delete), 1, dist_costs$delete)
        replace_costs = ifelse(is.null(dist_costs$replace), 1, dist_costs$replace)

        correct_out <- lapply(seq_along(cleanBc),
            function(i) {
                d <- cleanBc[[i]]
                seq_v <- d$barcode_seq
                count_v <- d$count

                # run levenshtein clustering
                seq_correct(
                    seq_v, 
                    count_v, 
                    parameter_df[i, "barcode_n"], 
                    parameter_df[i, "distance"],
                    2,
                    insert_costs,
                    delete_costs,
                    replace_costs
                )
            }
        )

        ## Prepare output
        cleanBc <- lapply(correct_out,
            function(d) {
                ##  The result is default ordered
                data.frame(d$seq_freq[order(d$seq_freq$count, decreasing = TRUE), ])
            }
        )

        names(cleanBc) <- rownames(barcodeObj$metadata)
    } else {
        stop("dist_method or merge_method is not valid.")
    }

    barcodeObj$cleanBc <- cleanBc
    barcodeObj
}

#' Filter on UMI-barcode tags
#'
#' When the UMI is used and identified by bc_extract, this function applies
#' the filtering on the UMI-barcode tags.  This function read in data in
#' messyBc element and create a cleanBc element in BarcodeObj.
#'
#' @param barcodeObj A BarcodeObj
#' @param depth A single or a vector of numeric, specifying the UMI-sequence
#' tags counts threshold. Only the barcode with UMI-barcode larger than the
#' threshold are considered barcodes.
#' @param doFish A single or a vector of bool value, if TRUE, the "fishing"
#' process will be applied to re-counting the UMI with  barcodes, which are not
#' satisfied depth threshold.
#' @param isUniqueUMI A single or a vector of bool value. In the case that a UMI
#' relates to several barcodes. If you believe that the UMI is absolute unique,
#' then only the dominant sequence is chosen for each UMI.
#' @return A BarcodeObj
#' @details When invoke this function, the order of each steps are:
#' \enumerate{
#'   \item (optional) Find dominant sequence in each UMI.
#'   \item UMI-barcode depth filtering.
#'   \item (optional) Fising the UMI with low UMI-barcode depth.
#' }
#'
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'    seq = c(
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
#'    freq = c(
#'        30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'        )
#'    )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Use UMI information to filter the barcode <= 5 UMI-barcode tags
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#'
#' @export
bc_cure_umi <- function(
    barcodeObj
    , depth = 2
    , doFish = FALSE
    , isUniqueUMI = FALSE
    ) {

    messyBc <- barcodeObj$messyBc

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata),
        umi_depth = depth,
        doFish = doFish,
        isUniqueUMI = isUniqueUMI
    )

    ## When UMI used, count the UMI-barcode
    cleanBc <- lapply(seq_along(messyBc), function(i) {
        d = messyBc[[i]]
        d0 = data.table(d)
        # filter by umi depth

        if (isUniqueUMI) {
            # dominant sequence is chosen for each UMI.
            d0 <- d0[, 
                .(barcode_seq = barcode_seq[which.max(count)], 
                    count = max(count)), by = umi_seq]
        } 

        d1 = d0[count > parameter_df[i, "umi_depth"]]

        if (doFish) {
            # including umi with true barcodes do not meet depth threshold
            d1 <- d0[barcode_seq %in% d1$barcode_seq, .(count = .N), by = barcode_seq]
        } else {
            d1 <- d1[, .(count = .N), by = barcode_seq]
        }
        d1
    })

    names(cleanBc) <- rownames(barcodeObj$metadata)
    barcodeObj$cleanBc <- cleanBc
    barcodeObj
}


