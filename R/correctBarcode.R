# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information

#' Filters barcodes by counts
#'
#' bc_cure_depth filters barcodes by applying the filter on the read counts or
#' the UMI counts.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A single or vector of numeric, specifying the threshold
#' of the threshold of minimum count for a sequence to kept. If the input is a
#' vector, the vector length should be the same to the sample number, and a 
#' numeric per sample.
#' @param isUpdate A logical value. If TRUE, the cleanBc element in BarcodeObj
#' will be used preferentially, otherwise the messyBc element will be used. If no
#' cleanBc is available, messyBc will be used instead.
#' @return A BarcodeObj object with cleanBc element updated.
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
#' # Remove barcodes with depth < 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#'
#' # Use UMI information to filter the barcode < 5 UMI-barcode tags
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
                by = barcode_seq][count >= parameter_df[i, "depth"]]
        }
    )

    ## save the result
    names(cleanBc) <- rownames(barcodeObj$metadata)
    barcodeObj$cleanBc <- cleanBc
    #   barcodeObj$cleanProc = cleanProc
    barcodeObj
}

#' Merges barcodes by editing distance
#'
#' bc_cure_cluster performs the clustering of barcodes by editing distance,
#' then merge the barcodes with similar sequence. This function is only
#' applicable for the BarcodeObj object with a cleanBc element.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param dist_thresh: a single integer or vector of integers with the length of sample
#' count, specifying the editing distance threshold of merging two similar
#' barcode sequences. If the input is a vector, each value in vector is for
#' one sample according to the sample order in BarcodeObj object.
#' @param dist_method: A  character string, specifying the distance algorithm for
#' evaluating barcodes similarity. It can be "hamm" for Hamming distance or
#' "leven" for Levenshtein distance.
#' @param merge_method: A character string specifying the algorithm used to perform the
#' clustering merging of barcodes. Currently only "greedy" is available, in this
#' case, the least abundant barcode is preferentially merged to the most
#' abundant ones.
#' @param barcode_n: A single integer or vector of integers, specifying the max number
#' of sequences expected. When the most abundant barcode number reaches this
#' number the merging finishes, and all the rest sub-abundant sequences are
#' discarded.
#' @param dist_costs: A list, the cost of the events when calculating distance between
#' two barcode sequences, applicable when Levenshtein distance is applied. The
#' names of vector have to be “insert”, “delete” and “replace”, specifying the
#' weight of insertion, deletion, replacement events respectively. The default
#' cost for each event is 1.
#' @return A BarcodeObj object with cleanBc element updated.
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
#' # Remove barcodes with depth < 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' 
#' # Do the clustering, merge the less abundent barcodes to the more abundent
#' # one by hamming distance <= 1 
#' bc_cure_cluster(bc_cured, dist_thresh = 1)
#' 
#' # Levenshtein distance <= 1
#' bc_cure_cluster(bc_cured, dist_thresh = 2, dist_method = "leven",
#'     dist_costs = list("insert" = 2, "replace" = 1, "delete" = 2))
#' 
#' ###
#' @export
bc_cure_cluster <- function(
    barcodeObj
    , dist_thresh = 1
    , dist_method = "hamm"
    , merge_method = "greedy"
    , barcode_n = 10000
    , dist_costs = list("replace" = 1, "insert" = 1, "delete" = 1)
    ) {
    # TODO: Add more clustering methods

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata)
        , distance =dist_thresh 
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

#' Filters on UMI-barcode tags counts when UMI used
#'
#' When the UMI is used, bc_cure_umi applies the filtering on the UMI-barcode tags counts. 
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A single or a vector of numeric, specifying the UMI-sequence
#' tags count threshold. Only the barcodes with UMI-barcode tag count larger than
#' the threshold are considered true barcodes.
#' @param doFish A single or a vector of logical value. If TRUE, for barcodes
#' with UMI read depth above the threshold, “fish” for identical barcodes with
#' UMI read depth below the threshold. The consequence of “doFish” will not
#' increase the number of identified barcodes, but the UMI counts will increase
#' due to including the low depth UMI barcodes. 
#' @param isUniqueUMI A single or a vector of logical value. When a UMI
#' relates to several barcodes. If you believe that the UMI is absolute unique,
#' then only the dominant sequence is chosen for each UMI.
#' @return A BarcodeObj object with cleanBc element updated (or created).
#' @details When invoke this function, the order of each steps are:
#' \enumerate{
#'   \item (optional when isUniqueUMI is TRUE) Find dominant sequence in each UMI.
#'   \item UMI-barcode depth filtering.
#'   \item (optional when doFish is TRUE) Fishing the UMI with low UMI-barcode depth.
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
#' # Use UMI information to remove the barcode < 5 UMI-barcode tags
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

        d1 = d0[count >= parameter_df[i, "umi_depth"]]

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


