# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information

# bc_filter_barcode = function(count, sequence) {
#     sequence = sequence[count != 0]
#     count = count[count != 0]
#     k = 2
#     x_m = mean(count)
#     count = log2(count + 1)
#     x_sub = count[count > log2(x_m)]
#     s_sub = sequence[count > log2(x_m)]
#     d_sub = data.table(count = 2^x_sub, barcode_seq = s_sub)
#     weight_log_reads = x_sub
#     result = Ckmeans.1d.dp::Ckmeans.1d.dp(x_sub, k, y = weight_log_reads, method = "linear")
#     d_res = d_sub[result$cluster == 2]
#     d_merge = merge(x_known_barcod, d_res, all = T)
#     count_sum = d_merge$count %>% sum(na.rm = T)
#     d_merge[, cell_out := count / count_sum * cell_number]
#     d_merge
# }
 
bc_find_depth_cutoff_point <- function(count, count_lower_bound = median(count[count > 0])) {

    count <- count[count != 0]
    k <- 2

    count_log <- log2(count + 1)

    count_sub_i <- count_log > log2(count_lower_bound + 1)
    count_log_sub <- count_log[count_sub_i]
    weight_log_reads <- count_log_sub

    result <- Ckmeans.1d.dp::Ckmeans.1d.dp(count_log_sub, k, y = weight_log_reads, method = "linear")
    
    # print(count)
    min(count[count_sub_i][result$cluster == 2])
}

#' Finds barcode count cutoff point
#'
#' Finds the cutoff point for the barcode count filtering based on the barcode
#' count distribution.
#'
#' The one dimension kmeans clustering is applied for identify the 
#' "true barcode" based on read count. The the algorithm detail is:
#' 1. Remove the barcodes with count below the median of counts.
#' 2. Transform the count by log2(x+1).
#' 3. Apply the 1 dimension clustering to the logarized count, with
#' the cluster number of 2 and weights of the logarized count.
#' 4. Choose the minimum count value in the cluster with higher count as
#' cutoff point.
#'
#' For more info about 1 dimension kmeans used here please refer to
#' \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}, which has been used here.
#'
#' @param barcodeObj A \code{BarcodeObj} object.
#' @param useCleanBc A logical value, if \code{TRUE}, the \code{cleanBc} element in the
#' \code{BarcodeObj} object will be used, otherwise the \code{messyBc} element will be used.
#' @return a numeric \code{vector} of the cutoff point.
#' @examples
#' 
#' data(bc_obj)
#' 
#' bc_auto_cutoff(bc_obj)
#' 
#' @export
bc_auto_cutoff <- function(barcodeObj, useCleanBc=TRUE) {

    if (is.null(barcodeObj$cleanBc) | !useCleanBc) {
        x <- barcodeObj$messyBc
        message("-message----\nbc_auto_cutoff: messyBc is used.\n------------")
    } else {
        x <- barcodeObj$cleanBc
        message("-message----\nbc_auto_cutoff: cleanBc is used.\n------------")
    }

    res <- vapply(x, function(x_i) {
        bc_find_depth_cutoff_point(x_i$count)
    }, c(1.0))
    names(res) <- bc_names(barcodeObj)
    res
}

#' Filters barcodes by counts
#'
#' bc_cure_depth filters barcodes by the read counts or the UMI counts.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A numeric or a vector of numeric, specifying the threshold of 
#' minimum count for a barcode to kept. If the input is a
#' vector, if the vector length is not the same to the sample number the element
#' will be repeatedly used. And when the depth argument is a number with negative
#' value, automatic cutoff point will be chosen by \code{bc_auto_cutoff} function
#' for each samples. See \code{\link[CellBarcode]{bc_auto_cutoff}} for details.
#' @param isUpdate A logical value. If TRUE, the \code{cleanBc} element in \code{BarcodeObj}
#' will be used preferentially, otherwise the \code{messyBc} element will be used. If no
#' cleanBc is available, \code{messyBc} will be used instead.
#' @return A \code{BarcodeObj} object with \code{cleanBc} element updated or
#' created.
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
#' # Use UMI information, filter the barcode < 5 UMI
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

    if (depth[1] < 0)
        depth <- NULL

    if (isUpdate) {
        message("-message----\nbc_cure_depth: isUpdate is TRUE, update the cleanBc.\n------------")
    } else {
        message("-message----\nbc_cure_depth: isUpdate is FALSE, use messyBc as input.\n------------")
    }
    
    if (is.null(barcodeObj$cleanBc) | !isUpdate) {
        cleanBc <- barcodeObj$messyBc
    } else {
        cleanBc <- barcodeObj$cleanBc
    }

    if (is.null(depth)) {
        message("-message----\nbc_cure_depth: Null depth or negative provided, apply auto depth threshold.\n------------")
        depth <- vapply(cleanBc, function(x_i) {
                bc_find_depth_cutoff_point(x_i$count)
        }, c(1.0))
    }

    bc_meta(barcodeObj, "depth_cutoff") <- depth

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata)
        , depth = depth
    )

    
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
#' \code{bc_cure_cluster} performs clustering of barcodes by editing distance,
#' and merging the barcodes with similar sequence. This function is only
#' applicable for the BarcodeObj object with a \code{cleanBc} element.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param dist_thresh A single integer or vector of integers with the length of sample
#' count, specifying the editing distance threshold of merging two similar
#' barcode sequences. If the input is a vector, each value in the vector relates to one 
#' sample according to the sample order in \code{BarcodeObj} object.
#' @param dist_method A  character string, specifying the distance algorithm used for
#' evaluating barcodes similarity. It can be "hamm" for Hamming distance or
#' "leven" for Levenshtein distance.
#' @param merge_method A character string specifying the algorithm used to perform the
#' clustering merging of barcodes. Currently only "greedy" is available, in this
#' case, the least abundant barcode is preferentially merged to the most
#' abundant ones.
#' @param count_threshold An integer, read depth threshold to consider a
#' barcode as a true barcode, when when a barcode with count higher than this
#' threshold it will not be merged into more abundant barcode.
#' @param dist_costs A list, the cost of the events of distance algorithm, 
#' applicable when Levenshtein distance is applied. The
#' names of vector have to be \code{insert}, \code{delete} and \code{replace}, specifying the
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
    , count_threshold = 1000
    , dist_costs = list("replace" = 1, "insert" = 1, "delete" = 1)
    ) {
    # TODO: Add more clustering methods

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj$metadata)
        , distance =dist_thresh 
        , count_threshold =count_threshold 
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
                    parameter_df[i, "count_threshold"], 
                    parameter_df[i, "distance"],
                    1
                )
            }
        )

        ## Prepare output
        cleanBc <- lapply(correct_out,
            function(d) {
                ##  The result is default ordered
                d$seq_freq[order(d$seq_freq$count, decreasing = TRUE), ]
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
                    parameter_df[i, "count_threshold"], 
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
                d$seq_freq[order(d$seq_freq$count, decreasing = TRUE), ]
            }
        )

        names(cleanBc) <- rownames(barcodeObj$metadata)
    } else {
        stop("dist_method or merge_method is not valid.")
    }

    barcodeObj$cleanBc <- cleanBc
    barcodeObj
}

#' Filters UMI-barcode tag by counts 
#'
#' When the UMI is applied, \code{bc_cure_umi} can filter the UMI-barcode tags by counts. 
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A numeric or a vector of numeric, specifying the UMI-barcode
#' tag count threshold. Only the barcodes with UMI-barcode tag count larger than
#' the threshold are kept. 
#' @param doFish A logical value, if true, for barcodes with UMI read depth
#' above the threshold, “fish” for identical barcodes with UMI read depth below
#' the threshold. The consequence of \code{doFish} will not increase the number of
#' identified barcodes, but the UMI counts will increase due to including the
#' low depth UMI barcodes. 
#' @param isUniqueUMI A logical value, In the case that a UMI
#' relates to several barcodes, if you believe that the UMI is absolute unique,
#' then only the UMI-barcodes tags with highest count are chosen for each UMI.
#' @return A \code{BarcodeObj} object with \code{cleanBc} element updated (or created).
#' @details When invoke this function, it processes the data with following
#' steps:
#' \enumerate{
#'   \item (if isUniqueUMI is TRUE) Find dominant sequence in each UMI.
#'   \item UMI-barcode depth filtering.
#'   \item (if doFish is TRUE) Fishing the UMI with low UMI-barcode depth.
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


