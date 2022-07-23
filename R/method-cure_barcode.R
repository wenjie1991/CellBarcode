

## private functions
###########################
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


#' @rdname bc_auto_cutoff
#' @exportMethod bc_auto_cutoff
setMethod("bc_auto_cutoff", c("BarcodeObj"), function(barcodeObj, useCleanBc=TRUE) {

    if (is.null(barcodeObj@cleanBc) | !useCleanBc) {
        x <- barcodeObj@messyBc
        message("------------\nbc_auto_cutoff: messyBc is used.\n------------")
    } else {
        x <- barcodeObj@cleanBc
        message("------------\nbc_auto_cutoff: cleanBc is used.\n------------")
    }

    res <- vapply(x, function(x_i) {
        bc_find_depth_cutoff_point(x_i$count)
    }, c(1.0))
    names(res) <- bc_names(barcodeObj)
    res
})

#' @rdname bc_cure_depth
#' @exportMethod bc_cure_depth
setMethod("bc_cure_depth", c("BarcodeObj"), function(
    barcodeObj
    , depth = 0
    , isUpdate = TRUE
    ){

    if (depth[1] < 0)
        depth <- NULL

    if (isUpdate) {
        message("------------\nbc_cure_depth: isUpdate is TRUE, update the cleanBc.\n------------")
    } else {
        message("------------\nbc_cure_depth: isUpdate is FALSE, use messyBc as input.\n------------")
    }
    
    if (is.null(barcodeObj@cleanBc) | !isUpdate) {
        cleanBc <- barcodeObj@messyBc
    } else {
        cleanBc <- barcodeObj@cleanBc
    }

    if (is.null(depth)) {
        message("------------\nbc_cure_depth: Null depth or negative provided, apply auto depth threshold.\n------------")
        depth <- vapply(cleanBc, function(x_i) {
                bc_find_depth_cutoff_point(x_i$count)
        }, c(1.0))
    }

    bc_meta(barcodeObj, "depth_cutoff") <- depth

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj@metadata)
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
    names(cleanBc) <- rownames(barcodeObj@metadata)
    barcodeObj@cleanBc <- cleanBc
    #   barcodeObj@cleanProc = cleanProc
    barcodeObj
})


#' @rdname bc_cure_cluster
#' @exportMethod bc_cure_cluster
setMethod("bc_cure_cluster", c("BarcodeObj"), function(
    barcodeObj
    , dist_thresh = 1
    , depth_fold_threshold = 1
    , dist_method = "hamm"
    , cluster_method = "greedy"
    , count_threshold = 1e7
    , dist_costs = list("replace" = 1, "insert" = 1, "delete" = 1)
    ) {
    # TODO: Add more clustering methods

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj@metadata)
        , distance =dist_thresh 
        , depth_fold_threshold = depth_fold_threshold
        , count_threshold =count_threshold 
    )

    cleanBc <- barcodeObj@cleanBc

    if (dist_method == "hamm" & cluster_method == "greedy") {
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
                    parameter_df[i, "depth_fold_threshold"], 
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

        names(cleanBc) <- rownames(barcodeObj@metadata)
        # TODO: The cleanProc data structure
        #     names(cleanProc) <- rownames(barcodeObj@metadata)
    } else if (dist_method == "leven" & cluster_method == "greedy") {

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
                    parameter_df[i, "depth_fold_threshold"], 
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

        names(cleanBc) <- rownames(barcodeObj@metadata)
    } else {
        stop("dist_method or cluster_method is not valid.")
    }

    barcodeObj@cleanBc <- cleanBc
    barcodeObj
})


#' @rdname bc_cure_umi
#' @exportMethod bc_cure_umi
setMethod("bc_cure_umi", c("BarcodeObj"), function(
    barcodeObj
    , depth = 2
    , doFish = FALSE
    , isUniqueUMI = FALSE
    ) {

    messyBc <- barcodeObj@messyBc

    parameter_df <- data.frame(
        sample_names = rownames(barcodeObj@metadata),
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

    names(cleanBc) <- rownames(barcodeObj@metadata)
    barcodeObj@cleanBc <- cleanBc
    barcodeObj
})
