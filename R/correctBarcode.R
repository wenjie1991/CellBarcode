# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information

#' Filter barcode by count
#'
#' The function applies the filter on the reads count or the UMI counts.
#'
#' @param barcodeObj A BarcodeObj
#' @param depth Numeric, only a sequence with counts larger than the threshold
#' will be kept.
#' @param isUpdate Bool. The function will use cleanBc if it exists, when it's
#' false, the messyBc will be used instead cleanBc to perform the filtering. In
#' the case that no cleanBc is available, messyBc is used.
#' @return A BarcodeObj
#'
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'  seq = c(
#'    "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AATTCGATCGATCGAAGAGATCGATCGATC",
#'    "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'    "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AAATCGATCGATCGAAGAGATCGATCGATC",
#'    "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'    "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'    "GGATCGATCGATCGAAGAGATCGATCGATC",
#'    "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'    "GGTTCGATCGATCGACGAGATCGATCGATC",
#'    "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'    ),
#'  freq = c(
#'    30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'    )
#'  )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'   pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth <= 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#' 
#' # Do the clustering, integre the less abundent barcodes to the more abundent
#' # one by hamming distance <= 1
#' bc_cure_cluster(bc_cured, distance = 1)
#'
#' # Use UMI information to filter the barcode <= 5 UMI-barcode tags
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#' ##
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
      d[, .(count = sum(count)), by = barcode_seq][count > parameter_df[i, "depth"]]
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
#' This function performs the clustering to merge the barcodes with similar sequence.
#' This function only appliable to the BarcodeObj with cleanBc element.
#'
#' @param barcodeObj A BarcodeObj
#' @param distance Integer, the editing distance threshold if two sequence is
#' similar enough to be merged
#' @param dist_method A String, currently it only be "hamm". The hamming
#' distance is used to evaluate the similarity between two barcodes.
#' @param merge_method A String, currently only "greedy" is available. Comparing
#' the least abundent barcode to the most abundent ones, merge the least
#' abundent to the most abundent ones.
#' @param barcode_n Integer, the max sequences number should be. When the most
#' abundent barcodes are satisfied this number the merging finished, and all the
#' rest sub-abundent sequences are removed. 
#' @return A BarcodeObj
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'  seq = c(
#'    "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AATTCGATCGATCGAAGAGATCGATCGATC",
#'    "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'    "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AAATCGATCGATCGAAGAGATCGATCGATC",
#'    "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'    "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'    "GGATCGATCGATCGAAGAGATCGATCGATC",
#'    "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'    "GGTTCGATCGATCGACGAGATCGATCGATC",
#'    "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'    ),
#'  freq = c(
#'    30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'    )
#'  )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'   pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth <= 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#' 
#' # Do the clustering, integre the less abundent barcodes to the more abundent
#' # one by hamming distance <= 1 
#' bc_cure_cluster(bc_cured, distance = 1)
#'
#' # Use UMI information to filter the barcode <= 5 UMI-barcode tags
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#'
#' @export
bc_cure_cluster <- function(
  barcodeObj
  , distance = 1
  , dist_method = "hamm"
  , merge_method = "greedy"
  , barcode_n = 10000
  ) {
  # TODO: Add more clustering methods

  parameter_df <- data.frame(
    sample_names = rownames(barcodeObj$metadata)
    , distance = distance
    , barcode_n = barcode_n
  )

  cleanBc <- barcodeObj$cleanBc

  if (dist_method == "hamm" & merge_method == "greedy") {
    correct_out <- lapply(1:length(cleanBc),
      function(i) {
        d <- cleanBc[[i]]
        seq_v <- d$barcode_seq
        count_v <- d$count

        # run hamming clustering
        seq_correct(
          seq_v, 
          count_v, 
          parameter_df[i, "barcode_n"], 
          parameter_df[i, "distance"]
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

    ## The correction log
    #     cleanProc <- lapply(correct_out,
    #       function(d) {
    #         d$link_table
    #       }
    #     )

    names(cleanBc) <- rownames(barcodeObj$metadata)
    # TODO: The cleanProc data structure
    #     names(cleanProc) <- rownames(barcodeObj$metadata)
  } else {
    stop("dist_method or merge_method is not valid.")
  }

  barcodeObj$cleanBc <- cleanBc
  barcodeObj
}

#' Filter on UMI-barcode 
#'
#' When the UMI is used and extracted by bc_extract, this function applies the
#' filtering on the UMI-barcode tags.  This function read in data in messyBc
#' element and create a cleanBc element in BarcodeObj.
#'
#' @param barcodeObj A BarcodeObj
#' @param depth Numeric, the UMI-sequence tags counts threshold. Only the
#' barcode with UMI-barcode larger than the threshold are considered potential
#' true barcodes.
#' @param doFish Bool value, the "fishing" process is re-counting the UMI with
#' potential true barcodes, and including those UMI-barcodes tags which are not
#' satisfied depth threshold.
#' @param isUniqueUMI Bool value, In the case that a UMI are with several
#' barcodes. If you believe that the UMI is absolute unique, then only the
#' dominant sequence is chosen for each UMI.
#' @return A BarcodeObj
#' @details The priority of potential steps are:
#'   1. (optional) Find dominant sequence in each UMI.
#'   2. UMI-barcode depth filtering.
#'   3. (optional) Fising the UMI with low UMI-barcode depth.
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'  seq = c(
#'    "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AATTCGATCGATCGAAGAGATCGATCGATC",
#'    "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'    "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'    "AAATCGATCGATCGAAGAGATCGATCGATC",
#'    "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'    "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'    "GGATCGATCGATCGAAGAGATCGATCGATC",
#'    "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'    "GGTTCGATCGATCGACGAGATCGATCGATC",
#'    "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'    ),
#'  freq = c(
#'    30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'    )
#'  )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'   pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth <= 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#' 
#' # Do the clustering, integre the less abundent barcodes to the more abundent
#' # one by hamming  distance <= 1 
#' bc_cure_cluster(bc_cured, distance = 1)
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
  cleanBc <- lapply(1:length(messyBc), function(i) {
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


