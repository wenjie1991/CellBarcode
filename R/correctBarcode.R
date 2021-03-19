# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information

#' Correct the PCR and sequencing mutation in the barcode sequence
#'
#' @export
bc_cure = function(x, ...) UseMethod("bc_cure", x)

#' Correct the PCR and sequencing mutation in the barcode sequence
#'
#' @param barcodeObj A BarcodeObj
#' @param depth_threshold A integer, only a barcode with sequencing depth larger than that will be processed
#' @param hamming_dist_threshold A integer, the hamming distance threshold of two distinct barcode
#' @param with_umi A bool value, True when UMI is used
#' @return A BarcodeObj
#' @export
bc_cure.BarcodeObj = function(
  barcodeObj
  , depth_threshold = 2
  , barcode_count_threshold = 10000
  , hammer_dist_threshold = 0
  , with_umi = F
  , doFish = T
  , isUniqueUMI = F
  ) {

  count = barcode_seq = NULL  # due to NOTE in check
  # TODO: Use the barcode distribution and sequence base pair quality to correct

  messyBc = barcodeObj$messyBc

  if (with_umi) {
    ## When UMI used, count the UMI-barcode
    messyBc = lapply(messyBc, function(d) {
      d0 = data.table(d)
      d1 = d0[count > depth_threshold]
      if (isUniqueUMI) {
        d1 = d1[count > depth_threshold, .(barcode_seq = barcode_seq[which.max(count)], count = max(count)), by = umi_seq]
      } else {
        d1 = d1[count > depth_threshold]
      }
      if (doFish) {
        d1 = d0[barcode_seq %in% d1$barcode_seq, .(count = .N), by = barcode_seq]
      } else {
        d1 = d1[, .(count = .N), by = barcode_seq]
      }
      d1
    })
  } else {
    ## If no UMI, count the reads
    messyBc = lapply(messyBc,
      function(d) {
        d = data.table(d)
        ## TODO: If the output is empty (with 0 row) ...
        d[, .(count = sum(count)), by = barcode_seq][count > depth_threshold]
      }
    )
  }

  ## Do the correction
  if (hammer_dist_threshold > 0) {
    correct_out = lapply(messyBc,
      function(d) {
        seq_v = d$barcode_seq
        count_v = d$count
        seq_correct(seq_v, count_v, barcode_count_threshold, hammer_dist_threshold)
      }
    )
    ## Prepare output
    cleanBc = lapply(correct_out,
      function(d) {
        data.frame(d$seq_freq[order(d$seq_freq$count, decreasing = T), ])
      }
    )

    ## The correction log
    cleanProc = lapply(correct_out,
      function(d) {
        d$link_table
      }
    )

    names(cleanBc) = names(messyBc)
    names(cleanProc) = names(messyBc)

  } else {

    cleanBc = messyBc
    cleanProc = NULL

  }


  ## save the result
  barcodeObj$cleanBc = cleanBc
  barcodeObj$cleanProc = cleanProc
  barcodeObj
}
