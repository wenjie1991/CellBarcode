# TODO: Remove barcode reads distribution graph with log x-axis
# Correlation between two samples w/o distribution information



#' Cleaning sequences
#'
#' @param barcodeObj A BarcodeObj
#' @param depth_threshold A Numeric, only a sequence with reads depth larger than that will be processed. If a UMI is used, the depth_threshold is the reads threshold for the UMI-sequence amplicon.
#' @param with_umi A bool value, True when UMI is applied.
#' @param umi_depth A Numeric, the minimum UMI-sequence amplicon counts needed for consider a true sequence. 
#' @param doFish A bool value, the UMI-sequence amplicon unsatisfied the depth_threshold are counted when when the sequence is in the final sequence list which satisfied the umi_depth threshold.
#' @param isUniqueUMI A bool value, True if UMI is considered absolute unique, if so only the dominant sequence is chosen for each UMI.
#' @param hammer_dist_threshold A integer, the less abundent sequence are compared to the most abundent ones and the hammering distance is calculated, the less abundent seuqneces with hammering distance less than the value will be removed.
#' @param barcode_count_threshold A integer, appliable when hammer_dist_threshold > 0, the max sequences number should be, when the number of sequence are satisfied all the rest less abundent sequences are removed. 
#' @return A BarcodeObj
#'
#' @export
bc_cure = function(
  barcodeObj
  , depth_threshold = 0 
  , with_umi = F
  , umi_depth = 2
  , doFish = F
  , isUniqueUMI = F
  , hammer_dist_threshold = 0
  , barcode_count_threshold = 10000
  ) {


  count = barcode_seq = NULL  # due to NOTE in check
  # TODO: Use the barcode distribution and sequence base pair quality to correct

  messyBc = barcodeObj$messyBc

  parameter_df = data.frame(
    sample_names = names(messyBc),
    depth_threshold = depth_threshold,
    with_umi = with_umi,
    umi_depth = umi_depth,
    doFish = doFish,
    isUniqueUMI = isUniqueUMI,
    barcode_count_threshold = barcode_count_threshold
    )

  if (with_umi) {
    ## When UMI used, count the UMI-barcode

    messyBc = lapply(1:length(messyBc), function(i) {
      d = messyBc[[i]]
      d0 = data.table(d)
      d1 = d0[count > parameter_df[i, "umi_depth"]]
      if (isUniqueUMI) {
        d1 = d1[count > parameter_df[i, "umi_depth"], .(barcode_seq = barcode_seq[which.max(count)], count = max(count)), by = umi_seq]
      } else {
        d1 = d1[count > parameter_df[i, "umi_depth"]]
      }
      if (doFish) {
        d1 = d0[barcode_seq %in% d1$barcode_seq, .(count = .N), by = barcode_seq]
      } else {
        d1 = d1[, .(count = .N), by = barcode_seq]
      }
      d1
    })
  } 
  ## If no UMI, count the reads directly
  messyBc = lapply(1:length(messyBc),
    function(i) {
      d = data.table(messyBc[[i]])
      ## TODO: If the output is empty (with 0 row) ...
      d[, .(count = sum(count)), by = barcode_seq][count > parameter_df[i, "depth_threshold"]]
    }
  )

  ## TODO: the barcode sequence clustering and correction
  ## Do the correction
  if (hammer_dist_threshold > 0) {
    correct_out = lapply(1:length(messyBc),
      function(i) {
        d = messyBc[[i]]
        seq_v = d$barcode_seq
        count_v = d$count
        seq_correct(seq_v, count_v, parameter_df[i, "barcode_count_threshold"], hammer_dist_threshold)
      }
    )
    ## Prepare output
    cleanBc = lapply(correct_out,
      function(d) {
        ##  The result is default ordered
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

   ## The result is default ordered
    cleanBc = lapply(messyBc, function(x) { x[order(count, decreasing = T)] })
    #     cleanProc = NULL
  }

  ## save the result
  names(cleanBc) = rownames(barcodeObj$meta_data)
  barcodeObj$cleanBc = cleanBc
  #   barcodeObj$cleanProc = cleanProc
  barcodeObj
}

