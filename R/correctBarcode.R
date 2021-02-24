#' Stair plot used to find a optimized cutoff point
stairPlot = function(x, ...) UseMethod("stairPlot", x)

#' Stair plot used to find a optimized cutoff point
#'
#' @param barcodeObj A BarcodeObj with extracted barcode.
#' @return NULL
#' @export
stairPlot.BarcodeObj = function(barcodeObj, sample = 1, max = 1000) {
  # TODO: The stair plot only apply to messyBc iterm, do we need to apply it to cleanBc
  #       Return the optimized cutoff point.
  messyBc = barcodeObj$messyBc
  d = data.table::data.table(messyBc[[sample]])
  d = d[, .(count = sum(count)), by = barcode_seq][order(count)]
  y = cumsum(d$count)
  x = d$count
  graphics::plot(x, log10(y), type = 's', xlab = "Reads count", ylab = "log10 Accumlated Reads count", xlim = c(0, max))
}

#' Correct the PCR and sequencing mutation in the barcode sequence
cureBc = function(x, ...) UseMethod("cureBc", x)

#' Correct the PCR and sequencing mutation in the barcode sequence
#'
#' @param barcodeObj A BarcodeObj
#' @param depth_threshold A integer, only a barcode with sequencing depth larger than that will be processed
#' @param hamming_dist_threshold A integer, the hamming distance threshold of two distinct barcode
#' @param with_umi A bool value, True when UMI is used
#' @return A BarcodeObj
#' @export
cureBc.BarcodeObj = function(
  barcodeObj
  , depth_threshold = 0
  , barcode_count_threshold = 1000
  , hammer_dist_threshold = 1
  , with_umi = F
  ) {
  # TODO: Use the barcode distribution and sequence base pair quality to correct

  messyBc = barcodeObj$messyBc

  if (with_umi) {
    ## When UMI used, count the UMI-barcode
    messyBc = lapply(messyBc,
      function(d) {
        d = data.table(d)
        d = d[count > depth_threshold, .(count = .N), by = .(barcode_seq)]
      }
    )
  } else {
    ## If no UMI, count the reads
    messyBc = lapply(messyBc,
      function(d) {
        d = data.table(d)
        d = d[, .(count = sum(count)), by = barcode_seq][count > depth_threshold]
      }
    )
  }
  ## Do the correction
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

  ## The correction log, in another words, who correction has been done
  cleanProc = lapply(correct_out,
    function(d) {
      d$link_table
    }
  )

  names(cleanBc) = names(messyBc)
  names(cleanProc) = names(messyBc)

  ## save the result
  barcodeObj$cleanBc = cleanBc
  barcodeObj$cleanProc = cleanProc
  barcodeObj
}
