#' Extract barcode from reads
#'
#' @param barcodeObj A BarcodeObj
#' @param pattern A string. The regular expression to match the barcode which capture pattern
#' @param maxLDist A integer. The mismatch threshold for barcode matching
#' @param pattern_type A vector. It defines the barcode (and UMI) capture pattern
#' @param costs A vector. Define the weight for each mismatch events
#' @return A barcodeObj
#'
#' @export
bc_extract = function(x, ...) UseMethod("bc_extract", x)

# report_barcode_extraction = function(sample_name, reads_num, barcode_reads) {
#   if (is.null(sample_name)) {
#     cat("In this sample, ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads ,".\n", sep = "")
#   } else {
#     cat("In sample ", sample_name, ", ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads_num ,".\n", sep = "")
#   }
# }

#' Extract barcode from ShortReadQ
#'
#' @export
bc_extract.ShortReadQ = function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = T) {
  reads_freq = ShortRead::tables(x, n=Inf)$top
  reads_seq = names(reads_freq)
  m = utils::aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = costs)
  seq_l = regmatches(reads_seq, m)

  seq_v = plyr::laply(seq_l, function(x) { x[1] })
  barcode_v = plyr::laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v = plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
  if (ordered) {
    d = d[order(count, decreasing = T)]
  } 
  if (verbose) {
    # BUG: The reads_num is not the total reads number
    #  uences    report_barcode_extraction(sample_name = sample_name, reads_num = sum(reads_freq), barcode_reads = sum(d$count, na.omit = T))
  }
  stats::na.omit(d)
}

#' Extract barcode from DNAStringSet
#'
#' @export
bc_extract.DNAStringSet = function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = T) {
  bc_extract.ShortReadQ(x, pattern = pattern, sample_name = sample_name, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = verbose)
}

#' Extract barcode from data.frame
#'
#' @export
bc_extract.data.frame = function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = T) {
  sequences = x$seq
  freq = x$freq
  x = DNAStringSet(rep(sequences, freq))

  bc_extract(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = verbose)
}

#' Extract barcode from character vector
#'
#' @export
bc_extract.character = function(file, pattern = "", sample_name = basename(file), meta_data = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = T) {
  if (length(file) > 1) {
    if (is.null(sample_name)) {
      sample_name = names(x)
    }
    if (is.null(sample_name)) {
      sample_name = 1:length(x)
    }
    if (is.null(meta_data)) {
      meta_data = data.frame(sample_name = sample_name)
      rownames(meta_data) = sample_name
    }

    messyBc = lapply(1:length(file), function(i) {
      bc_extract(ShortRead::readFastq(file[i]), sample_name = sample_name[i], pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = verbose)
    })

    names(messyBc) = sample_name
    output = list(messyBc = messyBc, meta_data = meta_data)
    class(output) = "BarcodeObj"
    return(output)
  } else {
    bc_extract(ShortRead::readFastq(file), sample_name = sample_name[i], pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = verbose)
  }
}

#' Extract barcode from integer vector
#'
#' @export
bc_extract.integer = function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = T) {
  # TODO: Table result
  x = DNAStringSet(rep(names(x), x))

  bc_extract(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = verbose)
}

#' Extract barcode from list
#'
#' @export
bc_extract.list = function(x, pattern = "", sample_name = NULL, meta_data = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = T, verbose = F) {
  if (is.null(sample_name)) {
    sample_name = names(x)
  }
  if (is.null(sample_name)) {
    sample_name = 1:length(x)
  }
  if (is.null(meta_data)) {
    meta_data = data.frame(sample_name = sample_name)
    rownames(meta_data) = sample_name
  }

  parallel::mclapply(1:length(x), function(i) {
    bc_extract(x[[i]], pattern = pattern, sample_name = sample_name[i], maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered, verbose = F)
  }) -> messyBc

  names(messyBc) = sample_name
  output = list(messyBc = messyBc, meta_data = meta_data)
  class(output) = "BarcodeObj"
  output
}

