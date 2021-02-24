#' Extract barcode from reads
#'
#' @param barcodeObj A BarcodeObj
#' @param pattern A string. The regular expression to match the barcode which capture pattern
#' @param maxLDist A integer. The mismatch threshold for barcode matching
#' @param pattern_type A vector. It defines the barcode (and UMI) capture pattern
#' @param costs A vector. Define the weight for each mismatch events
#' @return A barcodeObj
#' @export
extractBc = function(x, ...) UseMethod("extractBc", x)


extractBc.ShortReadQ = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  reads_freq = ShortRead::tables(x, n=Inf)$top
  reads_seq = names(reads_freq)
  m = utils::aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = list(sub = 1))
  seq_l = regmatches(reads_seq, m)

  seq_v = laply(seq_l, function(x) { x[1] })
  barcode_v = laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v = laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
}

extractBc.DNAStringSet = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  reads_freq = ShortRead::tables(x, n=Inf)$top
  reads_seq = names(reads_freq)
  m = utils::aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = list(sub = 1))
  seq_l = regmatches(reads_seq, m)

  seq_v = laply(seq_l, function(x) { x[1] })
  barcode_v = laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v = laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
  stats::na.omit(d)
}


extractBc.data.frame = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  sequences = x$seq
  freq = x$freq
  x = DNAStringSet(rep(sequences, freq))

  extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs)
}

extractBc.character = function(file, sampleName = basename(file), pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  if (length(file) > 1) {
    bc_list = lapply(file, function(x) {
      extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs)
    })
    names(bc_list) = sampleName
    bc_list
  } else {
    x = ShortRead::readFastq(file)
    extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs)
  }
}

extractBc.integer = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  # TODO: Table result
  x = DNAStringSet(rep(names(x), x))

  extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs)
}

extractBc.list = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), sample_name = NULL, costs = list(sub = 1, ins = 99, del = 99)) {
  if (is.null(sample_name)) {
    sample_name = names(x)
  }
  mclapply(x, function(x_i) {
    extractBc(x_i, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs)
  }) -> messyBc
  names(messyBc) = sample_name
  output = list(messyBc = messyBc)
  class(output) = "BarcodeObj"
  output
}
