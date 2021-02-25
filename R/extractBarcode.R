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

report_barcode_extraction = function(sample_name, reads_num, barcode_reads) {
  if (is.null(sample_name)) {
    cat("In this sample, ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads ,".\n", sep = "")
  } else {
    cat("In sample ", sample_name, ", ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads_num ,".\n", sep = "")
  }
}

extractBc.ShortReadQ = function(x, pattern = "", sample_name = NULL, maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = T) {
  reads_freq = ShortRead::tables(x, n=Inf)$top
  reads_seq = names(reads_freq)
  m = utils::aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = list(sub = 1))
  seq_l = regmatches(reads_seq, m)

  seq_v = plyr::laply(seq_l, function(x) { x[1] })
  barcode_v = plyr::laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v = plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
  if (verbose) {
    report_barcode_extraction(sample_name = sample_name, reads_num = sum(reads_freq), barcode_reads = sum(d$count, na.omit = T))
  }
  stats::na.omit(d)
}

extractBc.DNAStringSet = function(x, pattern = "", sample_name = NULL, maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = T) {
  reads_freq = ShortRead::tables(x, n=Inf)$top
  reads_seq = names(reads_freq)
  m = utils::aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = list(sub = 1))
  seq_l = regmatches(reads_seq, m)

  seq_v = plyr::laply(seq_l, function(x) { x[1] })
  barcode_v = plyr::laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v = plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
  if (verbose) {
    report_barcode_extraction(sample_name = sample_name, reads_num = sum(reads_freq), barcode_reads = sum(d$count, na.omit = T))
  }
  stats::na.omit(d)
}


extractBc.data.frame = function(x, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = T) {
  sequences = x$seq
  freq = x$freq
  x = DNAStringSet(rep(sequences, freq))

  extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, verbose = verbose)
}

extractBc.character = function(file, sample_name = basename(file), pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = T) {
  if (length(file) > 1) {
    bc_list = lapply(1:length(file), function(i) {
      extractBc(file[i], sample_name = sample_name[i], pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, verbose = verbose)
    })
    as.BarcodeObj(bc_list, sample_name)
  } else {
    x = ShortRead::readFastq(file)
    extractBc(x, sample_name = sample_name, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, verbose = verbose)
  }
}

extractBc.integer = function(x, pattern = "", sample_name = NULL, maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = T) {
  # TODO: Table result
  x = DNAStringSet(rep(names(x), x))

  extractBc(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, verbose = verbose)
}

extractBc.list = function(x, pattern = "", sample_name = NULL, maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), verbose = F) {
  if (is.null(sample_name)) {
    sample_name = names(x)
  }
  parallel::mclapply(1:length(x), function(i) {
    extractBc(x[[i]], pattern = pattern, sample_name = sample_name[i], maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, verbose = F)
  }) -> messyBc
  names(messyBc) = sample_name
  output = list(messyBc = messyBc)
  class(output) = "BarcodeObj"
  output
}

as.BarcodeObj = function(x, sample_name = names(x)) {
  # TODO: 
  messyBc = x
  names(messyBc) = sample_name
  output = list(messyBc = messyBc)
  class(output) = "BarcodeObj"
  output
}
