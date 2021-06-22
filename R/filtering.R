# TODO: Trimming the sequence file remove adaptor

#' Filter the reads by the QC threshold
#'
#' @param x file location, ShortReadQ, DNAStringSet, data.frame, integer vector or list of above data type.
#' @param min_average_quality A number, sequence with average base smaller than the value are filtered.
#' @param min_read_length A integer, sequence with length shorter than the value are filtered.
#' @param N_threshold A integer, sequence with reads number less that the value are filtered.
#' @param sample_name A string vector, applicable when a list is used as input, it provides new sample_name for the output list.
#' @return Sequences pass the filter criterial
#' @examples
#' fq_file = system.file("extdata", "simple.fq", package="Bc")
#' sr = readFastq(fq_file)
#' ds = sr@sread
#'
#' bc_filterSeq(fq_file)
#' bc_filterSeq(sr)
#' bc_filterSeq(ds)
#' @export
bc_filterSeq = function(...) UseMethod("bc_filterSeq")

#' Filter the reads by the QC threshold for ShortReadQ
#'
#' @export
bc_filterSeq.ShortReadQ = function(x, min_average_quality = 30, min_read_length = 0, N_threshold = 0) {
  goodq <- srFilter(function(x) { apply(methods::as(quality(x), "matrix"), 1, mean, na.rm=TRUE) >= min_average_quality }, name="GoodQualityBases")
  goodlength <- srFilter(function(x) { width(ShortRead::sread(x)) >= min_read_length}, name="GoodReadLength")
  goodN <- nFilter(threshold = N_threshold)
  goodFinal = compose(goodlength, goodq, goodN)
  x[goodFinal(x)]
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.DNAStringSet = function(x, min_read_length = 0, N_threshold = 0) {
  goodlength <- srFilter(function(x) { width(x) >= min_read_length}, name="GoodReadLength")
  goodN <- nFilter(threshold = N_threshold)
  goodFinal = compose(goodlength, goodN)
  x[goodFinal(x)]
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.data.frame = function(x, min_read_length = 0, N_threshold = 0) {
  sequences = x$seq
  freq = x$freq
  x = DNAStringSet(rep(sequences, freq))

  bc_filterSeq(x, min_read_length = min_read_length, N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.character = function(x, min_average_quality = 30, min_read_length = 0, N_threshold = 0, sample_name = basename(x)) {
  # TODO: when sample_name length is not right ... 
  fq_list = lapply(x, ShortRead::readFastq)
  names(fq_list) = sample_name
  bc_filterSeq(fq_list, min_average_quality = min_average_quality, min_read_length = min_read_length, N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.integer = function(x, min_read_length = 0, N_threshold = 0) {
  # TODO: Use Rle?
  x = DNAStringSet(rep(names(x), x))
  bc_filterSeq(x, min_read_length = min_read_length, N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.list = function(x, min_average_quality = 30, min_read_length = 0, N_threshold = 0, sample_name = names(x)) {
  if (is(x[[1]], "ShortReadQ")) {
    output = lapply(x, bc_filterSeq, min_average_quality = min_average_quality, min_read_length = min_read_length, N_threshold = N_threshold)
    names(output) = sample_name
    output
  } else {
    output = lapply(x, bc_filterSeq, min_read_length = min_read_length, N_threshold = N_threshold)
    names(output) = sample_name
    output
  }
}
