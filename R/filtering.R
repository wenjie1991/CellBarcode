#' Filter the reads by the QC threshold
#'
#' @export
bc_filterSeq = function(...) UseMethod("bc_filterSeq")

#' Filter the reads by the QC threshold
#'
#' @export
bc_filterSeq.ShortReadQ = function(x, min_average_quality = 30, min_read_length = 0) {
  goodq <- srFilter(function(x) { apply(methods::as(quality(x), "matrix"), 1, mean, na.rm=TRUE) >= min_average_quality }, name="GoodQualityBases")
  goodlength <- srFilter(function(x) { width(ShortRead::sread(x)) >= min_read_length}, name="GoodReadLength")
  x = x[goodq(x)]
  x[goodlength(x)]
}

#' @export
bc_filterSeq.DNAStringSet = function(x, min_read_length = 0) {
  goodlength <- srFilter(function(x) { width(x) >= min_read_length}, name="GoodReadLength")
  x[goodlength(x)]
}

#' @export
bc_filterSeq.data.frame = function(x, min_read_length = 0) {
  x[nchar(x) >= min_read_length, ]
}

#' @export
bc_filterSeq.character = function(x, min_average_quality = 30, min_read_length = 0, sample_name = basename(x)) {
  # TODO: when sample_name length is not right ... 
  fq_list = lapply(x, ShortRead::readFastq)
  names(fq_list) = sample_name
  bc_filterSeq(fq_list, min_average_quality = min_average_quality, min_read_length = min_read_length)
}

#' @export
bc_filterSeq.integer = function(x, min_read_length = 0) {
  x[nchar(names(x)) >= min_read_length]
}

#' @export
bc_filterSeq.list = function(x, min_average_quality = 30, min_read_length = 0, sample_name = names(x)) {
  if (class(x[0]) == "ShortReadQ") {
    output = lapply(x, bc_filterSeq, min_average_quality = min_average_quality, min_read_length = min_read_length)
    names(output) = sample_name
    output
  } else {
    output = lapply(x, bc_filterSeq, min_read_length = min_read_length)
    names(output) = sample_name
    output
  }
}
