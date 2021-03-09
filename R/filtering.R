#' Filter the reads by the QC threshold
filterFq = function(...) UseMethod("filterFq")

#' Filter the reads by the QC threshold
filterFq.ShortReadQ = function(x, min_average_quality = 30, min_read_length = 0) {
  goodq <- srFilter(function(x) { apply(methods::as(quality(x), "matrix"), 1, mean, na.rm=TRUE) >= min_average_quality }, name="GoodQualityBases")
  goodlength <- srFilter(function(x) { width(ShortRead::sread(x)) >= min_read_length}, name="GoodReadLength")
  x = x[goodq(x)]
  x[goodlength(x)]
}

filterFq.DNAStringSet = function(x, min_read_length = 0) {
  goodlength <- srFilter(function(x) { width(x) >= min_read_length}, name="GoodReadLength")
  x[goodlength(x)]
}

filterFq.data.frame = function(x, min_read_length = 0) {
  x[nchar(x) >= min_read_length, ]
}

filterFq.character = function(x, min_read_length = 0) {
  x[nchar(x) >= min_read_length]
}

filterFq.integer = function(x, min_read_length = 0) {
  x[nchar(names(x)) >= min_read_length]
}

filterFq.list = function(x, min_average_quality = 30, min_read_length = 0) {
  if (class(x[0]) == "ShortReadQ") {
    lapply(x, filterFq, min_average_quality = min_average_quality, min_read_length = min_read_length)
  } else {
    lapply(x, filterFq, min_read_length = min_read_length)
  }
}
