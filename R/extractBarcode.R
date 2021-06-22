#' Extract barcode from reads
#'
#' @param x file location, ShortReadQ, DNAStringSet, data.frame, integer or list of above data types. 
#' @param file A character vector, Fastq file name.
#' @param pattern A string. The regular expression with capture to match the barcode which capture pattern
#' @param sample_name A string vector. When x is list or file directions, this parameter provides the sample names. If the sample_name is not given, then the list names of the x will be used, or the file base name will be used as sample name.
#' @param meta_data A data.frame with rownames of sample_name, the column is the meta data to keep the characteristics of the samples.
#' @param maxLDist A integer. The mismatch threshold for barcode matching
#' @param pattern_type A vector. It defines the barcode (and UMI) capture pattern
#' @param costs A list. Define the weight for each mismatch events while matching the barcode using pattern, with names of sub(substitution), ins(insertion), del(deletion). The default value is list(sub = 1, ins = 99, del = 99).
#' @param ordered A bool value. if the value is true, then the return barcodes (tages) are ordered by the reads counts.
#' @param ... Additional arguments
#' @return If x is a file location vector with more than 1 items or x is a list, the output will be a BarcodeObj. Otherwise a data.frame will be returned.
#' @examples
#' fq_file <- system.file("extdata", "simple.fq", package="Bc")
#'
#' library(ShortRead)
#'
#' sr <- readFastq(fq_file)  # ShortReadQ
#' ds <- sr@sread  # DNAStringSet
#' iv <- tables(ds, n = Inf)$top # integer vector
#' df <- data.frame(seq = names(iv), freq = as.character(iv)) # data.frame
#' l <- list(sample1 = ds, sample2 = ds) # list
#' 
#' ds
#' # extract barcode from fastq file
#' bc_extract(fq_file, pattern = "AAAAA(.*)CCCCC")
#'
#' # extract barcode from ShortReadQ object
#' bc_extract(sr, pattern = "AAAAA(.*)CCCCC")
#'
#' # extract barcode from DNAStringSet object
#' bc_extract(ds, pattern = "AAAAA(.*)CCCCC")
#'
#' # extract barcode from integer vector
#' bc_extract(iv, pattern = "AAAAA(.*)CCCCC")
#'
#' # extract barcode from data.frame 
#' bc_extract(df, pattern = "AAAAA(.*)CCCCC")
#'
#' # extract barcode from list of DNAStringSet
#' bc_extract(l, pattern = "AAAAA(.*)CCCCC")
#'
#' # Extract UMI and barcode
#' d1 <- data.frame(
#'   seq = c(
#'     "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'     "AATTCGATCGATCGAAGAGATCGATCGATC",
#'     "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'     "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'     "AAATCGATCGATCGAAGAGATCGATCGATC",
#'     "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'     "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'     "GGATCGATCGATCGAAGAGATCGATCGATC",
#'     "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'     "GGTTCGATCGATCGACGAGATCGATCGATC",
#'     "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'     ),
#'   freq = c(
#'     30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'     )
#'   ) 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' (bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'   pattern_type=c(UMI=1, barcode=2)))
#' (bc_obj <- bc_cure(bc_obj, depth=0, doFish=FALSE, with_umi=TRUE, umi_depth=5, isUniqueUMI=TRUE))
#' @export
bc_extract <- function(...) UseMethod("bc_extract")

# report_barcode_extraction <- function(sample_name, reads_num, barcode_reads) {
#   if (is.null(sample_name)) {
#     cat("In this sample, ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads ,".\n", sep = "")
#   } else {
#     cat("In sample ", sample_name, ", ", barcode_reads, "(", round(barcode_reads / reads_num * 100, 2),"%) reads with barcode are identified from total reads of ", reads_num ,".\n", sep = "")
#   }
# }

#' @rdname bc_extract
#' @export
bc_extract.ShortReadQ <- function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  reads_freq <- ShortRead::tables(x, n=Inf)$top
  reads_seq <- names(reads_freq)
  m <- utils::aregexec(pattern, reads_seq, fixed = FALSE, max.distance = maxLDist, costs = costs)
  seq_l <- regmatches(reads_seq, m)

  seq_v <- plyr::laply(seq_l, function(x) { x[1] })
  barcode_v <- plyr::laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
  if ("UMI" %in% names(pattern_type)) {
    umi_v <- plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    d <- data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
  } else {
    d <- data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
  }
  if (ordered) {
    d <- d[order(count, decreasing = TRUE)]
  } 
  stats::na.omit(d)
}

#' @rdname bc_extract
#' @export
bc_extract.DNAStringSet <- function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  bc_extract.ShortReadQ(x, pattern = pattern, sample_name = sample_name, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
}

#' @rdname bc_extract
#' @export
bc_extract.data.frame <- function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  sequences <- x$seq
  freq <- x$freq
  x <- DNAStringSet(rep(sequences, freq))

  bc_extract(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
}

#' @rdname bc_extract
#' @export
bc_extract.character <- function(file, pattern = "", sample_name = basename(file), meta_data = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  if (length(file) > 1) {
    if (is.null(sample_name)) {
      sample_name <- names(file)
    }
    if (is.null(sample_name)) {
      sample_name <- 1:length(file)
    }
    if (is.null(meta_data)) {
      meta_data <- data.frame(sample_name = sample_name)
      rownames(meta_data) <- sample_name
    }

    messyBc <- lapply(1:length(file), function(i) {
      bc_extract(ShortRead::readFastq(file[i]), sample_name = sample_name[i], pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
    })

    names(messyBc) <- sample_name
    output <- list(messyBc = messyBc, meta_data = meta_data)
    class(output) <- "BarcodeObj"
    return(output)
  } else {
    bc_extract(ShortRead::readFastq(file), sample_name = sample_name[i], pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
  }
}

#' @rdname bc_extract
#' @export
bc_extract.integer <- function(x, pattern = "", sample_name = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  # TODO: Table result
  x <- DNAStringSet(rep(names(x), x))

  bc_extract(x, pattern = pattern, maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
}

#' @rdname bc_extract
#' @export
bc_extract.list <- function(x, pattern = "", sample_name = NULL, meta_data = NULL, maxLDist = 0, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99), ordered = TRUE, ...) {
  if (is.null(sample_name)) {
    sample_name <- names(x)
  }
  if (is.null(sample_name)) {
    sample_name <- 1:length(x)
  }
  if (is.null(meta_data)) {
    meta_data <- data.frame(sample_name = sample_name)
    rownames(meta_data) <- sample_name
  } else {
    meta_data$sample_name <- rownames(sample_name)
  }

  parallel::mclapply(1:length(x), function(i) {
    bc_extract(x[[i]], pattern = pattern, sample_name = sample_name[i], maxLDist = maxLDist, pattern_type = pattern_type, costs = costs, ordered = ordered)
  }) -> messyBc

  names(messyBc) <- sample_name
  output <- list(messyBc = messyBc, meta_data = meta_data)
  class(output) <- "BarcodeObj"
  output
}

