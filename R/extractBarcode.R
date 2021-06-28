#' Extract barcode from reads
#' 
#' This function extracts the barcodes and (or) UMI from the seqeunces.
#' `pattern` and `pattern_type` argument is necessary, which provide the barcode
#' and UMI pattern and location within the sequences.
#'
#' @param x file location, ShortReadQ, DNAStringSet, data.frame, integer or list
#' of above data types. 
#' @param file A character vector, Fastq file name.
#' @param pattern A string. The regular expression with capture to match the
#' barcode which capture pattern
#' @param sample_name A string vector. When x is list or file directions, this
#' parameter provides the sample names. If the sample_name is not given, then
#' the list names of the x will be used, or the file base name will be used as
#' sample name.
#' @param metadata A data.frame with rownames of sample_name, the column is the
#' meta data to keep the characteristics of the samples.
#' @param maxLDist A integer. The mismatch threshold for barcode matching
#' @param pattern_type A vector. It defines the barcode (and UMI) capture
#' pattern
#' @param costs A list. Define the weight for each mismatch events while
#' matching the barcode using pattern, with names of sub(substitution),
#' ins(insertion), del(deletion). The default value is list(sub = 1, ins = 99,
#' del = 99).
#' @param ordered A bool value. if the value is true, then the return barcodes
#' (tages) are ordered by the reads counts.
#' @param ... Additional arguments
#' @details
#' The `pattern` argument is a regular expression, the capture operation `()`
#' shows the barcode or UMI. `pattern_type` argument is the annotation of the
#' capture pattern, shows which captured patteern is UMI or the barcode. In the
#' example:
#' ([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC
#' |--- starts with 3 base pairs UMI
#'            |--- constant sequence in the backbone
#'                        |--- flexible barcode seqeunces.
#'                                |--- 3' constant sequence
#' 
#' In the UMI pattern `[ACGT]{3}`, the `[ACGT]` means it can be one of the "A",
#' "C", "G" and "T", and the `{3}` means there are 12 `[ACGT]`. In the barcode
#' pattern `[ACGT]+`, the `+` denotes that there is at least one of the
#' `[ACGT]`.
#' 
#' @return 
#' This function returns a `BarcodeObj` object if the input is a list or a
#' vector of Fastq files, otherwise it returns a `data.frame`. In the later case
#' the `data.frame` has 5 columns:
#'   1. reads_seq: sequence of the reads before parsing.
#'   2. match_seq: the sequence among read matched by pattern.
#'   3. umi_seq (optional): UMI sequence, applicable when there is UMI in
#'      `pattern` and `pattern_type` argument.
#'   4. barcode_seq: barcode sequence.
#'   5. count: reads number.
#' 
#' @details In the output, the `match_seq` is part of `reads_seq`. The `umi_seq`
#' and `barcode_seq` are part of `match_seq`. Be attention that, the `reads_seq`
#' is the unique id for each row. The `match_seq`, `umi_seq` or `barcode_seq`
#' can be duplicated, due to the potential variation in the region outside of
#' `match_seq`. Please keep this in mind when you use data in `$messyBc`.
#'
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
#' # barcode backbone with UMI and barcode
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_extract(
#'   list(test = d1), 
#'   pattern, 
#'   sample_name=c("test"), 
#'   pattern_type=c(UMI=1, barcode=2))
#'
#' ###
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
bc_extract.data.frame <- function(
  x, 
  pattern = "", 
  sample_name = NULL, 
  maxLDist = 0, 
  pattern_type = c(barcode = 1), 
  costs = list(sub = 1, ins = 99, del = 99), 
  ordered = TRUE, ...) {

  # sequence frequecy
  reads_freq <- x$freq
  # sequence
  reads_seq <- x$seq

  # if maxLDist > 0, use utils::aregexec to match the barcode else use
  # stringr::str_match to perform the barcode match. The stringr::str_match is
  # faster
  if (maxLDist > 0) {
    # regular expression match
    m <- utils::aregexec(
      pattern, 
      reads_seq, 
      fixed = FALSE, 
      max.distance = maxLDist, 
      costs = costs
    )

    # captured sequence
    seq_l <- regmatches(reads_seq, m)

    # matched sequence
    seq_v <- plyr::laply(seq_l, function(x) { x[1] })
    # captured barcode
    barcode_v <- plyr::laply(seq_l, function(x) { 
      x[pattern_type["barcode"] + 1] 
    })

    if ("UMI" %in% names(pattern_type)) {
      umi_v <- plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
    }

  } else {
    # regular expression match
    m <- stringr::str_match(reads_seq, pattern)

    # capture sequence
    seq_v <- m[, 1]

    # matched barcode
    barcode_v <- m[, pattern_type['barcode'] + 1]

    # matched UMI
    if ("UMI" %in% names(pattern_type)) {
      umi_v <- m[, pattern_type['UMI'] + 1]
    }
  }

  # captured UMI
  if ("UMI" %in% names(pattern_type)) {
    d <- data.table(
      reads_seq = reads_seq, 
      match_seq = seq_v, 
      umi_seq = umi_v, 
      barcode_seq = barcode_v, 
      count = reads_freq
    )
  } else {
    d <- data.table(
      reads_seq = reads_seq, 
      match_seq = seq_v, 
      barcode_seq = barcode_v, 
      count = reads_freq
    )
  }

  # order the data by counts
  if (ordered) {
    d <- d[order(count, decreasing = TRUE)]
  } 
  stats::na.omit(d)
}


#' @rdname bc_extract
#' @export
bc_extract.ShortReadQ <- function(
  x, 
  pattern = "", 
  sample_name = NULL, 
  maxLDist = 0, 
  pattern_type = c(barcode = 1), 
  costs = list(sub = 1, ins = 99, del = 99), 
  ordered = TRUE, ...) {

  # sequence frequecy
  reads_freq <- ShortRead::tables(x, n=Inf)$top

  x = data.frame(
    freq = as.integer(reads_freq),
    seq = names(reads_freq)
  )

  bc_extract(x, 
    pattern = pattern, 
    maxLDist = maxLDist, 
    pattern_type = pattern_type, 
    costs = costs, 
    ordered = ordered
  )
}


#' @rdname bc_extract
#' @export
bc_extract.DNAStringSet <- function(
  x, 
  pattern = "", 
  sample_name = NULL, 
  maxLDist = 0, 
  pattern_type = c(barcode = 1), 
  costs = list(sub = 1, ins = 99, del = 99), 
  ordered = TRUE, ...) {

  bc_extract.ShortReadQ(
    x, 
    pattern = pattern, 
    sample_name = sample_name, 
    maxLDist = maxLDist, 
    pattern_type = pattern_type, 
    costs = costs, 
    ordered = ordered
  )
}

#' @rdname bc_extract
#' @export
bc_extract.integer <- function(
  x,
  pattern = "",
  sample_name = NULL,
  maxLDist = 0,
  pattern_type = c(barcode = 1),
  costs = list(sub = 1, ins = 99, del = 99),
  ordered = TRUE, ...) {

  x = data.frame(
    freq = as.integer(x),
    seq = names(x)
  )

  bc_extract(
    x,
    pattern = pattern,
    maxLDist = maxLDist,
    pattern_type = pattern_type,
    costs = costs,
    ordered = ordered)
}


#' @rdname bc_extract
#' @export
bc_extract.character <- function(file, 
  pattern = "", 
  sample_name = basename(file), 
  metadata = NULL, 
  maxLDist = 0, 
  pattern_type = c(barcode = 1), 
  costs = list(sub = 1, ins = 99, del = 99), 
  ordered = TRUE, ...) {

  # sample_name given
  # meta_data given
  #   no sample_name
  #   no rowname

  # if more than one fastq file as input
  if (length(file) > 1) {

    # use sample_name
    if (is.null(sample_name)) {
      if (!is.null(metadata$sample_name)) {
        # use metadata
        sample_name = metadata$sample_name
      } else {
        # use file name
        sample_name <- names(file)
      }
    }
    # still no sample_name use length of file
    if (is.null(sample_name)) {
      sample_name <- 1:length(file)
    }

    if (is.null(metadata)) {
      metadata <- data.frame(sample_name = sample_name)
    }
    if (is.null(metadata$sample_name)) {
      metadata$sample_name <- sample_name
    }
    rownames(metadata) <- sample_name

    if (length(metadata$sample_name) != length(file))
      stop("sample_name or metadata does not match sample number.")


    messyBc <- lapply(1:length(file), function(i) {
      bc_extract(
        ShortRead::readFastq(file[i]), 
        sample_name = sample_name[i], 
        pattern = pattern, 
        maxLDist = maxLDist, 
        pattern_type = pattern_type, 
        costs = costs, 
        ordered = ordered)
    })

    names(messyBc) <- sample_name
    output <- list(messyBc = messyBc, metadata = metadata)
    class(output) <- "BarcodeObj"
    return(output)
  } else {
    # if one fastq file as input
    bc_extract(
      ShortRead::readFastq(file), 
      sample_name = sample_name[i], 
      pattern = pattern,
      maxLDist = maxLDist,
      pattern_type = pattern_type,
      costs = costs,
      ordered = ordered)
  }
}


#' @rdname bc_extract
#' @export
bc_extract.list <- function(
  x,
  pattern = "",
  sample_name = NULL,
  metadata = NULL,
  maxLDist = 0,
  pattern_type = c(barcode = 1),
  costs = list(sub = 1, ins = 99, del = 99),
  ordered = TRUE, ...) {

  # use sample_name
  if (is.null(sample_name)) {
    if (!is.null(metadata)) {
      # use metadata
      sample_name = metadata$sample_name
    } else {
      # use list name
      sample_name <- names(x)
    }
  }

  # still no sample_name use length of list
  if (is.null(sample_name)) {
    sample_name <- 1:length(x)
  }

  if (is.null(metadata)) {
    metadata <- data.frame(sample_name = sample_name)
  }
  if (is.null(metadata$sample_name)) {
    metadata$sample_name <- sample_name
  }
  rownames(metadata) <- sample_name


  if (length(metadata$sample_name) != length(x))
    stop("sample_name or metadata does not match sample number.")


  lapply(1:length(x), function(i) {
    bc_extract(
      x[[i]],
      pattern = pattern,
      sample_name = sample_name[i],
      maxLDist = maxLDist,
      pattern_type = pattern_type,
      costs = costs,
      ordered = ordered)
  }) -> messyBc

  names(messyBc) <- sample_name
  output <- list(messyBc = messyBc, metadata = metadata)
  class(output) <- "BarcodeObj"
  output
}

