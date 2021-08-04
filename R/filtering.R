# TODO: Trimming the sequence file remove adaptor

#' Filter the sequences by quality
#' 
#' Apply quality filters to remove low quality sequence.
#'
#' @param x A single or a list of Fastq file, ShortReadQ, DNAStringSet,
#' data.frame, integer vector.
#' @param min_average_quality A single or a vector of numeric, indicating the
#' threshold of the minimum average base quality of a sequence. 
#' @param min_read_length A single or a vector of integer, specifying the length
#' threshold of a sequence. 
#' @param N_threshold A single or a vector of integer, specifying the maximum
#' \code{N} in a sequence.
#' @param sample_name A string vector, rename the samples.
#' @param ... Additional arguments
#' @return A ShortReadQ or DNAStringSet object with sequences passed the filters
#' @examples
#' library(ShortRead)
#' 
#' fq_file <- system.file("extdata", "simple.fq", package="Bc")
#'
#' # apply filter to fastq files
#' bc_filterSeq(fq_file)
#'
#' # Read fastq files get ShortReadQ object
#' sr <- readFastq(fq_file[1])
#' # apply filter to ShortReadQ
#' bc_filterSeq(sr)
#'
#' # get DNAStringSet object
#' ds <- sr@sread
#' # apply filter to DNAStringSet
#' bc_filterSeq(ds)
#'
#' ###
#' @export
bc_filterSeq <- function(x, ...) UseMethod("bc_filterSeq")

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.ShortReadQ <- function(
    x, 
    min_average_quality = 30, 
    min_read_length = 0,
    N_threshold = 0, ...) {

    # good base quality filter
    goodq <- ShortRead::srFilter(function(x) { 
        apply(methods::as(Biostrings::quality(x), "matrix"), 1, mean, na.rm=TRUE) >= 
            min_average_quality 
    }, name="GoodQualityBases")

    # good sequences length filter
    goodlength <- ShortRead::srFilter(function(x) { width(ShortRead::sread(x)) >= 
        min_read_length}, name="GoodReadLength")

    # 'N' nucleic filter
    goodN <- ShortRead::nFilter(threshold = N_threshold)

    goodFinal <- ShortRead::compose(goodlength, goodq, goodN)
    x[goodFinal(x)]
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.DNAStringSet <- function(
    x,
    min_read_length = 0,
    N_threshold = 0, ...) {

    goodlength <- ShortRead::srFilter(function(x) {
        width(x) >= min_read_length
    }, name="GoodReadLength")

    goodN <- ShortRead::nFilter(threshold = N_threshold)
    goodFinal <- ShortRead::compose(goodlength, goodN)
    x[goodFinal(x)]
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.data.frame <- function(x,
    min_read_length = 0,
    N_threshold = 0, ...) {

    sequences <- x$seq
    freq <- x$freq
    x <- Biostrings::DNAStringSet(rep(sequences, freq))

    bc_filterSeq(x, min_read_length = min_read_length, N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.character <- function(x,
    min_average_quality = 30,
    min_read_length = 0,
    N_threshold = 0,
    sample_name = basename(x), ...) {

    # TODO: when sample_name length is not right ... 
    fq_list <- lapply(x, ShortRead::readFastq)
    names(fq_list) <- sample_name

    bc_filterSeq(
        fq_list,
        min_average_quality = min_average_quality,
        min_read_length = min_read_length,
        N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.integer <- function(x, min_read_length = 0, N_threshold = 0, ...) {

    # TODO: Use Rle?
    x <- Biostrings::DNAStringSet(rep(names(x), x))
    bc_filterSeq(x, min_read_length = min_read_length, N_threshold = N_threshold)
}

#' @rdname bc_filterSeq
#' @export
bc_filterSeq.list <- function(
    x, 
    min_average_quality = 30,
    min_read_length = 0,
    N_threshold = 0,
    sample_name = names(x), ...) {


    if (is(x[[1]], "ShortReadQ") | is(x[[1]], "character")) {
        output <- lapply(
            x,
            bc_filterSeq,
            min_average_quality = min_average_quality,
            min_read_length = min_read_length,
            N_threshold = N_threshold)

        names(output) <- sample_name
        output
    } else {
        output <- lapply(
            x,
            bc_filterSeq,
            min_read_length = min_read_length,
            N_threshold = N_threshold)

        names(output) <- sample_name
        output
    }
}
