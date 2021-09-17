
bc_process_sample_name <- function(sample_name, metadata, input_names) {
    # input metadata,
    # input samle_name,
    # names of input data or the file name
    # metadata = sample_name > names or file name
    if (is.null(sample_name)) {
        if (!is.null(metadata)) {
            sample_name <- rownames(metadata)
        } else {
            sample_name <- input_names
        }
    } else {
        if (length(sample_name) != length(input_names)) {
            stop("Sample name length does not meet sample number.")
        }
    }
        
    if (!is.null(metadata)) {
        if (!all(sample_name == rownames(metadata))) {
            stop("Sample name does not match row name of metadata.")
        }
    }
    sample_name
}

bc_process_metadata <- function(sample_name, old_metadata, new_metadata) {
    res <- merge(old_metadata, new_metadata, by = 0, all=TRUE)
    res <- as.data.frame(res)
    rownames(res) <- res$Row.names
    res$Row.names <- NULL
    res[sample_name, ]
}

bc_extract_metadata <- function(x, sample_name) {
    d <- vapply(x, function(x_i) {
        res <- c(attr(x_i, "raw_read_count"), attr(x_i, "barcode_read_count"))
        res
    }, c("raw_read_count" = 1, "barcode_read_count" = 1))
    d <- as.data.frame(t(d))
    rownames(d) <- sample_name
    d
}


#' Extract barcode from sequences
#' 
#' \code{bc_extract} identifies the barcodes (and UMI) from the sequences using
#' regular expressions.  \code{pattern} and \code{pattern_type} arguments are
#' necessary, which provide the barcode (and UMI) pattern and their location
#' within the sequences.
#'
#' @param x A single or a list of fastq file, ShortReadQ, DNAStringSet,
#' data.frame, or named integer.
#' @param file A single or a vector of character string of fastq file name.
#' @param pattern A string, specifying the regular expression
#' with capture. It matchs the barcode (and UMI) with capture pattern.
#' @param sample_name A string vector, applicable when \code{x} is a list or fastq file
#' vector. This argument specifies the sample names. If not
#' provided, the function will look for sample name in the rownames of metadata,
#' the fastqfile name or the \code{list} names.
#' @param metadata A \code{data.frame} with sample names as the row names, and
#' each metadata record by column, specifying the sample characteristics. 
#' @param maxLDist A integer. The mismatch threshold for barcode matching, when
#' maxLDist is 0, the \code{\link[stringr]{str_match}}  is
#' invoked for barcode matching which is faster, otherwise
#' \code{\link[utils]{aregexec}} is invoked and the \code{costs} parameters can
#' be used to specifying the weight of the distance calculation.
#' @param pattern_type A vector. It defines the barcode (and UMI) capture
#' group. See Details.
#' @param costs A named list, applicable when maxLDist > 0, specifying the
#' weight of each mismatch events while extracting the barcodes.  The list
#' element name have to be \code{sub} (substitution), \code{ins} (insertion) and
#' \code{del} (deletion). The default value is \code{list(sub = 1, ins = 99, del
#' = 99)}.  See \code{\link[utils]{aregexec}} for more detail information.
#' @param ordered A logical value. If the value is true, the return barcodes
#' (UMI-barcode tags) are sorted by the reads counts.
#' @param ... Additional arguments
#' @details
#' The \code{pattern} argument is a regular expression, the capture operation \code{()}
#' identifying the barcode or UMI. \code{pattern_type} argument annotates 
#' capture, denoting the UMI or the barcode captured pattern. In the
#' example:
#' \preformatted{
#' ([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC
#' |---------| starts with 3 base pairs UMI.
#'            |----------| constant sequence in the backbone.
#'                        |-------| flexible barcode sequences.
#'                                 |---------| 3' constant sequence.
#' }
#'
#' In UMI part \code{[ACGT]{3}}, \code{[ACGT]} means it can be one of
#' the "A", "C", "G" and "T", and \code{{3}} means it repeats 3 times. 
#' In the barcode pattern \code{[ACGT]+}, the \code{+} denotes
#' that there is at least one of the \code{A} or \code{C} or \code{G} or \code{T.}
#' 
#' @return 
#' This function returns a BarcodeObj object if the input is a \code{list} or a
#' \code{vector} of Fastq files, otherwise it returns a \code{data.frame.} In the later case
#' the \code{data.frame} has 5 columns:
#' \enumerate{
#'   \item \code{reads_seq}: full sequence.
#'   \item \code{match_seq}: part of the full sequence matched by pattern.
#'   \item \code{umi_seq} (optional): UMI sequence, applicable when there is UMI in
#'      `pattern` and `pattern_type` argument.
#'   \item \code{barcode_seq}: barcode sequence.
#'   \item \code{count}: reads number.
#' }
#' 
#' The \code{match_seq} is part of \code{reads_seq}; The \code{umi_seq} and
#' \code{barcode_seq} are part of \code{match_seq}. The \code{reads_seq} is the
#' full sequence, and is unique id for each record (row), On the contrast,
#' \code{match_seq}, \code{umi_seq} or \code{barcode_seq} may duplicated between
#' rows.
#'
#' @examples
#' fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
#'
#' library(ShortRead)
#'
#' # barcode from fastq file
#' bc_extract(fq_file, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from ShortReadQ object
#' sr <- readFastq(fq_file)  # ShortReadQ
#' bc_extract(sr, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from DNAStringSet object
#' ds <- sr@sread  # DNAStringSet
#' bc_extract(ds, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from integer vector
#' iv <- tables(ds, n = Inf)$top # integer vector
#' bc_extract(iv, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from data.frame 
#' df <- data.frame(seq = names(iv), freq = as.integer(iv)) # data.frame
#' bc_extract(df, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from list of DNAStringSet
#' l <- list(sample1 = ds, sample2 = ds) # list
#' bc_extract(l, pattern = "AAAAA(.*)CCCCC")
#'
#' # Extract UMI and barcode
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'     )
#'   ) 
#' # barcode backbone with UMI and barcode
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_extract(
#'     list(test = d1), 
#'     pattern, 
#'     sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' ###
#' @export
bc_extract <- function(...) UseMethod("bc_extract")


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


    d <- stats::na.omit(d)

    # prepare the row_read_count & barcode_read_count
    attr(d, "raw_read_count") <- sum(reads_freq)
    attr(d, "barcode_read_count") <- sum(d$count)

    d
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
    sample_name = NULL,
    metadata = NULL, 
    maxLDist = 0, 
    pattern_type = c(barcode = 1), 
    costs = list(sub = 1, ins = 99, del = 99), 
    ordered = TRUE, ...) {

    # sample_name given
    # meta_data given
    #   no sample_name
    #   no rowname
    
    if (!is.null(metadata) & !is(metadata, "data.frame")) {
        stop("metadata should be a data.frame.")
    }

    # if more than one fastq file as input
    if (length(file) > 1) {

        input_names <- basename(file)
        sample_name <- bc_process_sample_name(sample_name, metadata, input_names)

        messyBc <- lapply(seq_along(file), function(i) {
            bc_extract(
                ShortRead::readFastq(file[i]), 
                sample_name = sample_name[i], 
                pattern = pattern, 
                maxLDist = maxLDist, 
                pattern_type = pattern_type, 
                costs = costs, 
                ordered = ordered)
        })

        new_metadata <- bc_extract_metadata(messyBc, sample_name)
        metadata <- bc_process_metadata(sample_name, metadata, new_metadata)

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

    input_names <- ifnullelse(names(x), seq_along(x))
    sample_name <- bc_process_sample_name(sample_name, metadata, input_names)

    lapply(seq_along(x), function(i) {
        bc_extract(
            x[[i]],
            pattern = pattern,
            sample_name = sample_name[i],
            maxLDist = maxLDist,
            pattern_type = pattern_type,
            costs = costs,
            ordered = ordered)
    }) -> messyBc

    new_metadata <- bc_extract_metadata(messyBc, sample_name)
    metadata <- bc_process_metadata(sample_name, metadata, new_metadata)

    names(messyBc) <- sample_name
    output <- list(messyBc = messyBc, metadata = metadata)


    class(output) <- "BarcodeObj"
    output
}

