process_sc_list <- function(l) {
    d <- l$barcode_df
    data.table::setDT(d)

    raw_count_dt <- l$raw_reads_df
    data.table::setDT(raw_count_dt)
    data.table::setkey(raw_count_dt, "cell_barcode")

    d <- d[, .(count = sum(count)), by = .(cell_barcode, umi, barcode_seq)][
        , .(cell_barcode, umi_seq = umi, barcode_seq, count)
        ][order(count, decreasing = TRUE)]

    # as.data.frame(d)
    messyBc <- plyr::dlply(d, .(cell_barcode), 
        function(x) {
            # attr(x, "raw_read_count") <- raw_count_dt[x$cell_barcode[1], count]
            # attr(x, "barcode_read_count") <- sum(x$count)
            x$cell_barcode <- NULL
            data.table(x)
        })

    attr(messyBc, "split_type") <- NULL
    attr(messyBc, "split_labels") <- NULL

    barcode_read_count_dt <- d[, .(count = sum(count)), by = cell_barcode]
    data.table::setkey(barcode_read_count_dt, "cell_barcode")

    metadata <- data.frame(
        raw_read_count = raw_count_dt[names(messyBc), count],
        barcode_read_count = barcode_read_count_dt[names(messyBc), count]
    )
    rownames(metadata) <- names(messyBc)

    output <- BarcodeObj(metadata = metadata, messyBc = messyBc)

    output
}

#' Parse 10X Genomic scRNASeq sam file (Experimental)
#'
#' \code{bc_extract_10X_sam} can extract cellular barcode, UMI and lineage barcode
#' sequences from 10X Genomics scRNASeq sam file. This function can not process
#' bam file directly, user needs to uncompress the bam file to get sam file in
#' order to run this function See example.
#'
#' @param sam A string, define the un-mapped sequences 
#' @param pattern A string, define the regular expression to match the barcode
#' sequence. The barcode sequence should be in the first catch. Please see the
#' documents of \code{\link[CellBarcode]{bc_extract}} and example for more information.
#' @param cell_barcode_tag A string, define the tag of 10X cell barcode field in sam
#' file. The default is "CR".
#' @param umi_tag A string, define the tag of UMI field in the sam file.
#' @details 
#' If the barcode sequence does not map to the reference genome. The user should
#' use samtools to get the un-mapped reads and save it as sam format for using
#' as the input. It can save a lot of time. The way to get the un-mapped reads:
#' \preformatted{
#' samtools view -f 4 input.bam > output.sam 
#' }
#' 
#' @return 
#' A data.frame with 4 columns:
#' \enumerate{
#'   \item \code{cell_barcode}: 10X cellular barcode.
#'   \item \code{umi}: UMI sequence.
#'   \item \code{barcode_seq}: lineage barcode.
#'   \item \code{count}: reads count.
#'   }
#' @examples
#' ## NOT run
#' # In the case that when the barcode sequence is not mapped to 
#' # reference genome, it will be much more efficiency to get 
#' # the un-mapped sequences as the input.
#'
#' ## Get un-mapped reads
#' # samtools view -f 4 input.bam > scRNASeq_10X.sam 
#'
#' sam_file <- system.file("extdata", "scRNASeq_10X.sam", package = "CellBarcode")
#'
#' bc_extract_10X_sam(
#'   sam = sam_file,
#'   pattern = "AGATCAG(.*)TGTGGTA",
#'   cell_barcode_tag = "CR",
#'   umi_tag = "UR"
#' )
#'
#' @rdname bc_extract_10X_sam
#' @export
bc_extract_10X_sam <- function(
    sam, 
    pattern, 
    cell_barcode_tag = "CR", 
    umi_tag = "UR"
    ) {
    if(!file.exists(sam)) {
        stop("The input bam file does not exist.")
    }
    sam <- path.expand(sam)
    l <- parse_10x_sam(sam, pattern)

    process_sc_list(l)
}

#' @export
bc_extract_10X_fastq <- function(
    fq1,
    fq2 = NULL,
    patternCellBarcode = NULL,
    patternUMI = NULL,
    patternBarcode = NULL,
    isOverlap = F
    ) {

    if (!file.exists(fq1)) {
        stop("The input fastq files do not exist.")
    } else {
        fq1 <- path.expand(fq1)
    }

    if (is.null(fq2) & file.exists(fq2)) {
        stop("The input fastq files do not exist.")
    } else {
        fq2 <- path.expand(fq2)
    }

    if (is.null(fq2)) {
        # TODO: l <- parse_10x_fastq(fq1, pattern)
        # TODO: analysis fq file
        d <- read_fastq_gz(fq1)
    } else {
        # join fq1 and fq2 and parse
        # TODO: overlap is true
        # TODO: analysis_fq file
        d <- read_fastq_gz2(fq1, fq2)
    }

    cb_x <- stringr::str_match(d$seq, patternCellBarcode)[, 2]
    umi_x <- stringr::str_match(d$seq, patternUMI)[, 2]
    barcode_x <- stringr::str_match(d$seq, patternBarcode)[, 2]

    d_raw <- data.table(cell_barcode = cb_x, umi = umi_x, barcode_seq  = barcode_x, count = d$freq)
    raw_count_dt <- d_raw[, .(count = sum(count)), by = cell_barcode]
    data.table::setkey(raw_count_dt, "cell_barcode")

    barcode_read_count_dt <- na.omit(d_raw)[, .(count = sum(count)), by = cell_barcode]
    data.table::setkey(barcode_read_count_dt, "cell_barcode")

    d <- na.omit(d_raw)[, .(count = sum(count)), by = .(cell_barcode, umi, barcode_seq)][
        , .(cell_barcode, umi_seq = umi, barcode_seq, count)
        ][order(count, decreasing = TRUE)]

    # as.data.frame(d)
    messyBc <- plyr::dlply(d, .(cell_barcode), 
        function(x) {
            # attr(x, "raw_read_count") <- raw_count_dt[x$cell_barcode[1], count]
            # attr(x, "barcode_read_count") <- sum(x$count)
            x$cell_barcode <- NULL
            data.table(x)
        })

    attr(messyBc, "split_type") <- NULL
    attr(messyBc, "split_labels") <- NULL

    barcode_read_count_dt <- d[, .(count = sum(count)), by = cell_barcode]
    data.table::setkey(barcode_read_count_dt, "cell_barcode")

    metadata <- data.frame(
        raw_read_count = raw_count_dt[names(messyBc), count],
        barcode_read_count = barcode_read_count_dt[names(messyBc), count]
    )
    rownames(metadata) <- names(messyBc)

    output <- BarcodeObj(metadata = metadata, messyBc = messyBc)

    output
}
