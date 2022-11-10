#' Parse 10X Genomic scRNASeq sam file (Experimental)
#'
#' \code{bc_extract_10XscSeq} can extract cellular barcode, UMI and lineage barcode
#' sequences from 10X Genomics scRNASeq bam file. This function can not process
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
#' bc_extract_10XscSeq(
#'   sam = sam_file,
#'   pattern = "AGATCAG(.*)TGTGGTA",
#'   cell_barcode_tag = "CR",
#'   umi_tag = "UR"
#' )
#'
#' @rdname bc_extract_10X_scSeq
#' @export
bc_extract_10XscSeq <- function(
    sam, 
    pattern, 
    cell_barcode_tag = "CR", 
    umi_tag = "UR"
    ) {
    if(!file.exists(sam)) {
        stop("The input bam file is not exist.")
    }
    l <- parse_10x_scSeq(sam, pattern)
    d <- l$barcode_df
    data.table::setDT(d)

    raw_count_dt <- l$raw_reads_df
    data.table::setDT(raw_count_dt)
    setkey(raw_count_dt, "cell_barcode")

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
        }
    )

    attr(messyBc, "split_type") <- NULL
    attr(messyBc, "split_labels") <- NULL

    barcode_read_count_dt <- d[, .(count = sum(count)), by = cell_barcode]
    setkey(barcode_read_count_dt, "cell_barcode")

    metadata <- data.frame(
        raw_read_count = raw_count_dt[names(messyBc), count],
        barcode_read_count = barcode_read_count_dt[names(messyBc), count]
    )
    rownames(metadata) <- names(messyBc)

    output <- BarcodeObj(metadata = metadata, messyBc = messyBc)

    output
}
