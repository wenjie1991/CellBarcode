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

#' Extract barcode from single-cell sequencing sam file
#'
#' \code{bc_extract_sc_sam} can extract cellular barcode, UMI, and lineage
#' barcode sequences from 10X Genomics scRNASeq sam file (or bam file have
#' similar data structure). This function can not process bam file directly,
#' users need to uncompress the bam file to get a sam file to run this
#' function See example.
#'
#' @param sam A string, define the un-mapped sequences 
#' @param bam A string, define the bam file, it will be converted to sam file 
#' @param pattern A string, define the regular expression to match the barcode
#' sequence. The barcode sequence should be in the first catch. Please see the
#' documents of \code{\link[CellBarcode]{bc_extract}} and example for more information.
#' @param cell_barcode_tag A string, define the tag of cellular barcode field in sam
#' file. The default is "CR".
#' @param umi_tag A string, define the tag of a UMI field in the sam file.
#' @details 
#'
#' Although the function `bc_extract_sc_bam` can process bam file directly,
#' some optimization is still working on, it will be much more efficient to use
#' `samtools` to get the sam file.
#'
#' What's more, if the barcode sequence does not map to the reference genome. The user should
#' use the samtools to get the un-mapped reads and save it as sam format for using
#' as the input. It can save a lot of time. The way to get the un-mapped reads:
#' \preformatted{
#' samtools view -f 4 input.bam > output.sam 
#' }
#' 
#' @return 
#' A BarcodeObj object with each cell as a sample.
#' @seealso \code{\link[CellBarcode]{bc_extract}},
#' \code{\link[CellBarcode]{bc_extract_sc_fastq}}
#'
#' @examples
#' ## NOT run
#' # In the case that when the barcode sequence is not mapped to 
#' # reference genome, it will be much more efficient to get 
#' # the un-mapped sequences as the input.
#'
#' ## Get un-mapped reads
#' # samtools view -f 4 input.bam > scRNASeq_10X.sam 
#'
#' sam_file <- system.file("extdata", "scRNASeq_10X.sam", package = "CellBarcode")
#'
#' bc_extract_sc_sam(
#'   sam = sam_file,
#'   pattern = "AGATCAG(.*)TGTGGTA",
#'   cell_barcode_tag = "CR",
#'   umi_tag = "UR"
#' )
#'
#' @rdname bc_extract_sc_sam
#' @export
bc_extract_sc_sam <- function(
    sam, 
    pattern, 
    cell_barcode_tag = "CR", 
    umi_tag = "UR"
    ) {
    if(!file.exists(sam)) {
        stop("The input sam file does not exist.")
    }
    sam <- path.expand(sam)
    l <- parse_10x_sam(sam, pattern)

    process_sc_list(l)
}

#' @rdname bc_extract_sc_sam 
#' @export
bc_extract_sc_bam <- function(
    bam,
    pattern,
    cell_barcode_tag = "CR",
    umi_tag = "UR"
    ) {
    ## check if the bam file exists
    if(!file.exists(bam)) {
        stop("The input bam file does not exist.")
    }

    # Define a temporary directory
    tmp_dir <- tempdir()

    # Define the path for the SAM file
    sam_file_path <- file.path(tmp_dir, "output.sam")

    # Create a BamFile object 
    # bam_file <- Rsamtools::BamFile(bam)

    # Read the BAM file and write it to a SAM file
    cat("Start to convert bam file to sam file.\n")
    cat("sam file path: ", sam_file_path, "\n")
    Rsamtools::asSam(bam, sub("\\.sam$", "", sam_file_path), overwrite = TRUE)

    out = bc_extract_sc_sam(sam_file_path, pattern, cell_barcode_tag, umi_tag)
    file.remove(sam_file_path)

    out
}

#' Extract barcode from single-cell sequencing fastq file
#'
#' \code{bc_extract_10X_fastq} can extract cellular barcode, UMI, and lineage barcode
#' sequences from 10X Genomics scRNASeq fastq file. This function can process
#' the barcodes in the scRNASeq fastq file or target amplified fastq files directly.
#'
#' @param fq1 A string, the fastq file contains the cellular barcode and lineage
#' barcode
#' @param fq2 A string, it is optional, it provides the second fastq file
#' contains the cellular barcode and lineage barcode. Two fastq files will be
#' concatenated for the barcode extraction
#' @param patternCellBarcode A string, defines the regular expression to match
#' the single cell cellular barcode sequence. The expected sequence should be in
#' the first catch. Please see the documents of
#' \code{\link[CellBarcode]{bc_extract}} and example for more information.
#' @param patternUMI A string, defines the regular expression to match the UMI
#' sequence. The expected sequence should be in the first catch. Please see the
#' documents of \code{\link[CellBarcode]{bc_extract}} and example for more
#' information.
#' @param patternBarcode the regular expression to match the lineage barcode. The
#' expected sequence should be in the first catch. Please see the documents of
#' \code{\link[CellBarcode]{bc_extract}} and example for more information.
#' @details 
#' 
#' It should take some effort to define the regular expression to match the
#' barcode sequence. Here I also provide the example to extract the barcode from
#' 10X Genomics scRNASeq results. It also can be used to extract the barcode from
#' other system.
#'
#' The function can process the barcodes in the scRNASeq fastq file or target 
#' amplified fastq files. For the 10X scRNASeq fastq file, the cellular barcode is 
#' in the first 16bp of the read1, the UMI is in the next 12bp, and the lineage
#' barcode is in the read2.
#' 
#' The usage of the function will be like this:
#' 
#' \preformatted{
#' bc_extract_sc_fastq(
#'    fq1 = "read1.fastq.gz",
#'    fq2 = "read2.fastq.gz",
#'    patternCellBarcode = "(.{16})",
#'    patternUMI = ".{16}(.{12})",
#'    patternBarcode = "CGAAGTATCAAG(.+)CCGTAGCAAG"
#' )
#' }
#' 
#' @return 
#' A BarcodeObj object with each cell as a sample.
#' @seealso \code{\link[CellBarcode]{bc_extract}}, \code{\link[CellBarcode]{bc_extract_sc_sam}},
#'
#' @rdname bc_extract_sc_fastq
#' @export
bc_extract_sc_fastq <- function(
    fq1,
    fq2 = NULL,
    patternCellBarcode = NULL,
    patternUMI = NULL,
    patternBarcode = NULL
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
