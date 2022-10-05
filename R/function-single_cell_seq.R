#' Parse 10X Genomic scRNASeq sam file
#'
#' \code{bc_extract_10XscSeq} can extract cell barcode, UMI and lineage barcode
#' sequences from 10X Genomics scRNASeq bam file. This funciont can not process
#' bam file directly, user needs to uncompress the bam file to get sam file in
#' order to run this function.
#'
#' @param sam A string, define the un-mapped sequences 
#' @param pattern A string, define the regular expression to match the barcode
#' sequence. The barcode sequence should be in the first catch. Please see the
#' \code{\link[CellBarcode]{bc_extract}} for detail.
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
#' bc_extract_10X_scSeq(
#'   sam = sam_file,
#'   pattern = "AGATCAG(.*)TGTGGTA",
#'   cell_barcode_tag = "CR",
#'   umi_tag = "UR"
#' )
#'
#' @rdname bc_extract_10X_scSeq
#' @export
bc_extract_10X_scSeq <- function(
    sam, 
    pattern, 
    cell_barcode_tag = "CR", 
    umi_tag = "UR"
    ) {
    if(!file.exists(sam)) {
        stop("The input bam file is not exist.")
    }
    d = parse_10x_scSeq(sam, pattern)
    setDT(d)
    d = d[, .(count = sum(count)), by = .(cell_barcode, umi, barcode_seq)]

    as.data.frame(d)
}
