#' Parse 10X Genomic scRNASeq bam file
#' 
#' @export
bc_extract_10X_bam <- function(bam, pattern) {
    if(!file.exists(bam)) {
        stop("The input bam file is not exist.")
    }
    parse_10x_bam(bam, pattern)
}
