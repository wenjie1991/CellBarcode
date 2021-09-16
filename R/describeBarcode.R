#' Evaluate barcode diversity
#' 
#' \code{bc_diversity} evaluates sequence diversity metrics using the barcodes
#' data in the \code{cleanBc} element of \code{BarcodeObj} object. It also generates Lorenz curve and
#' barcode frequency distribution graphs.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param plot A logical value, if TRUE, draw the Lorenz curve and barcode
#' distribution graphs. 
#' @param log_x A logical value, if TRUE, the \code{x} axis is logarized.
#' @return A data.frame with following columns:
#' \itemize{
#'   \item \code{total_reads}: total read number.
#'   \item \code{uniq_barcode}: how many barcodes in the dataset.
#'   \item \code{shannon_index}: Shannon's diversity index or Shannon–Wiener index.
#'   \item \code{equitability_index}: Shannon's equitability.
#'   \item \code{bit_index}: Shannon bit information.
#' }
#'
#' @details
#' Followings are the metrics used for evaluating the barcode diversity:
#'
#' \emph{Richness}: The unique barcodes number \eqn{R}, it evaluates the richness of
#' the barcodes.
#'
#' \emph{Shannon index}: Shannon diversity index is weighted geometric
#' average of the proportion \eqn{p} of barcodes.
#' \deqn{ H' = - \sum_{i=1}^{R}p_ilnp_i }
#'
#' \emph{Equitability index}: Shannon equitability \eqn{E_H} characterize the
#' evenness of the barcodes, it is a value between 0 and 1, with 1 being
#' complete evenness.
#' \deqn{ E_H = H' / H'_{max} = H / ln(R) }
#'
#' \emph{Bit}:
#' Shannon entropy \eqn{H}, with a units of bit, 
#' \deqn{ H = - \sum_{i=1}^{R}p_ilog_2p_i }
#'
#' @examples
#' data(bc_obj)
#' 
#' # filter barcode by depth
#' bc_obj <- bc_cure_depth(bc_obj)
#'
#' # evaluate barcode diversity
#' bc_diversity(bc_obj)
#' @export
bc_diversity <- function(barcodeObj, plot = TRUE, log_x=TRUE) {

    # BUG: Could not work.

    # TODO: output UMI relate information when UMI is available
    #       - reads per UMI barcode
    #       - UMI number
    #       - UMI barcode Number
    #       - ...
    count <- barcode_seq <- NULL

    # x is BarcodeObj
    if (!("cleanBc" %in% names(barcodeObj))) {
        stop("Please run bc_cure_*() first.")
    }

    d <- barcodeObj$cleanBc
    d <- rbindlist(d, idcol = "sample_name")

    if (plot) {
        g1 <- plot_lorenz_curve(d)
        g2 <- plot_reads_per_barcode_distribution(d)

        if (log_x) {
            g1 <- g1 + scale_x_log10()
            g2 <- g2 + scale_x_log10()
        }
        # TODO: use common legend?
        egg::ggarrange(plots = list(g1, g2), nrow = 1, ncol = 2)
    }

    d[, .(
        total_barcode_reads = sum(count)
        , uniq_barcode = length(count)
        , shannon_index = calc_shannon_index(count)
        , equitability_index = calc_equitability_index(count)
        , bit_index = calc_bit_info(count)
        ), by = "sample_name"] %>% as.data.frame
}

plot_lorenz_curve <- function(x) {
    p_cum <- sample_name <- NULL

    # x is data.table with columns: count, barcode_seq, sample_name
    d <- x[, .(
        p_cum = cumsum(sort(count / sum(count), decreasing = TRUE)),
        x = seq_along(count)), by = sample_name]

    g <- ggplot(d) + 
        aes(x = x, y = p_cum, color = sample_name, fill = sample_name) + 
        geom_line() + 
        geom_abline(slope = 1 / max(d$x), alpha = 0.3, color = 'black') + theme_bw()

    g <- g + labs(x = "Barcodes", y = "Cumulative reads per sample", color = "Sample")
    g
}

plot_reads_per_barcode_distribution <- function(x) {

    # x is data.table with columns: count, barcode_seq, sample_name
    g <- ggplot(x) + aes(x = count, color = sample_name) + geom_density() + theme_bw()
    g <- g + labs(x = "Reads Count", y = "Density", color = "Sample", fill = "Sample")
    g
}

calc_gini_index <- function(x) {
    # x is data.table with columns: count, barcode_seq
    # TODO: 
}

calc_shannon_index <- function(x) {
    # x is numeric vector that contains the frequency of each barcode
    count_sum <- sum(x)
    p <- x / count_sum
    sum(-(p * log(p)))
}

calc_correct_barcode_number <- function(x) {
    exp(calc_shannon_index(x))
}

calc_equitability_index <- function(x) {
    calc_shannon_index(x) / log(length(x))
}

calc_bit_info <- function(x) {
    # x is data.table with columns: count, barcode_seq
    count_sum <- sum(x)
    p <- x / count_sum
    sum(-(p * log2(p)))
}
