#' Calculate barcode diversity
#' 
#' @param barcodeObj A BarcodeObj
#' @param plot A bool, if draw the lorenz curve and barcode distribution curve
#' @return A data.frame with following columns
#'   total_reads: total resds number
#'   shannon_index: Shannon information
#'   equitability_index: Equitability index
#'   bit_index: Shannon bit information
#'   complexity: e^shannon_index
#' @examples
#' data(bc_obj)
#'
#' bc_obj = bc_cure(bc_obj)
#' bc_diversity(bc_obj)
#' @export
bc_diversity = function(barcodeObj, plot = T) {

  # BUG: Could not work.

  # TODO: output UMI relate information when UMI is available
  #       - reads per UMI barcode
  #       - UMI number
  #       - UMI barcode Number
  #       - ...
  count = barcode_seq = NULL

  # x is BarcodeObj
  if (!("cleanBc" %in% names(barcodeObj))) {
    stop("Please run bc_cure() first.")
    # d = lapply(barcodeObj$messyBc, function(barcodeObj) { barcodeObj[, .(count = sum(count)), by = barcode_seq] })
    # names(d) = names(barcodeObj$messyBc)
  }

  d = barcodeObj$cleanBc
  d = rbindlist(d, idcol = "sample_name")

  if (plot) {
    g1 = plot_lorenz_curve(d)
    g2 = plot_reads_per_barcode_distribution(d)
    # TODO: use common legend?
    egg::ggarrange(plots = list(g1, g2), nrow = 1, ncol = 2) %>% print
  }

  d[, .(
    total_reads = sum(count)
    , uniq_barcode = length(count)
    , shannon_index = calc_shannon_index(count)
    , equitability_index = calc_equitability_index(count)
    , bit_index = calc_bit_info(count)
    , complexity = calc_correct_barcode_number(count)
    ), by = "sample_name"] %>% as.data.frame
}

plot_lorenz_curve = function(x) {
  p_cum = sample_name = NULL

  # x is data.table with columns: count, barcode_seq, sample_name
  d = x[, .(p_cum = cumsum(sort(count / sum(count), decreasing = T)), x = 1:length(count)), by = sample_name]
  g = ggplot(d) + aes(x = x, y = p_cum, color = sample_name, fill = sample_name) + geom_line() + geom_abline(slope = 1 / max(d$x), alpha = 0.3, color = 'black') + theme_bw()
  g = g + labs(x = "Barcodes", y = "Cumulative reads per sample", color = "Sample")
  g
}

plot_reads_per_barcode_distribution = function(x) {

  # x is data.table with columns: count, barcode_seq, sample_name
  g = ggplot(x) + aes(x = count, color = sample_name) + geom_density() + theme_bw()
  g = g + labs(x = "Reads Count", y = "Density", color = "Sample", fill = "Sample")
  g
}

calc_gini_index = function(x) {
  # x is data.table with columns: count, barcode_seq
  # TODO: 
}

calc_shannon_index = function(x) {
  # x is numeric vector that contains the frequency of each barcode
  count_sum = sum(x)
  p = x / count_sum
  sum(-(p * log(p)))
}

calc_correct_barcode_number = function(x) {
  exp(calc_shannon_index(x))
}

calc_equitability_index = function(x) {
  calc_shannon_index(x) / log(length(x))
}

calc_bit_info = function(x) {
  # x is data.table with columns: count, barcode_seq
  count_sum = sum(x)
  p = x / count_sum
  sum(-(p * log2(p)))
}
