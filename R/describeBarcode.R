#' @export
bc_diversity.BarcodeObj = function(x, clean = NULL, plot = T) {

  # BUG: Could not work.

  # TODO: output UMI relate information when UMI is available
  #       - reads per UMI barcode
  #       - UMI number
  #       - UMI barcode Number
  #       - ...
  count = barcode_seq = NULL

  # x is BarcodeObj
  if (is.null(clean)) {
    if ("cleanBc" %in% names(x)) {
      clean = TRUE
    } else {
      clean = FALSE
    }
  }
  if (clean) {
    d = x$cleanBc
  } else {
    d = lapply(x$messyBc, function(x) { x[, .(count = sum(count)), by = barcode_seq] })
    names(d) = names(x$messyBc)
  }

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
    ), by = "sample_name"]
}

plot_lorenz_curve = function(x) {
  p_cum = sample_name = NULL

  # x is data.table with columns: count, barcode_seq, sample_name
  d = x[, .(p_cum = cumsum(sort(count / sum(count), decreasing = T)), x = 1:length(count)), by = sample_name]
  g = ggplot(d) + aes(x = x, y = p_cum, color = sample_name, fill = sample_name) + geom_line() + geom_abline(slope = 1 / max(d$x), alpha = 0.3, color = 'black') + theme_bw()
  g = g + labs(x = "Barcodes", y = "Cumulative reads per sample", color = "Sample")
  g
}

plot_stack_barplot = function(x) {
  # TODO: Too many barcodes, slow the drawing process
  sample_name = count = barcode_seq = NULL
  # x is data.table with columns: count, barcode_seq, sample_name
  
  barcodes_n_max = x$barcode_seq %>% unique %>% length #max(x[, .N, by = sample_name]$N)
  g = ggplot(x) + aes(x = sample_name, y = count, fill = barcode_seq) + geom_bar(stat = "identity", position = "fill")
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

# TODO: How to deal the network graph
#' Draw Network graph to show the similarity of barcode
#'
# @export
bc_plotNetwork = function(x, ...) UseMethod("bc_plotNetwork", x)

#' Draw Network graph to show the similarity of barcode
#'
#' @param barcodeObj A BarcodeObj
#' @param sample1 A string, sample name
#' @param sample2 A string or NULL, the second sample while comparing the barcode between two samples
#' @param topN A integer, only the barcodes with top N depth or UMI-barcode count are used
#' @param type A integer, "raw" or "clean" correspond to the uncorrected or corrected barcde data
#' @param UMI_depth_threshold_n A integer, only the barcode with read depth or UMI-barcode count greater than the value will be used
#' @return NULL
# @export
bc_plotNetwork.BarcodeObj = function(barcodeObj, sample1 = NULL, sample2 = NULL, topN = 100, type = "raw", UMI_depth_threshold_n = 5, ...) {

  count = barcode_seq = NULL  # due to NOTE in check

  if (type == "clean" & !is.null(barcodeObj$cleanBc)) {
    d = barcodeObj$cleanBc
  } else if (type == "raw" & !is.null(barcodeObj$messyBc)) {
    d = barcodeObj$messyBc
  } else {
    stop("No barcodes found.")
  }

  if (is.null(sample2)) {
    d1 = data.table(d[[sample1]])
    if (type == "raw") {
      d1 = d1[count > UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
    }
    d1 = d1[order(count, decreasing = T)][1:min(topN, nrow(d1))]
    clusterSeq(d1$barcode_seq, weight = d1$count, sample1 = sample1, ...)
  } else if (!is.null(sample1) & !is.null(sample2)) {
    d1 = data.table(d[[sample1]])
    d2 = data.table(d[[sample2]])
    if (type == "raw") {
      d1 = d1[count > UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
      d2 = d2[count > UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
    }
    d1 = d1[order(count, decreasing = T)][1:min(topN, nrow(d1))]
    d2 = d2[order(count, decreasing = T)][1:min(topN, nrow(d2))]
    clusterSeq(d1$barcode_seq, d2$barcode_seq, weight = c(d1$count, d2$count), sample1 = sample1, sample2 = sample2, ...)
  } else {
    stop("Error sample names.")
  }
}

clusterSeq = function(seqs, seqs2 = NA, weight = 1, isLabel = T, highLight = NA, highLightLink=F, sample1, sample2 = NA) {
  if (is.na(seqs2[1])) {
    drawNetwork(calcDistance(seqs), weight, isLabel, highLight)
  } else {
      drawNetwork(
        calcDistance(c(seqs, seqs2)),
        weight,
        isLabel,
        highLight,
        col = rep(
          c(sample1, sample2),
          c(length(seqs), length(seqs2))
          ),
        highLightLink = highLightLink
      )
  }
}

calcDistance = function(seqs, weight = "", label = NULL) {

  distance = utils::adist(seqs)
  if (length(label) == nrow(distance)) {
    rownames(distance) = make.names(label)
  } else {
    rownames(distance) = paste0(weight, seqs)
  }
  distance
}

drawNetwork = function(distance, weight, isLabel, highLight, col = NA, highLightLink = F, sample1, sample2 = NA) {

  x = y = xend = yend = label = size_Log10 = edge_alpha = edge_size = vertex_col = edge_col = NULL

  if (is.na(col[1])) {
    col = rep(sample1, length(nrow(distance)))
  }
  distance[distance > 1] = 3
  distance[distance == 0] = 2
  distance[distance == 3] = 0
  net = network::network(distance, directed = F)
  if (isLabel) {
    network::set.vertex.attribute(net, "label", rownames(distance))
  }
  network::set.vertex.attribute(net, "size_Log10", log10(weight))
  network::set.vertex.attribute(net, "vertex_col", col)
  network::`%e%`(net, "dist") = distance
  network::`%e%`(net, "edge_col") = c("Dist1", "Dist0")[network::`%e%`(net, "dist")]
  network::`%e%`(net, "edge_alpha") = c("Dist1" = 0.3, "Dist0" = 0.9)[network::`%e%`(net, "edge_col")]
  network::`%e%`(net, "edge_size") = c("Dist1" = 0, "Dist0" = 1)[network::`%e%`(net, "edge_col")]

  ggnet = data.table::data.table(ggnetwork::ggnetwork(net))

  set.seed(1)
  if (!is.na(highLight[1])) {
    ggplot2::ggplot(ggnet) + ggplot2::aes(x = x, y = y, xend = xend, yend = yend, label = label) +
      ggnetwork::geom_edges(ggplot2::aes(linetype = edge_col), alpha = 0.2, ncp = 0) +
      ggnetwork::geom_nodes(ggplot2::aes(size = size_Log10,  color = label %in% highLight), alpha = 0.7) +
      ggnetwork::theme_blank() + labs(shape = "Sample", color = "Sample", size = "Count", linetype = "Hammer Dist")
  } else {
    if (highLightLink) {
      ggplot2::ggplot(ggnet) + ggplot2::aes(x = x, y = y, xend = xend, yend = yend, label = label) +
        ggnetwork::geom_edges(aes(linetype = edge_col, col = edge_col, alpha = edge_alpha, size = edge_size), ncp = 0, lineend="square") +
        ggnetwork::geom_nodes(aes(size = size_Log10,  shape = vertex_col), alpha = 0.1) +
        ggnetwork::theme_blank() + labs(shape = "Sample", color = "Sample", size = "Count", linetype = "Hammer Dist")
    } else {
      ggplot(ggnet) + ggplot2::aes(x = x, y = y, xend = xend, yend = yend, label = label) +
        ggnetwork::geom_edges(ggplot2::aes(linetype = edge_col), alpha = 0.2, ncp = 0) +
        ggnetwork::geom_nodes(ggplot2::aes(size = size_Log10,  shape = vertex_col, color = vertex_col), alpha = 0.7) +
        ggnetwork::theme_blank() + labs(shape = "Sample", color = "Sample", size = "Count", linetype = "Hammer Dist")
    }
  }
}
