#' Readout the barcode count
countBc = function(x, ...) UseMethod("countBc", x)

#' Readout the barcode count
#'
#' @param barcodeObj A BarcodeObj
#' @return A list, the item is a data.frame, with column names of `barcode_seq` and `count`, the list name is sample names
#' @export
countBc.BarcodeObj = function(barcodeObj) {
  barcodeObj$cleanBc
}

#' Draw Network graph to show the similarity of barcode
plotNetwork = function(x, ...) UseMethod("plotNetwork", x)

#' Draw Network graph to show the similarity of barcode
#'
#' @param barcodeObj A BarcodeObj
#' @param sample1 A string, sample name
#' @param sample2 A string or NULL, the second sample while comparing the barcode between two samples
#' @param topN A integer, only the barcodes with top N depth or UMI-barcode count are used
#' @param type A integer, "raw" or "clean" correspond to the uncorrected or corrected barcde data
#' @param UMI_depth_threshold_n A integer, only the barcode with read depth or UMI-barcode count greater than the value will be used
#' @return NULL
#' @export
plotNetwork.BarcodeObj = function(barcodeObj, sample1 = NULL, sample2 = NULL, topN = 100, type = "raw", UMI_depth_threshold_n = 5, ...) {

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

## no export
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

## no export
calcDistance = function(seqs, weight = "", label = NULL) {

  distance = utils::adist(seqs)
  if (length(label) == nrow(distance)) {
    rownames(distance) = make.names(label)
  } else {
    rownames(distance) = paste0(weight, seqs)
  }
  distance
}

## no export
drawNetwork = function(distance, weight, isLabel, highLight, col = NA, highLightLink = F, sample1, sample2 = NA) {
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
