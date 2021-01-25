# library(plyr)
# library(stringr)
# library(readr)
# library(ggplot2)
# library(network)
# library(ShortRead)
# library(ggpubr)
# library(data.table)
# library(magrittr)

# Design
# library(Bc)
# 
# BarcodeObj -> List(rawBc, meta, rawBc, messyBc, cleanBc)
# 
# readFq() -> class BarcodeObj(meta, rawBc )
# 
# runQc() -> class BarcodeObj(meta, rawBc)
# 
# filterFq.Barcode() -> class BarcodeObj(meta, rawBc)
# 
# extractBc.Barcode() -> class BarcodeObj(meta, rawBc, messyBc)
# 
# cureBc.Barcode() -> class BarcodeObj(meta, rawBc, messyBc, cleanBc)
# 
# countBc.Barcode() -> calss BarcodeObj(meta, rawBc, messyBc, cleanBc, ?)
# 
# readFq -> runQc -> filterFq -> extractBc -> cureBc -> cureBc
# barcodeObj[, 1]
# barcodeObj[, "sample_name1"]
# barcodeObj[, 1:2]
# barcodeObj[, c("sample_name1", "sample_name2")]
# barcodeObj["AACT", ]
# subset(barcodeObj, name = c("sample_name1", "sample_name2"))
# subset(barcodeObj, barcode_white_list = c("AAAA", "TTT"), bacode_black = c("GGG", "CCC"))
# barcodeObj1 + barcodeObj2
# barcodeObj - black_list
# barcodeObj * white_list
# samplenames(bc_obj)
# bc2genBarCode(bc_obj)

# TODO:
# filterFq() -> output the filtering detail
# make function to cross talk with genBaRcode and other packages


#' Subset operation of BarcodeObj object
#' 
#' @param x A vector or string. Select the barcode subset from BarcodeObj
#' @param y A vector or string. Select the sample subset from BarcodeObj
#' @return A BarcodeObj
#' @export
#' @examples
#' # Not run 
#' bc_obj[c("AACCTT", "AACCTT"), ]
#' bc_obj[, "sample1_rep1"]
#' bc_obj[, c("sample1_rep1", "sample1_rep2")]
"[.BarcodeObj" = function(barcodeObj, x = NULL, y = NULL) {
  if (is.null(x) & is.null(y)) {
    error("Error message. TODO.")
  }
  if (is.null(x)) {
    return(subset(d, sample = y))
  }
  if (is.null(y)) {
    return(subset(d, barcode = x))
  }
  return(subset(d, sample = y, barcode = x))
  return(barcodeObj)
}

#' Subset operation of BarcodeObj object
#' 
#' @param sample A vector or string. Select the sample subset from BarcodeObj
#' @param barcode A vector or string. Select the barcode subset from BarcodeObj
#' @param black_list A vector or string. Remove the the barcode in the list.
#' @return A BarcodeObj
#' @export
#' @examples
#' # Not run 
#' subset(bc_obj, barcode = c("AACCTT", "AACCTT"))
#' subset(bc_obj, sample = "sample1_rep1", barcode = c("AACCTT", "AACCTT"))
#' subset(bc_obj, sample = c("sample1_rep1", "sample1_rep2"), barcode = c("AACCTT", "AACCTT"))
subset.BarcodeObj = function(barcodeObj, sample = NULL, barcode = NULL, black_list = NULL) {
  # TODO: The funciton only can apply the operation to the `messyBc` and `cleanBc`. We need to make it 
  # capable to apply the selection to all information in the object.
  if (is.null(sample) & is.null(barcode)) {
    error("Error message. TODO.")
  }
  if (!is.null(barcode)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc %<>% lapply(function(d) {
        d[reads_seq %in% barcode]
      })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc %<>% lapply(function(d) {
        d[reads_seq %in% barcode]
      })
    }
  }
  if (!is.null(black_list)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc %<>% lapply(function(d) {
        d[!(reads_seq %in% black_list)]
      })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc %<>% lapply(function(d) {
        d[!(reads_seq %in% black_list)]
      })
    }
  }
  if (!is.null(sample)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc %<>% magrittr::extract(sample)
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc %<>% magrittr::extract(sample)
    }
  }
  return(barcodeObj)
}

#' Combined the BarcodeObj
#'
#' @param x A BarcodeObj
#' @param y A BarcodeObj
#' @return A BarcodeObj
#' @export
"+.BarcodeObj" = function(x, y) {
  # TODO: Apply the merge to all parts of the data
  #       How to deal when two BarcodeObj have the same samples
  x$messyBc %<>% append(y$messyBc)
  x$cleanBc %<>% append(y$cleanBc)
  x
}

#' Remove the barcode in black list
#' 
#' @param x A BarcodeObj
#' @param y A string or vector, the barcode black list to revmoe
#' @return A BarcodeObj
"-.BarcodeObj" = function(x, y) {
  subset(x, black_list = y)
}

#' Get the sample names from BarcodeObj
#'
#' @export
samplenames = function(x, ...) UseMethod("samplenames", x)

#' Get the sample names from BarcodeObj
#'
#' @param barcodeObj BarcodeObj
#' @return A string vector
#' @export
samplenames.BarcodeObj = function(barcodeObj) {
  barcodeObj$sampleNames
}


#' Select the barcode by white list
#'
#' @param x A BarcodeObj
#' @param y A string or vector, the barcode white list
#' @return A BarcodeObj
"*.BarcodeObj" = function(x, y) {
  subset(x, barcode = y)
}


#' Read the fastq file
#' 
#' @param fastqFiles A vector or string. The fastq file location
#' @param sampleNames A vector or string, which should be the same length to the `fastqFiles`. The sample name corespond to the input fastq files
#' @return A BarcodeObj
#' @export
#' @examples
readFq = function(fastqFiles, sampleNames) {
  ## Prepare the barcodeObj
  barcodeObj = list()
  class(barcodeObj) = "BarcodeObj"
  # TODO: if no sampleNames, generate the name from fastqFiles.
  meta = data.frame(sample_name = sampleNames)
  rownames(meta) = sampleNames

  ## read in fastq files
  rawBc = lapply(fastqFiles, ShortRead::readFastq)
  names(rawBc) = sampleNames

  ## generate reads quality information
  # TODO: Create a new function to check the quality of fastq file.
  qaResult = ShortRead::qa(fastqFiles)
  meta$raw_read_count = qaResult[["readCounts"]][, "read"]
  ## Base per cycle
  raw_baseCall = qaResult[["perCycle"]][["baseCall"]]
  raw_baseCall$sample_name = sampleNames[match(raw_baseCall$lane, basename(fastqFiles))]
  raw_baseCall$lane = NULL
  ## Base quality per cycle
  raw_quality = qaResult[["perCycle"]][["quality"]]
  raw_quality$sample_name = sampleNames[match(raw_quality$lane, basename(fastqFiles))]
  raw_quality$lane = NULL
  ## high frequency sequence
  raw_frequentSequences = qaResult[["frequentSequences"]]
  raw_frequentSequences$sample_name = sampleNames[match(raw_frequentSequences$lane, basename(fastqFiles))]
  raw_frequentSequences$lane = NULL
  ## save the qc information
  qc = list()
  qc$raw_baseCall = raw_baseCall
  qc$raw_quality = raw_quality
  qc$raw_frequentSequences = raw_frequentSequences

  ## assemble data into barcodeObj
  barcodeObj$rawBc = rawBc
  barcodeObj$sampleNames = sampleNames
  barcodeObj$qc = qc
  barcodeObj$meta = meta
  barcodeObj
}


#' Visualize the Quality information of Fastq file
qcSummary = function(x, ...) UseMethod("qcSummary", x)

#' Visualize the Quality information of Fastq file
#' 
#' @param barcodeObj
#' @return NULL
#' @examples
qcSummary.BarcodeObj = function(barcodeObj) {
  # TODO:
  # Raw data:
  #   Most frequency reads
  #   Reads length distribution
  #   nucleic acide per base
  #   quality per base
  d = data.table(barcodeObj$qc$raw_baseCall)
  d[, base_num := sum(Count), by = .(Cycle, sample_name)]
  d[, base_percent := Count / base_num, by = .(Base, Cycle, sample_name)]
  p1 = ggplot(d, aes(Cycle, sample_name, fill = Base, alpha = base_percent)) + geom_tile(color = "black") + labs(y = "Sample Name") + theme_bw()

  d = data.table(barcodeObj$qc$raw_quality)
  d[, quality_median := median(Score), by = .(Cycle, sample_name)]
  p2 = ggplot(d, aes(Cycle, sample_name, fill = quality_median)) + geom_tile(color = "black") + labs(y = "Sample Name", fill = "Median Base Quality") + theme_bw()

  ggpubr::ggarrange(p1, p2, nrow = 2, ncol=1) %>% print
}

#' Filter the reads by the QC threshold 
filterFq = function(x, ...) UseMethod("filterFq", x)

#' Filter the reads by the QC threshold 
#'
#' @param barcodeObj A BarcodeObj
#' @param min_quality A number. The minimum reads base quality threshold
#' @param min_mean_quality A number. The minimum average reads base quality threshold
#' @param min_length A number. The minimum length of the reads
#' @export
#' @return A BarcodeObj
filterFq.BarcodeObj = function(barcodeObj, min_quality = 20, min_mean_quality = 30, min_length = 30) {
  # TODO: The information about the removed reads
  lapply(barcodeObj$rawBc, function(fq) {
    is_min_quality = apply(as(ShortRead::quality(fq), "matrix"), 1, min, na.rm = TRUE) >= min_quality
    is_mean_quality = apply(as(ShortRead::quality(fq), "matrix"), 1, mean, na.rm = TRUE) >= min_mean_quality
    is_min_length = width(ShortRead::sread(fq)) >= min_length
    fq[is_min_quality & is_mean_quality & is_min_length]
}) -> barcodeObj$rawBc
  barcodeObj
}


#' Extract barcode from reads
#' 
#' @param barcodeObj A BarcodeObj
#' @param pattern A string. The regular expression to match the barcode which capture pattern
#' @param maxLDist A integer. The mismatch threshold for barcode matching
#' @param pattern_type A vector. It defines the barcode (and UMI) capture pattern
#' @param costs A vector. Define the weight for each mismatch events
#' @return A barcodeObj
#' @export
extractBc = function(x, ...) UseMethod("extractBc", x)
extractBc.BarcodeObj = function(barcodeObj, pattern = "", maxLDist = 2, pattern_type = c(barcode = 1), costs = list(sub = 1, ins = 99, del = 99)) {
  mclapply(barcodeObj$rawBc, function(bc_reads) {
    reads_freq = ShortRead::tables(bc_reads, n=Inf)$top
    reads_seq = names(reads_freq)
    m = aregexec(pattern, reads_seq, fixed = F, max.distance = maxLDist, costs = list(sub = 1))
    seq_l = regmatches(reads_seq, m)

    seq_v = laply(seq_l, function(x) { x[1] })
    barcode_v = laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
    if ("UMI" %in% names(pattern_type)) {
      umi_v = laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
      d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
    } else {
      d = data.table(reads_seq = reads_seq, match_seq = seq_v, barcode_seq = barcode_v, count = reads_freq)
    }
    d %>% na.omit
}) -> messyBc
  names(messyBc) = names(bc_obj$rawBc)
  barcodeObj$messyBc = messyBc
  barcodeObj
}

#' Print out the summary of BarcodeObj
#' 
#' @param barcodeObj A BarcodeObj
#' @return The general information about the BarcodeObj
#' @export
print.BarcodeObj = function(barcodeObj) {
  cat("Bonjour le monde. This is a baby barcode Object.\n----------\n")
  cat("It contains:\n")
  ## The items in the list
  cat(names(barcodeObj) %>% paste(collapse = "  "), "\n----------\n")
  if (!is.null(barcodeObj$messyBc)) {
    ## How many samples in messyBc item
    cat("$messyBc:", length(barcodeObj$messyBc), "Samples for uncleaned barcodes\n")
    for (i in names(barcodeObj$messyBc)) {
      ## The number of barcode in each sample
      cat("    In sample", paste0("$", i), "there are:", nrow(barcodeObj$messyBc[[i]]), "barcodes\n")
    }
    cat("\n")
  }
  if (!is.null(barcodeObj$cleanBc)) {
    ## How many samples in cleanBc iterm
    cat("$cleanBc: ", length(barcodeObj$cleanBc), "Samples for uncleaned barcodes\n")
    for (i in names(barcodeObj$cleanBc)) {
      ## The number of barcode in each sample
      cat("    In sample", paste0("$", i), "there are:", nrow(barcodeObj$cleanBc[[i]]), "barcodes\n")
    }
    cat("\n")
  }
}

view = function(x, ...) UseMethod("view", x)
view.BarcodeObj = function(barcodeObj) {
  print.default(barcodeObj)
}


#' Stair plot used to find a optimized cutoff point
stairPlot = function(x, ...) UseMethod("stairPlot", x)

#' Stair plot used to find a optimized cutoff point
#'
#' @param barcodeObj A BarcodeObj with extracted barcode.
#' @return NULL
#' @export
stairPlot.BarcodeObj = function(barcodeObj, sample = 1) {
  # TODO: The stair plot only apply to messyBc iterm, do we need to apply it to cleanBc
  #       Return the optimized cutoff point.
  messyBc = barcodeObj$messyBc
  d = messyBc[[1]] %>% data.table
  d = d[, .(count = sum(count)), by = barcode_seq][order(count)]
  y = cumsum(d$count)
  x = d$count
  plot(x, log10(y), type = 's', xlab = "Reads count", ylab = "Accumlated Reads count")
}

#' Correct the PCR and sequencing mutation in the barcode sequence
cureBc = function(x, ...) UseMethod("cureBc", x)

#' Correct the PCR and sequencing mutation in the barcode sequence
#' 
#' @param barcodeObj A BarcodeObj
#' @param depth_threshold A integer, only a barcode with sequencing depth larger than that will be processed
#' @param hamming_dist_threshold A integer, the hamming distance threshold of two distinct barcode
#' @param with_umi A bool value, True when UMI is used
#' @return A BarcodeObj
#' @export
cureBc.BarcodeObj = function(
  barcodeObj
  , depth_threshold = 0
  , barcode_count_threshold = 1000
  , hammer_dist_threshold = 1
  , with_umi = F
  ) {
  # TODO: Use the barcode distribution and sequence base pair quality to correct

  messyBc = barcodeObj$messyBc

  if (with_umi) {
    ## When UMI used, count the UMI-barcode
    messyBc = lapply(messyBc,
      function(d) {
        d = data.table(d)
        d = d[count > depth_threshold, .(count = .N), by = .(barcode_seq)]
      }
    )
  } else {
    ## If no UMI, count the reads
    messyBc = lapply(messyBc,
      function(d) {
        d = data.table(d)
        d = d[, .(count = sum(count)), by = barcode_seq][count > depth_threshold]
      }
    )
  }
  ## Do the correction
  correct_out = lapply(messyBc,
    function(d) {
      seq_v = d$barcode_seq
      count_v = d$count
      seq_correct(seq_v, count_v, barcode_count_threshold, hammer_dist_threshold)
    }
  )

  ## Prepare output
  cleanBc = lapply(correct_out, 
    function(d) {
      data.frame(d$seq_freq[order(d$seq_freq$count, decreasing = T), ])
    }
  )

  ## The correction log, in another words, who correction has been done
  cleanProc = lapply(correct_out,
    function(d) {
      d$link_table
    }
  )

  names(cleanBc) = names(messyBc)
  names(cleanProc) = names(messyBc)

  ## save the result
  barcodeObj$cleanBc = cleanBc
  barcodeObj$cleanProc = cleanProc
  barcodeObj
}

#' Readout the barcode count
countBc = function(x, ...) UseMethod("countBc", x)

#' Readout the barcode count
#' 
#' @param barcodeObj A BarcodeObj
#' @return A list, the item is a data.frame, with column names of `barcode_seq` and `count`, the list name is sample names
#' @export
#' @examples
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
    calcDistance(seqs) %>% drawNetwork(weight, isLabel, highLight)
  } else {
    calcDistance(c(seqs, seqs2)) %>%
      drawNetwork(
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

  distance = adist(seqs)
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
  net %e% "dist" = distance
  net %e% "edge_col" = c("Dist1", "Dist0")[net %e% "dist"]
  net %e% "edge_alpha" = c("Dist1" = 0.3, "Dist0" = 0.9)[net %e% "edge_col"]
  net %e% "edge_size" = c("Dist1" = 0, "Dist0" = 1)[net %e% "edge_col"]

  ggnet = ggnetwork::ggnetwork(net) %>% data.table

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

## Foreign

#' @title Generate genBarCode object from Bc BarcodeObj
#'
#' @description Use the uncorrected barcode to generate the genBarCode data object
#' @param barcodeObj A BarcodeObj
#' @param BC_backbone A string, the backbone of the barcode
#' @param resDir A string, the location to store the data of genBarCode object
#' @return A list, contains genBarCode object
bc2genBarCode = function(barcodeObj, BC_backbone = NULL, resDir = getwd()) {

  messyBc = BarcodeObj$messyBc
  lapply(1:length(messyBc), function(i) {
    d = messyBc[[i]][, .(read_count = sum(count)), by = barcode_seq][order(count, decreasing = T)][, .(read_count = count, barcode = factor(barcode_seq))] %>% 
      as.data.frame
    return(methods::new(Class = "BCdat", reads = d, results_dir = resDir,
      label = names(messyBc)[i], BC_backbone = BC_backbone))
   }) 
}
