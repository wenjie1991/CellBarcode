# library(ShortRead)
# library(data.table)
# library(ggplot2)
# library(plyr)
# library(magrittr)
# library(readr)
# library(stringr)
# library(ggpubr)
# library(network)

# library(Bc)
# 
# BarcodeObj -> List(rawBc, meta, rawBc, messyBc, cleanBc)
# 
# readFq() -> class BarcodeObj -> rawBc
# 
# runQc() -> class BarcodeObj -> meta
# 
# filterFq.Barcode() -> class BarcodeObj -> rawBc
# 
# extractBc.Barcode() -> class BarcodeObj -> messyBc
#
# cureBc.Barcode() -> class BarcodeObj -> cleanBc
# 
# countBc.Barcode() -> calss BarcodeObj -> meta

readFq = function(fastqFiles, sampleNames) {
    barcodeObj = list()
    qc = list()
    meta = data.frame(sample_name = sampleNames)
    rownames(meta) = sample_name_vector

    rawBc = lapply(fastq_file_vector, readFastq)
    names(rawBc) = sample_name_vector

    qaResult = qa(fastq_file_vector)
    meta$raw_read_count = qaResult[["readCounts"]][, "read"]

    d = qaResult[["perCycle"]]
    raw_baseCall = qaResult[["perCycle"]][["baseCall"]]
    raw_baseCall$sample_name = sample_name_vector[match(raw_baseCall$lane, basename(fastq_file_vector))]
    raw_baseCall$lane = NULL

    raw_quality = qaResult[["perCycle"]][["quality"]]
    raw_quality$sample_name = sample_name_vector[match(raw_quality$lane, basename(fastq_file_vector))]
    raw_quality$lane = NULL

    raw_frequentSequences = qaResult[["frequentSequences"]] 
    raw_frequentSequences$sample_name = sample_name_vector[match(raw_frequentSequences$lane, basename(fastq_file_vector))]
    raw_frequentSequences$lane = NULL

    qc$raw_baseCall = raw_baseCall
    qc$raw_quality = raw_quality
    qc$raw_frequentSequences = raw_frequentSequences

    barcodeObj$rawBc = rawBc
    barcodeObj$sampleNames = sample_name_vector
    barcodeObj$qc = qc
    barcodeObj$meta = meta
    class(barcodeObj) = "BarcodeObj"
    barcodeObj
}

# TODO: 
# Raw data:
#   Most frequency reads
#   Reads length distribution
#   nucleic acide per base
#   quality per base
qcSummary = function(x, ...) UseMethod("qcSummary", x)
qcSummary.BarcodeObj = function(barcodeObj) {
    d = data.table(barcodeObj$qc$raw_baseCall)
    d[, base_num := sum(Count), by = .(Cycle, sample_name)]
    d[, base_percent := Count / base_num, by = .(Base, Cycle, sample_name)]
    p1 = ggplot(d, aes(Cycle, sample_name, fill = Base, alpha = base_percent)) + geom_tile(color = "black") + labs(y = "Sample Name") + theme_bw()

    d = data.table(barcodeObj$qc$raw_quality)
    d[, quality_median := median(Score), by = .(Cycle, sample_name)]
    p2 = ggplot(d, aes(Cycle, sample_name, fill = quality_median)) + geom_tile(color = "black") + labs(y = "Sample Name", fill = "Median Base Quality") + theme_bw()

    ggarrange(p1, p2, nrow = 2, ncol=1) %>% print
}


filterFq = function(x, ...) UseMethod("filterFq", x)
filterFq.BarcodeObj = function(barcodeObj, min_quality = 20, min_mean_quality = 30, min_length = 30) {
    lapply(barcodeObj$rawBc, function(fq) {
        is_min_quality = apply(as(ShortRead::quality(fq), "matrix"), 1, min, na.rm = TRUE) >= min_quality
        is_mean_quality = apply(as(ShortRead::quality(fq), "matrix"), 1, mean, na.rm = TRUE) >= min_mean_quality
        is_min_length = width(sread(fq)) >= min_length
        fq[is_min_quality & is_mean_quality & is_min_length]
    }) -> barcodeObj$filteredBc
    barcodeObj
}

extractBc = function(x, ...) UseMethod("extractBc", x)
extractBc.BarcodeObj = function(barcodeObj, pattern = "", maxLDist = 2) {
    mclapply(barcodeObj$rawBc, function(bc_reads) {
        reads_freq = tables(bc_reads, n=Inf)$top
        reads_seq = names(reads_freq)
        # pattern = "([ACGC]{16})CCTCGAG(.+)GTAGCTACTACCGTAGCAAGCTCGAGAGTAGAC"
        pattern_type = c(UMI = 1, barcode = 2)
        m = aregexec(pattern, reads_seq, fixed = F, costs = list(sub = maxLDist))
        seq_l = regmatches(reads_seq, m)

        seq_v = laply(seq_l, function(x) { x[1] })
        umi_v = laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
        barcode_v = laply(seq_l, function(x) { x[pattern_type["barcode"] + 1] })
        d = data.table(reads_seq = reads_seq, match_seq = seq_v, umi_seq = umi_v, barcode_seq = barcode_v, count = reads_freq)
        d %>% na.omit
    }) -> messyBc
    names(messyBc) = names(bc_obj$rawBc)
    barcodeObj$messyBc = messyBc
    barcodeObj
}

print.BarcodeObj = function(barcodeObj) {
    print("Bonjour le monde. This is a baby barcode Object.")
}

view = function(x, ...) UseMethod("view", x)
view.BarcodeObj = function(barcodeObj) {
    print.default(barcodeObj)
}

# Rcpp::sourceCpp("./lib.cpp")
## test ##
# seq_v = c("AAAAAA", "AABAAA", "BAAAAA", "AAAAAB", "AAAAAA", "AAAAAA", "BBBB")
# count_v = c(3, 10, 11, 2, 4, 8, 1)
# data.table(seq_v, count_v)[order(count_v)]
# seq_correct(seq_v, count_v, 11, 2)
## test end ##

cureBc = function(x, ...) UseMethod("cureBc", x)
cureBc.BarcodeObj = function(
    barcodeObj
    , UMI_depth_threshold_n = 0
    , expect_barcode_n = Inf
    , barcode_distribution = NULL
    , barcode_count_threshold = 1000
    , hammer_dist_threshold = 1
    , with_umi = T
    ) {

    messyBc = barcodeObj$messyBc 

    if (with_umi) {
        ## TODO:
        ## clean the UMI and generate the barcode count 

        messyBc = lapply(messyBc,
            function(d) {
                d = data.table(d)
                d = d[count >= UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
            }
        )
    }

    cleanBc = lapply(messyBc,
        function(d) {
            seq_v = d$barcode_seq
            count_v = d$count
            out = seq_correct(seq_v, count_v, barcode_count_threshold, hammer_dist_threshold)
            out[order(out$count, decreasing = T), ]
        }
    )

    barcodeObj$cleanBc = cleanBc

    barcodeObj
}

countBc = function(x, ...) UseMethod("countBc", x)
countBc.BarcodeObj = function(barcodeObj) {
    barcodeObj$cleanBc
}

plotNetwork = function(x, ...) UseMethod("plotNetwork", x)
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
            d1 = d1[count >= UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
        }
        d1 = d1[order(count, decreasing = T)][1:min(topN, nrow(d1))]
        clusterSeq(d1$barcode_seq, weight = d1$count, sample1 = sample1, ...)
    } else if (!is.null(sample1) & !is.null(sample2)) {
        d1 = data.table(d[[sample1]])
        d2 = data.table(d[[sample2]])
        if (type == "raw") {
            d1 = d1[count >= UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
            d2 = d2[count >= UMI_depth_threshold_n, .(count = .N), by = .(barcode_seq)]
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


calcDistance = function(seqs, weight = "", label = NULL) {
    
    distance = adist(seqs) 
    if (length(label) == nrow(distance)) {
        rownames(distance) = make.names(label)
    } else {
        rownames(distance) = paste0(weight, seqs)
    }
    distance
}

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


