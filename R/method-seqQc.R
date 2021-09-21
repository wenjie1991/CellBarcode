## Internal functions
###########################
get_base_quality_per_cycle <- function(bstringset) {

    # Output: data.table Cycle, Mean, Median, P5, P95
    quality_m <- bstringset
    quality_encode <- ShortRead::encoding(quality_m)
    m <- ShortRead::alphabetByCycle(quality_m)[-1, ]
    quality_v <- quality_encode[rownames(m)]

    y <- lapply(seq_len(ncol(m)), function(i) {
        S4Vectors::Rle(quality_v, m[, i])
})

    quality_stat <- vapply(y, function(x) {
        quantile_result <- stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
        names(quantile_result) <- c("P5", "P25", "Median", "P75", "P95")
        c(Mean = mean(x), quantile_result)
}, c(Mean = 0, P5 = 0, P25 = 0, Median = 0, P75 = 0, P95 = 0)) 
    quality_stat <- t(quality_stat)
    as.data.frame(data.table(Cycle = seq(1, nrow(quality_stat)), quality_stat))
}

get_base_freq_per_cycle <- function(dnastringset) {

    # Output: data.table Cycle, Count, Base
    base_m <- ShortRead::alphabetByCycle(dnastringset)
    base_m <- t(base_m[c("A", "C", "G", "T"), ])
    base_m <- data.table(Cycle = seq(seq_len(nrow(base_m))), base_m)
    as.data.frame(data.table::melt(
            base_m, 
            id.vars = "Cycle", 
            measure.vars = c("A", "C", "G", "T"), 
            value.name = "Count", 
            variable.name = "Base"))
}


bc_seqqc <- function(x) {

    output <- ShortRead::tables(x)

    if (is(x, "ShortReadQ")) {
        output$base_quality_per_cycle <- get_base_quality_per_cycle(x@quality)
        output$base_freq_per_cycle <- get_base_freq_per_cycle(x@sread)
    } else {
        output$base_freq_per_cycle <- get_base_freq_per_cycle(x)
    }

    output$summary <- c(
        total_read = length(x),
        p5_read_length = quantile(width(x), 0.05)[[1]],
        median_read_length = quantile(width(x), 0.50)[[1]],
        p95_read_length = quantile(width(x), 0.95)[[1]]
        )

    output <- BarcodeQc(
        top = output$top,
        distribution = output$distribution,
        base_quality_per_cycle = output$base_quality_per_cycle,
        base_freq_per_cycle = output$base_freq_per_cycle,
        summary = output$summary
        )

    output
}

plot_reads_depth_distribution <- function(distribution) {

    nReads <- nOccurrences <- NULL

    # distribution is data.frame with two columns "nOccurrences" and "nReads".
    d <- distribution
    g <- ggplot(d) + aes(y = nReads, x = nOccurrences) + geom_point() + 
        theme_bw() + scale_y_log10() + scale_x_log10()
    g
}

plot_base_percentage_distribution <- function(base_freq_per_cycle) {

    Cycle <- Count <- Base <- NULL

    # base_freq_per_cycle is data.frame with three columns "Cycle", "Count" and
    # "Base"
    d <- data.table(base_freq_per_cycle)
    d[, all_base_per_cycle := sum(Count), by = .(Cycle)]
    d[, base_percent := Count / all_base_per_cycle, by = .(Cycle, Base)]
    g <- ggplot(d) + aes(x = Cycle, y = base_percent, color = Base) + 
        geom_point() + geom_line() + theme_bw()
    g
}

plot_base_quality_distribution <- function(base_quality_per_cycle) {

    Cycle <- Median <- P95 <- P75 <- P25 <- P5 <- NULL

    # base_quality_per_cycle is data.table with columns Cycle, Mean, Median, P5,
    # P95
    d <- base_quality_per_cycle
    g <- ggplot(d) + 
        aes(
            x = Cycle,
            middle  = Median,
            upper = P75,
            lower = P25,
            ymax = P95,
            ymin = P5,
            group = Cycle) + 
        geom_boxplot(stat = "identity") + theme_bw()
    g
}


#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("ShortReadQ"), function(x) {
    # output: top, distribution (nOccurrences, nReads), base_quality_per_cycle,
    # base_freq_per_cycle
    bc_seqqc(x)
})


#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("DNAStringSet"), function(x) {
    # output: top, distribution (nOccurrences, nReads), base_freq_per_cycle
    bc_seqqc(x)
})

#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("data.frame"), function(x) {
    ## Input: data.frame(seq, freq)
    ## convert data.frame to DNAStringSet
    sequences <- x$seq
    freq <- x$freq
    x <- Biostrings::DNAStringSet(rep(sequences, freq))

    bc_seqqc(x)
})


#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("integer"), function(x) {
    x <- Biostrings::DNAStringSet(rep(names(x), x))

    bc_seqqc(x)
})


#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("character"), function(x, sample_name = basename(x)) {
    # TODO: use qa to do fast check for the quality
    ## generate reads quality information
    # TODO: Create a new function to check the quality of fastq file.
    if (length(x) == 1) {
        bc_seq_qc(ShortRead::readFastq(x))
    } else {
        qc_list <- lapply(x, function(f) {
            bc_seq_qc(ShortRead::readFastq(f))
        })
        names(qc_list) <- sample_name
        output <- BarcodeQcSet(qc_list=qc_list)
        output
    }
})


#' @rdname bc_seq_qc
#' @exportMethod bc_seq_qc
setMethod("bc_seq_qc", c("list"), function(x, sample_name = names(x)) {
    qc_list <- lapply(x, bc_seq_qc)
    if (is.null(sample_name)) {
        sample_name <- seq_along(x)
    }
    names(qc_list) <- sample_name
    output <- BarcodeQcSet(qc_list=qc_list)
    output
})

#' @rdname bc_seq_qc
#' @exportMethod bc_plot_seqQc
setMethod("bc_plot_seqQc", c("BarcodeQc"), function(x) {
    # BarcodeQc has: top, distribution (nOccurrences, nReads),
    # base_quality_per_cycle, base_freq_per_cycle
    # TODO: user can select which figure to draw
    g_list <- list()
    if ("distribution" %in% slotNames(x)) {
        g_list <- append(
            g_list,
            list(plot_reads_depth_distribution(x@distribution)))
    }
    if ("base_freq_per_cycle" %in% slotNames(x)) {
        g_list <- base::append(
            g_list, list(plot_base_percentage_distribution(x@base_freq_per_cycle)))
    }
    if ("base_quality_per_cycle" %in% slotNames(x)) {
        g_list <- base::append(g_list,
            list(plot_base_quality_distribution(x@base_quality_per_cycle)))
    }
    g_num <- length(g_list)
    col_n <- ceiling(g_num / 2)
    egg::ggarrange(plots = g_list, nrow = 2, ncol=col_n)
})

#' @rdname bc_seq_qc
#' @exportMethod bc_plot_seqQc
setMethod("bc_plot_seqQc", c("BarcodeQcSet"), function(x) {

    Count <- Cycle <- fileName <- base_num <- Base <- base_percent <- Median <- 
        NULL

    # TODO:
    # Raw data:
    #   Most frequency reads
    #   Reads length distribution
    #   nucleic acide per base
    #   quality per base

    if (length(x@qc_list) == 1) {
        bc_plot_seqQc(x@qc_list[[1]])
    } else {
        g_list <- list()

        d <- lapply(
            seq_along(x),
            function(i) {
                data.table(x@qc_list[[i]]@base_freq_per_cycle, fileName = names(x@qc_list)[i]) 
            }) %>% rbindlist()

        d[, base_num := sum(Count), by = .(Cycle, fileName)]
        d[, base_percent := Count / base_num, by = .(Base, Cycle, fileName)]

        p1 <- ggplot(d, aes(Cycle, fileName, fill = Base, alpha = base_percent)) + 
            geom_tile(color = "black") + labs(y = "Sample Name") + theme_bw()
        g_list <- append(g_list, list(p1))

        if (!is.null(names(x@qc_list[[1]]@base_quality_per_cycle))) {
            d <- lapply(seq_along(x), function(i) { 
                data.table(x@qc_list[[i]]@base_quality_per_cycle, sample_name = names(x@qc_list)[i]) 
            }) %>% rbindlist()

            p2 <- ggplot(d, aes(Cycle, sample_name, fill = Median)) + 
                geom_tile(color = "black") + 
                labs(y = "Sample Name", fill = "Median Base Quality") + theme_bw()

            g_list <- append(g_list, list(p2))
        }

        row_n <- length(g_list)
        egg::ggarrange(plots = g_list, nrow = row_n, ncol=1)
    }
})


#' @rdname bc_summary_seqQc
#' @exportMethod bc_summary_seqQc
setMethod("bc_summary_seqQc", c("BarcodeQcSet"), function(x) {
    res <- lapply(x@qc_list, function(x_i) {
        as.list(x_i@summary)
    }) %>% rbindlist(idcol = TRUE)
    names(res)[1] <- "sample_name"
    res
})


#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("ShortReadQ"), function(
    x, 
    min_average_quality = 30, 
    min_read_length = 0,
    N_threshold = 0) {

    # good base quality filter
    goodq <- ShortRead::srFilter(function(x) { 
        apply(as(Biostrings::quality(x), "matrix"), 1, mean, na.rm=TRUE) >= 
            min_average_quality 
    }, name="GoodQualityBases")

    # good sequences length filter
    goodlength <- ShortRead::srFilter(function(x) { width(ShortRead::sread(x)) >= 
        min_read_length}, name="GoodReadLength")

    # 'N' nucleic filter
    goodN <- ShortRead::nFilter(threshold = N_threshold)

    goodFinal <- ShortRead::compose(goodlength, goodq, goodN)
    x[goodFinal(x)]
})


#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("DNAStringSet"), function(
    x,
    min_read_length = 0,
    N_threshold = 0) {

    goodlength <- ShortRead::srFilter(function(x) {
        width(x) >= min_read_length
    }, name="GoodReadLength")

    goodN <- ShortRead::nFilter(threshold = N_threshold)
    goodFinal <- ShortRead::compose(goodlength, goodN)
    x[goodFinal(x)]
})


#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("data.frame"), function(x,
    min_read_length = 0,
    N_threshold = 0) {

    sequences <- x$seq
    freq <- x$freq
    x <- Biostrings::DNAStringSet(rep(sequences, freq))

    bc_seq_filter(x, min_read_length = min_read_length, N_threshold = N_threshold)
})

#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("character"), function(x,
    min_average_quality = 30,
    min_read_length = 0,
    N_threshold = 0,
    sample_name = basename(x)) {

    # TODO: when sample_name length is not right ... 
    fq_list <- lapply(x, ShortRead::readFastq)
    names(fq_list) <- sample_name

    bc_seq_filter(
        fq_list,
        min_average_quality = min_average_quality,
        min_read_length = min_read_length,
        N_threshold = N_threshold)
})

#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("integer"), function(x, min_read_length = 0, N_threshold = 0) {

    # TODO: Use Rle?
    x <- Biostrings::DNAStringSet(rep(names(x), x))
    bc_seq_filter(x, min_read_length = min_read_length, N_threshold = N_threshold)
})

#' @rdname bc_seq_filter
#' @exportMethod bc_seq_filter
setMethod("bc_seq_filter", c("list"), function(
    x, 
    min_average_quality = 30,
    min_read_length = 0,
    N_threshold = 0,
    sample_name = names(x)) {


    if (is(x[[1]], "ShortReadQ") | is(x[[1]], "character")) {
        output <- lapply(
            x,
            bc_seq_filter,
            min_average_quality = min_average_quality,
            min_read_length = min_read_length,
            N_threshold = N_threshold)

        names(output) <- sample_name
        output
    } else {
        output <- lapply(
            x,
            bc_seq_filter,
            min_read_length = min_read_length,
            N_threshold = N_threshold)

        names(output) <- sample_name
        output
    }
})

#' Show BarcodeQc object
#'
#' Show the summary of BarcodeQc object for pretty print.
#'
#' @param object A BarcodeQc object
#' @return Formated summary text.
#'
#' @rdname show
#' @exportMethod show
setMethod("show", c("BarcodeQc"), function(object) {
    cat("Sequnece QC, summary:\n",
        paste(
            "    ",
            paste(
                names(object@summary), object@summary, sep=": "
                ), 
            collapse="\n"), 
        "\n", sep="")
})

#' Show BarcodeQcSet object
#'
#' Show the summary of BarcodeQcSet object for pretty print.
#'
#' @param object A BarcodeQcSet object
#' @return Formated summary text.
#'
#' @rdname show
#' @exportMethod show
setMethod("show", c("BarcodeQcSet"), function(object) {
    cat("The sequence QC set, use `[]` to select sample:\n",
        paste(
            "    ",
            names(object@qc_list), collapse="\n"
            ), 
        "\n", sep="")
})

#' Subset the BarcodeQcSet
#'
#' Subset the BarcodeQcSet
#'
#' @param x A BarcodeQcSet object
#' @param i A integer vector or a character vector, specifying the selected
#' samples.
#' @param drop a logical value, if TRUE, when only one sample is selected, the
#' output will be a BarcodeQc object.
#' @return A BarcodeQcSet or BarcodeQc
#' @examples
#'
#' example_data <- system.file("extdata", "mef_test_data", package = "CellBarcode")
#' fq_files <- dir(example_data, "fastq.gz", full=TRUE)
#' qc_noFilter <- bc_seq_qc(fq_files) 
#' qc_noFilter[1:3]
#'
#' @exportMethod [
setMethod("[", c("BarcodeQcSet"), function(x, i, drop=TRUE) {
    if (length(i) == 1) {
        x@qc_list[[i]]
    } else {
        BarcodeQcSet(qc_list=x@qc_list[i])
    }
})

#' @rdname bc_names
#' @exportMethod bc_names
setMethod("bc_names", c("BarcodeQcSet"), function(x) {
    names(x@qc_list)
})

#' @rdname bc_names
#' @exportMethod bc_names
setMethod("bc_names<-", c("BarcodeQcSet"), function(x, value) {
    names(x@qc_list)
})
