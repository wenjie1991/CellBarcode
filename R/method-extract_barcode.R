## private functions
###########################

bc_process_sample_name <- function(sample_name, metadata, input_names) {
    # input metadata,
    # input samle_name,
    # names of input data or the file name
    # metadata = sample_name > names or file name
    if (is.null(sample_name)) {
        if (!is.null(metadata)) {
            sample_name <- rownames(metadata)
        } else {
            sample_name <- input_names
        }
    } else {
        if (length(sample_name) != length(input_names)) {
            stop("Sample name length does not meet sample number.")
        }
    }
        
    if (!is.null(metadata)) {
        if (!all(sample_name == rownames(metadata))) {
            stop("Sample name does not match row name of metadata.")
        }
    }
    sample_name
}

bc_process_metadata <- function(sample_name, old_metadata, new_metadata) {
    res <- merge(old_metadata, new_metadata, by = 0, all=TRUE)
    res <- as.data.frame(res)
    rownames(res) <- res$Row.names
    res$Row.names <- NULL
    res[sample_name, ]
}

bc_extract_metadata <- function(x, sample_name) {
    d <- vapply(x, function(x_i) {
        res <- c(attr(x_i, "raw_read_count"), attr(x_i, "barcode_read_count"))
        res
    }, c("raw_read_count" = 1, "barcode_read_count" = 1))
    d <- as.data.frame(t(d))
    rownames(d) <- sample_name
    d
}

#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "data.frame", function(
    x, 
    pattern = "", 
    sample_name = NULL, 
    maxLDist = 0, 
    pattern_type = c(barcode = 1), 
    costs = list(sub = 1, ins = 99, del = 99), 
    ordered = TRUE) {

    # sequence frequecy
    reads_freq <- x$freq
    # sequence
    reads_seq <- x$seq

    # if maxLDist > 0, use utils::aregexec to match the barcode else use
    # stringr::str_match to perform the barcode match. The stringr::str_match is
    # faster
    if (maxLDist > 0) {
        # regular expression match
        m <- utils::aregexec(
            pattern, 
            reads_seq, 
            fixed = FALSE, 
            max.distance = maxLDist, 
            costs = costs
        )

        # captured sequence
        seq_l <- regmatches(reads_seq, m)

        # matched sequence
        seq_v <- plyr::laply(seq_l, function(x) { x[1] })
        # captured barcode
        barcode_v <- plyr::laply(seq_l, function(x) { 
            x[pattern_type["barcode"] + 1] 
        })

        if ("UMI" %in% names(pattern_type)) {
            umi_v <- plyr::laply(seq_l, function(x) { x[pattern_type["UMI"] + 1] })
        }

    } else {
        # regular expression match
        m <- stringr::str_match(reads_seq, pattern)

        # capture sequence
        seq_v <- m[, 1]

        # matched barcode
        barcode_v <- m[, pattern_type['barcode'] + 1]

        # matched UMI
        if ("UMI" %in% names(pattern_type)) {
            umi_v <- m[, pattern_type['UMI'] + 1]
        }
    }

    # captured UMI
    if ("UMI" %in% names(pattern_type)) {
        d <- data.table(
            reads_seq = reads_seq, 
            match_seq = seq_v, 
            umi_seq = umi_v, 
            barcode_seq = barcode_v, 
            count = reads_freq
        )
    } else {
        d <- data.table(
            reads_seq = reads_seq, 
            match_seq = seq_v, 
            barcode_seq = barcode_v, 
            count = reads_freq
        ) 
    }

    # order the data by counts
    if (ordered) {
        d <- d[order(count, decreasing = TRUE)]
    } 


    d <- stats::na.omit(d)

    # prepare the row_read_count & barcode_read_count
    attr(d, "raw_read_count") <- sum(reads_freq)
    attr(d, "barcode_read_count") <- sum(d$count)

    d
})


#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "ShortReadQ", function(
    x, 
    pattern = "", 
    sample_name = NULL, 
    maxLDist = 0, 
    pattern_type = c(barcode = 1), 
    costs = list(sub = 1, ins = 99, del = 99), 
    ordered = TRUE) {

    # sequence frequecy
    reads_freq <- ShortRead::tables(x, n=Inf)$top

    x = data.frame(
        freq = as.integer(reads_freq),
        seq = names(reads_freq)
    )

    bc_extract(x, 
        pattern = pattern, 
        maxLDist = maxLDist, 
        pattern_type = pattern_type, 
        costs = costs, 
        ordered = ordered
    )
})


#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "DNAStringSet", function(
    x, 
    pattern = "", 
    sample_name = NULL, 
    maxLDist = 0, 
    pattern_type = c(barcode = 1), 
    costs = list(sub = 1, ins = 99, del = 99), 
    ordered = TRUE) {

    # sequence frequecy
    reads_freq <- ShortRead::tables(x, n=Inf)$top

    x = data.frame(
        freq = as.integer(reads_freq),
        seq = names(reads_freq)
    )


    bc_extract(
        x, 
        pattern = pattern, 
        sample_name = sample_name, 
        maxLDist = maxLDist, 
        pattern_type = pattern_type, 
        costs = costs, 
        ordered = ordered
    )
})

#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "integer", function(
    x,
    pattern = "",
    sample_name = NULL,
    maxLDist = 0,
    pattern_type = c(barcode = 1),
    costs = list(sub = 1, ins = 99, del = 99),
    ordered = TRUE) {

    x = data.frame(
        freq = as.integer(x),
        seq = names(x)
    )

    bc_extract(
        x,
        pattern = pattern,
        maxLDist = maxLDist,
        pattern_type = pattern_type,
        costs = costs,
        ordered = ordered)
})


#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "character", function(
    x, 
    pattern = "", 
    sample_name = NULL,
    metadata = NULL, 
    maxLDist = 0, 
    pattern_type = c(barcode = 1), 
    costs = list(sub = 1, ins = 99, del = 99), 
    ordered = TRUE) {

    # sample_name given
    # meta_data given
    #   no sample_name
    #   no rowname
    
    if (!is.null(metadata) & !is(metadata, "data.frame")) {
        stop("metadata should be a data.frame.")
    }

    # if more than one fastq file as input
    if (length(x) > 1) {

        input_names <- basename(x)
        sample_name <- bc_process_sample_name(sample_name, metadata, input_names)

        messyBc <- lapply(seq_along(x), function(i) {
            bc_extract(
                ShortRead::readFastq(x[i]), 
                sample_name = sample_name[i], 
                pattern = pattern, 
                maxLDist = maxLDist, 
                pattern_type = pattern_type, 
                costs = costs, 
                ordered = ordered)
        })

        new_metadata <- bc_extract_metadata(messyBc, sample_name)
        metadata <- bc_process_metadata(sample_name, metadata, new_metadata)

        names(messyBc) <- sample_name
        # output <- list(messyBc = messyBc, metadata = metadata)
        # class(output) <- "BarcodeObj"
        output <- BarcodeObj(metadata=metadata, messyBc=messyBc)
        return(output)
    } else {
        # if one fastq file as input
        bc_extract(
            ShortRead::readFastq(x), 
            sample_name = sample_name[i], 
            pattern = pattern,
            maxLDist = maxLDist,
            pattern_type = pattern_type,
            costs = costs,
            ordered = ordered)
    }
})


#' @rdname bc_extract
#' @exportMethod bc_extract
setMethod("bc_extract", "list", function(
    x,
    pattern = "",
    sample_name = NULL,
    metadata = NULL,
    maxLDist = 0,
    pattern_type = c(barcode = 1),
    costs = list(sub = 1, ins = 99, del = 99),
    ordered = TRUE) {

    input_names <- ifnullelse(names(x), seq_along(x))
    sample_name <- bc_process_sample_name(sample_name, metadata, input_names)

    lapply(seq_along(x), function(i) {
        bc_extract(
            x[[i]],
            pattern = pattern,
            sample_name = sample_name[i],
            maxLDist = maxLDist,
            pattern_type = pattern_type,
            costs = costs,
            ordered = ordered)
    }) -> messyBc

    new_metadata <- bc_extract_metadata(messyBc, sample_name)
    metadata <- bc_process_metadata(sample_name, metadata, new_metadata)

    names(messyBc) <- sample_name

    # output <- list(messyBc = messyBc, metadata = metadata)
    # class(output) <- "BarcodeObj"
    output <- BarcodeObj(metadata=metadata, messyBc=messyBc)
    output
})
