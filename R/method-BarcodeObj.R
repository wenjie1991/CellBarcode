## internal functions
###########################
count_BarcodeObj <- function(barcodeObj) {
    summary_res <- list()

    # number of barcodes for each samples in messyBc
    summary_res$messyBc_barcode_n <- lapply(barcodeObj@messyBc, nrow)
    names(summary_res$messyBc_barcode_n) <- rownames(barcodeObj@metadata)

    # number of barcodes for each samples in cleanBc
    if (!is.null(barcodeObj@cleanBc)) {
        summary_res$cleanBc_barcode_n <- lapply(barcodeObj@cleanBc, nrow)
        names(summary_res$cleanBc_barcode_n) <- rownames(barcodeObj@metadata)
    }

    summary_res
}

#' @rdname bc_subset
#' @exportMethod bc_subset
setMethod("bc_subset", c("BarcodeObj"), function(
    barcodeObj,
    sample = NULL,
    barcode = NULL,
    black_list = NULL,
    is_sample_quoted_exp = FALSE) {

    # check sample names consistency
    # check_sample_name(barcodeObj)

    metadata <- barcodeObj@metadata

    # preprocess sample selecting expression
    if (is_sample_quoted_exp) {
        sample_call <- sample
    } else {
        sample_call <- substitute(sample)
    }

    # TODO: The function only can apply the operation to the `messyBc` `cleanBc`
    # and `metadata`. We need to make it capable to apply the selection to all
    # information in the object.

    ## TODO: How do handle messyBc

    # select barcodes
    if (!is.null(barcode)) {
        # select barcodes in messyBc
        if (!is.null(barcodeObj@messyBc)) {
            barcodeObj@messyBc <- lapply(barcodeObj@messyBc, function(d) {
                d[barcode_seq %in% barcode]
            })
        }
        # select barcodes in cleanBc
        if (!is.null(barcodeObj@cleanBc)) {
            barcodeObj@cleanBc <- lapply(barcodeObj@cleanBc, function(d) {
                d[barcode_seq %in% barcode]
            })
        }
    }

    # remove barcodes
    if (!is.null(black_list)) {
        # remove barcodes in messyBc
        if (!is.null(barcodeObj@messyBc)) {
            barcodeObj@messyBc <- lapply(barcodeObj@messyBc, function(d) {
                d[!(barcode_seq %in% black_list)]
            })
        }

        # remove barcodes in cleanBc
        if (!is.null(barcodeObj@cleanBc)) {
            barcodeObj@cleanBc <- lapply(barcodeObj@cleanBc, function(d) {
                d[!(barcode_seq %in% black_list)]
            })
        }
    }

    # select samples
    if (!is.null(sample_call)) {
        
        # evaluate the sample argument
        sample_i <- eval(sample_call, metadata, parent.frame())

        # subset metadata
        barcodeObj@metadata <- metadata[sample_i, , drop = FALSE]

        # subset messyBc
        if (!is.null(barcodeObj@messyBc)) {
            barcodeObj@messyBc <-
                barcodeObj@messyBc[rownames(barcodeObj@metadata)]
        }

        # subset cleanBc
        if (!is.null(barcodeObj@cleanBc)) {
            barcodeObj@cleanBc <-
                barcodeObj@cleanBc[rownames(barcodeObj@metadata)]
        }
    }


    # check_sample_name(barcodeObj)
    return(barcodeObj)
})


#' @export
"[.BarcodeObj" <- function(barcodeObj, barcode = NULL, sample = NULL, ...) {
    # do not evaluate the expression
    y_call <- sample
    # invoke bc_subset to done the job
    return(
        bc_subset(
            barcodeObj,
            sample = y_call,
            barcode = barcode,
            is_sample_quoted_exp = TRUE
        )
    )
}

#' @rdname bc_subset
#' @exportMethod bc_merge
setMethod("bc_merge", c("BarcodeObj", "BarcodeObj"), function(barcodeObj_x, barcodeObj_y) {
    barcodeObj_x + barcodeObj_y
})

#' @rdname bc_subset
#' @export
"+.BarcodeObj" <- function(barcodeObj_x, barcodeObj_y) {
    # TODO: Apply the merge to all parts of the data
    #       How to deal when two BarcodeObj have the same samples, the same
    #       samples will merged

    # merge metadata
    suffixes <- paste0(".", 
        c(
        deparse(substitute(barcodeObj_x)), 
        deparse(substitute(barcodeObj_y))
        ))
    metadata_x <- barcodeObj_x@metadata
    metadata_y <- barcodeObj_y@metadata
    metadata_xy <-
        merge(metadata_x, metadata_y, by = 0, all = TRUE, suffixes=suffixes, no.dups = TRUE)
    rownames(metadata_xy) <- metadata_xy$Row.names
    metadata_xy$Row.names <- NULL

    # merge messyBc
    # if the messyBc do not have the same header, do not merge them
    flag_remove_umi <- FALSE
    if (!all(names(barcodeObj_x@messyBc) == names(barcodeObj_y@messyBc))) {
        message("------------\n+.BarcodeObj: You are merge data with UMI to data without UMI. The UMI info are discarded.\n------------")
        flag_remove_umi <- TRUE
    }

    barcodeObj_x@messyBc <-
        lapply(rownames(metadata_xy), function(sample_name) {
            d_x <- barcodeObj_x@messyBc[[sample_name]]
            d_y <- barcodeObj_y@messyBc[[sample_name]]
            if (flag_remove_umi) {
                if (!is.null(d_x$umi_seq))
                    d_x$umi_seq <- NULL
                if (!is.null(d_y$umi_seq))
                    d_y$umi_seq <- NULL
            }

            d_merged <- rbind(d_x, d_y)
            var_by <- setdiff(names(d_merged), "count")
            d_merged[, .(count = sum(count)), by = var_by]
        })

    names(barcodeObj_x@messyBc) <- rownames(metadata_xy)

    # merge cleanBc
    if (!is.null(barcodeObj_x@cleanBc) & !is.null(barcodeObj_y@cleanBc)) {
        barcodeObj_x@cleanBc <-
            lapply(rownames(metadata_xy), function(sample_name) {
                d_x <- barcodeObj_x@cleanBc[[sample_name]]
                d_y <- barcodeObj_y@cleanBc[[sample_name]]

                d_merged <- rbind(d_x, d_y)
                var_by <- setdiff(names(d_merged), "count")
                d_merged[, .(count = sum(count)), by = var_by]
            })

        names(barcodeObj_x@cleanBc) <- rownames(metadata_xy)
    } else if (is.null(barcodeObj_x@cleanBc) + is.null(barcodeObj_y@cleanBc) == 1) {
        message("------------\n+.BarcodeObj: One of the BarcodesObj does not have cleanBc, discard the cleanBc while merging.\n------------")
        barcodeObj_x@cleanBc <- NULL
        bc_meta(barcodeObj_x, "depth_cutoff") <- NULL
    }

    barcodeObj_x@metadata <- metadata_xy
    barcodeObj_x
}

#' @rdname bc_subset
#' @export
"-.BarcodeObj" <- function(barcodeObj, black_list) {
    bc_subset(barcodeObj, black_list = black_list)
}

#' @rdname bc_subset
#' @export
"*.BarcodeObj" <- function(barcodeObj, white_list) {
    bc_subset(barcodeObj, barcode = white_list)
}

#' @rdname bc_barcodes
#' @exportMethod bc_barcodes
setMethod("bc_barcodes", "BarcodeObj", function(barcodeObj, unlist = TRUE) {
    if (is.null(barcodeObj@cleanBc)) {
        stop("No cleanBc found, please run bc_cure_* first.")
    }

    d <- lapply(barcodeObj@cleanBc, function(x) {
        x$barcode_seq
        })
    if (unlist) {
        unlist(d) %>% as.character %>% unique
    } else {
        names(d) = names(barcodeObj@cleanBc)
        d
    }
})

#' @rdname bc_names
#' @exportMethod bc_names
setMethod("bc_names", c("BarcodeObj"), function(x){
    rownames(x@metadata)
})


#' @rdname bc_names
#' @exportMethod bc_names<-
setMethod("bc_names<-", c("BarcodeObj", "character"), function(x, value) {
    # if the sample names are not consistent stop
    # check_sample_name(barcodeObj)

    # check if the new names fit the sample number
    if (length(unique(value)) != nrow(x@metadata))
        stop("The given sample names do not have the same length with sample number. Or the sample names are not unique.")

    # If exists messyBc renew the name
    if (!is.null(x@messyBc)) 
        names(x@messyBc) <- value
    
    # If exists cleanBc renew the name
    if (!is.null(x@cleanBc)) 
        names(x@cleanBc) <- value

    # renew sample name in metadata
    rownames(x@metadata) <- value
    # barcodeObj@metadata <- value
    # check_sample_name(barcodeObj)
    x
})

#' @rdname bc_meta
#' @exportMethod bc_meta
setMethod("bc_meta", c("BarcodeObj"), function(barcodeObj) {
    check_sample_name(barcodeObj)
    barcodeObj@metadata
})

#' @rdname bc_meta
#' @exportMethod bc_meta<-
setMethod("bc_meta<-", c("BarcodeObj", "ANY", "ANY"), function(barcodeObj, key = NULL, value) {
    # check sample names consistency
    check_sample_name(barcodeObj)

    # if no key is given, update the metadata
    if (is.null(key)) {
        if (!is(value, "data.frame"))
            stop("The input data is not data.frame")

        # if new value matches the sample number
        if (length(value) != 1 & nrow(value) != nrow(barcodeObj@metadata))
            stop("The given meta data does not have correct length.")

        barcodeObj@metadata <- value
    } else {
        # if new value matches the sample number
        if (length(value) != 1 & length(value) != nrow(barcodeObj@metadata))
            stop("The given meta data does not have correct length.")

        barcodeObj@metadata[[key]] <- value
    }

    # check sample names consistency
    check_sample_name(barcodeObj)

    barcodeObj
})

#' @rdname bc_2df
#' @exportMethod bc_2df
setMethod("bc_2df", c("BarcodeObj"), function(barcodeObj) {
    bc_2dt(barcodeObj) %>% as.data.frame
})


#' @rdname bc_2df
#' @exportMethod bc_2dt
setMethod("bc_2dt", c("BarcodeObj"), function(barcodeObj) {
    # TODO: enable messyBc?

    if (is.null(barcodeObj@cleanBc)) {
        stop("No cleanBc found, please run bc_cure_* first.")
    }

    d <- barcodeObj@cleanBc %>% rbindlist(idcol = TRUE)
    names(d)[1] <- "sample_name"
    return(d)
})

#' @rdname bc_2df
#' @exportMethod bc_2matrix
setMethod("bc_2matrix", c("BarcodeObj"), function(barcodeObj) {
    # TODO: enable messyBc?

    if (is.null(barcodeObj@cleanBc)) {
        stop("No cleanBc found, please run bc_cure_* first.")
    }

    # TODO: use sparse matrix?
    # long shape data to wide shape data
    # with barcode_seq by row and sample by column
    d <- bc_2dt(barcodeObj) %>% data.table::dcast(barcode_seq ~ sample_name,
        value.var = "count",
        fill = 0)
    m <- data.frame(d[,-1]) %>% data.matrix
    rownames(m) = d[, barcode_seq]
    m
})

#' Formats BarcodeObj object
#'
#' Format the summary of BarcodeObj object for pretty print.
#'
#' @param x A BarcodeObj object
#' @return Formated summary text.
#'
#' @examples
#' data(bc_obj)
#'
#' # format BarcodeObj for pretty print
#' format(bc_obj)
#'
#' ###
#' @rdname format
#' @exportMethod format
setMethod("format", c("BarcodeObj"), function(x) {
    summary_res <- count_BarcodeObj(x)

    # elements list in BarcodeObj
    # subjects <- paste(slotNames(barcodeObj), collapse = "  ")

    # number of samples
    messyBc_n <- length(summary_res$messyBc_barcode_n)

    messyBc_info <-
        lapply(seq_along(summary_res$messyBc_barcode_n), function(i) {
            # sample name
            sample_name <- names(summary_res$messyBc_barcode_n)[i]
            # barcode number
            n <- summary_res$messyBc_barcode_n[i]
            stringr::str_glue("    In sample ${sample_name} there are: {n} Tags")
            }) %>% unlist %>% paste(collapse = "\n")

    # elements in metadata
    metadata_info <- paste(names(x@metadata), collapse = "  ")
    metadata_n <- length(names(x@metadata))


    res <- stringr::str_glue(
"Bonjour le monde, This is a BarcodeObj.
----------
It contains: 
----------
@metadata: {metadata_n} field(s) available:
{metadata_info}
----------
@messyBc: {messyBc_n} sample(s) for raw barcodes:
{messyBc_info}"
    )


    # cleanBc
    if (!is.null(summary_res$cleanBc_barcode_n)) {
        # sample number
        cleanBc_n <- length(summary_res$cleanBc_barcode_n)
        cleanBc_info <-
            lapply(seq_along(summary_res$cleanBc_barcode_n),
                function(i) {
                    # sample name
                    sample_name <-
                        names(summary_res$cleanBc_barcode_n)[i]
                    # sample number
                    n <- summary_res$cleanBc_barcode_n[i]
                    stringr::str_glue("    In sample ${sample_name} there are: {n} barcodes")
                }) %>% unlist %>% paste(collapse = "\n")
        res_cleanBc <- stringr::str_glue(
"----------
@cleanBc: {cleanBc_n} samples for cleaned barcodes
{cleanBc_info}"
        )
        res <- stringr::str_c(res, "\n", res_cleanBc)
    }
    res
})

#' Show BarcodeObj object
#'
#' Show the summary of BarcodeObj object for pretty print.
#'
#' @param object A BarcodeObj or a BarcodeQcSet object
#' @return Formated summary text.
#'
#' @examples
#' data(bc_obj)
#'
#' # show BarcodeObj for pretty print
#' bc_obj
#'
#' ###
#' @rdname show
#' @exportMethod show
setMethod("show", c("BarcodeObj"), function(object) {
    cat(format(object), "\n")
})



