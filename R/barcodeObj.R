check_sample_name <- function(barcodeObj) {
  meta_name <- rownames(barcodeObj$meta_data)
  if (all(meta_name != names(barcodeObj$messyBc))) stop("The messyBc sample names are not consistent with metadta")
  if (!is.null(barcodeObj$cleanBc)) {
    if (all(meta_name != names(barcodeObj$cleanBc))) stop("The cleanBc sample names are not consistent with metadta")
  }
}


#' Subset operation of BarcodeObj object
#'
#' @param barcodeObj A BarcodeObj
#' @param sample A vector of integer or string, select the barcode subset from BarcodeObj
#' @param barcode A vector or an expression, select the sample subset from BarcodeObj
#' @param black_list A character vector. Remove the the barcode in the list.
#' @param white_list A character vector. The barcode reference list.
#' @param is_sample_quoted_exp A bool value. If the sample parameter is quoted expression or not
#' @param barcodeObj_x A BarcodeObj
#' @param barcodeObj_y A BarcodeObj
#' @param ... Additional arguments
#' @return A BarcodeObj
#'
#' @examples
#' data(bc_obj)
#'
#' bc_obj
#'
#' # Select barcodes
#' bc_subset(bc_obj, barcode = c("AACCTT", "AACCTT"))
#' bc_obj[c("AGAG", "AAAG"), ]
#'
#' # Select by meta data
#' bc_meta(bc_obj)$phenotype <- c("l", "b")
#' bc_meta(bc_obj)
#' bc_subset(bc_obj, phenotype == "l")
#'
#' bc_obj[, phenotype == "l"]
#'
#' # Select samples by sample name
#' bc_obj[, "test1"]
#' bc_obj[, c("test1", "test2")]
#' bc_subset(bc_obj, sample = "sample1_rep1", barcode = c("AACCTT", "AACCTT"))
#' 
#' # Apply black list
#' bc_subset(bc_obj, sample = c("sample1_rep1", "sample1_rep2"), barcode = c("AACCTT"))
#'
#' # Join two samples with different barcode sets
#' bc_obj["AGAG", "test1"] + bc_obj["AAAG", "test1"]
#'
#' # Join two samples with overlap barcodes
#' bc_obj_join <- bc_obj["AGAG", "test1"] + bc_obj["AGAG", "test1"]
#' bc_obj_join
#' # The same barcode will merged after applying bc_cure()
#' bc_cure(bc_obj_join)
#'
#' # Remove barcodes
#' bc_obj
#' bc_obj - "AAAG"
#' 
#' # Select barcode in white list
#' bc_obj
#' bc_obj * "AAAG"
#' @export bc_subset
#' @rdname bc_subset
bc_subset <- function(barcodeObj, sample = NULL, barcode = NULL, black_list = NULL, is_sample_quoted_exp = FALSE) {

  reads_eq <- NULL # due to NOTE in check

  meta_data <- barcodeObj$meta_data
  if (is_sample_quoted_exp) {
    sample_call <- sample
  } else {
    sample_call <- substitute(sample)
  }

  # TODO: The funciton only can apply the operation to the `messyBc` and `cleanBc`. We need to make it
  # capable to apply the selection to all information in the object.

  ## TODO: How do handle messyBc
  if (!is.null(barcode)) {
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <- lapply(barcodeObj$messyBc, function(d) { d[barcode_seq %in% barcode] })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- lapply(barcodeObj$cleanBc, function(d) { d[barcode_seq %in% barcode] })
    }
  }
  if (!is.null(black_list)) {
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <- lapply(barcodeObj$messyBc, function(d) { d[!(barcode_seq %in% black_list)] })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- lapply(barcodeObj$cleanBc, function(d) { d[!(barcode_seq %in% black_list)] })
    }
  }
  if (!is.null(sample_call)) {
    check_sample_name(barcodeObj)

    sample_i <- eval(sample_call, meta_data, parent.frame())

    barcodeObj$meta_data <- meta_data[sample_i, , drop=FALSE]

    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <-  barcodeObj$messyBc[rownames(barcodeObj$meta_data)]
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- barcodeObj$cleanBc[rownames(barcodeObj$meta_data)]
    }
  }
  return(barcodeObj)
}

#' @export
#' @rdname bc_subset
"[.BarcodeObj" <- function(barcodeObj, barcode = NULL, sample = NULL, ...) {

  y_call <- substitute(sample)
  return(bc_subset(barcodeObj, sample = y_call, barcode = barcode, is_sample_quoted_exp = TRUE))
}

#' @rdname bc_subset
#' @export
"+.BarcodeObj" <- function(barcodeObj_x, barcodeObj_y) {
  # TODO: Apply the merge to all parts of the data
  #       How to deal when two BarcodeObj have the same samples, the same samples will merged

  x <- barcodeObj_x
  y <- barcodeObj_y

  meta_data_x <- x$meta_data
  meta_data_y <- y$meta_data

  meta_data_xy <- merge(meta_data_x, meta_data_y, by = "sample_name", all = TRUE)
  rownames(meta_data_xy) <- meta_data_xy$sample_name

  x$messyBc <- lapply(rownames(meta_data_xy), function(sample_name) {
    rbind(x$messyBc[[sample_name]], y$messyBc[[sample_name]])
  })
  names(x$messyBc) <- meta_data_xy$sample_name

  if (!is.null(x$cleanBc) & !is.null(y$cleanBc)) {

    x$cleanBc <- lapply(rownames(meta_data_xy), function(sample_name) {
      rbind(x$cleanBc[[sample_name]], y$cleanBc[[sample_name]])
    })
    names(x$cleanBc) <- meta_data_xy$sample_name
  }
  x$meta_data <- meta_data_xy
  x
}

#' @rdname bc_subset
#' @export
"-.BarcodeObj" <- function(barcodeObj, black_list) {

  x <- barcodeObj
  y <- black_list

  bc_subset(x, black_list = y)
}

#' @rdname bc_subset
#' @export
"*.BarcodeObj" <- function(barcodeObj, white_list) {

  x <- barcodeObj
  y <- white_list

  bc_subset(x, barcode = y)
}

#' Get barcode sequences
#' 
#' @param barcodeObj A BarcodeObj.
#' @param unlist A bool value, if TRUE, then the function will return a vector of barcodes; else a list will be returned, each element is corresponding to a sample.
#' @return A character vector or a list
#' @examples
#' data(bc_obj)
#'
#' bc_barcodes(bc_obj)
#'
#' bc_barcodes(bc_obj, unlist = FALSE)
#' @export
bc_barcodes <- function(barcodeObj, unlist = TRUE) {
  d <- lapply(barcodeObj$cleanBc, function(x) { x$barcode_seq }) 
  if (unlist) {
    unlist(d) %>% as.character %>% unique
  } else {
    names(d) = names(barcodeObj$cleanBc)
    d
  }
}

#' Get sample names
#'
#' @param barcodeObj A BarcodeObj
#' @param value A string vector, the new sample names
#' @return A character vector
#' @examples
#' data(bc_obj)
#'
#' bc_names(bc_obj)
#' bc_names(bc_obj) <- c("new1", "new2")
#' @export
bc_names <- function(barcodeObj) {
  check_sample_name(barcodeObj)
  rownames(barcodeObj$meta_data)
}

#' @rdname bc_names
#' @export
"bc_names<-" <- function(barcodeObj, value) {
  check_sample_name(barcodeObj)
  if (length(value) != nrow(barcodeObj$meta_data)) stop("The sample names do not have correct length.")
  if (!is.null(barcodeObj$messyBc)) {
    names(barcodeObj$messyBc) <- value
  } 
  if (!is.null(barcodeObj$cleanBc)) {
    names(barcodeObj$cleanBc) <- value
  }
  rownames(barcodeObj$meta_data) <- value
  barcodeObj$meta_data$sample_name <- value
  barcodeObj
}

#' Get meta data
#'
#' @param barcodeObj A BarcodeObj
#' @param key A string, the meta data record name
#' @param value A string vector, the meta data record value
#' @return A BarcodeObj
#' @examples
#' data(bc_obj)
#'
#' bc_meta(bc_obj)
#' bc_meta(bc_obj)$phenotype <- c("l", "b")
#' bc_meta(bc_obj, key = "sample_type") <- c("l", "b")
#' bc_meta(bc_obj)
#' meta_data <- data.frame(
#'   sample_name <- c("test1", "test2"),
#'   phenotype <- c("l", "b")
#'   )
#' bc_meta(bc_obj) <- meta_data
#' bc_meta
#' @export
bc_meta <- function(barcodeObj) {
  check_sample_name(barcodeObj)
  barcodeObj$meta_data
}

#' @rdname bc_meta
#' @export
"bc_meta<-" <- function(barcodeObj, key = NULL, value) {
  check_sample_name(barcodeObj)
  if (is.null(key)) {
    if (nrow(value) != nrow(barcodeObj$meta_data)) stop("The given meta data does not hvae correct length.")
    barcodeObj$meta_data <- value
  } else {
    if (length(value) != nrow(barcodeObj$meta_data)) stop("The given meta data does not hvae correct length.")
    barcodeObj$meta_data[[key]] <- value
  }
  barcodeObj
}


#' Transform BarcodeObj into other data type
#' 
#' @param barcodeObj A BarcodeObj
#' @return data.frame
#' @examples
#' data(bc_obj)
#'
#' bc_obj <- bc_cure(bc_obj)
#'
#' bc_2df(bc_obj)
#' bc_2dt(bc_obj)
#' bc_2matrix(bc_obj)
#' @export
bc_2df <- function(barcodeObj) {
  bc_2dt(barcodeObj) %>% as.data.frame
}

#' @rdname bc_2df
#' @return data.table
#' @export
bc_2dt <- function(barcodeObj) {
  # TODO: enable messyBc
  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found.")
  }
  d <- barcodeObj$cleanBc %>% rbindlist(idcol = TRUE)
  names(d)[1] <- "sample_name"
  return(d)
}


#' @rdname bc_2df
#' @return matrix
#' @export
bc_2matrix <- function(barcodeObj) {
  # TODO: enable messyBc
  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found.")
  }
  # TODO: use sparse matrix
  d <- bc_2dt(barcodeObj) %>% data.table::dcast(barcode_seq ~ sample_name, value.var = "count", fill = 0) 
  m <- data.frame(d[, -1]) %>% data.matrix
  rownames(m) = d[, barcode_seq]
  m
} 

summary_BarcodeObj <- function(barcodeObj) {

  summary_res <- list()
  summary_res$messyBc_barcode_n <- lapply(barcodeObj$messyBc, nrow)
  names(summary_res$messyBc_barcode_n) <- rownames(barcodeObj$meta_data)
  if (!is.null(barcodeObj$cleanBc)) {
    summary_res$cleanBc_barcode_n <- lapply(barcodeObj$cleanBc, nrow)
    names(summary_res$cleanBc_barcode_n) <- rownames(barcodeObj$meta_data)
  }

  summary_res
}

#' Format the summary of BarcodeObj
#'
#' @param x A BarcodeObj
#' @param ... Additional arguments
#'
#' @examples
#' data(bc_obj)
#'
#' format(bc_obj)
#' print(bc_obj)
#' @export
format.BarcodeObj <- function(x, ...) {

  summary_res <- summary_BarcodeObj(x)

  subjects <- paste(names(x), collapse = "  ")
  messyBc_n <- length(summary_res$messyBc_barcode_n)
  messyBc_info <- lapply(seq_along(summary_res$messyBc_barcode_n), function(i) { 
    sample_name <- names(summary_res$messyBc_barcode_n)[i]
    n <- summary_res$messyBc_barcode_n[i]
    stringr::str_glue("    In sample ${sample_name} there are: {n} Tags")
  }) %>% unlist %>% paste(collapse = "\n")

  meta_data_info <- paste(names(x$meta_data), collapse = "  ")
  meta_data_n <- length(names(x$meta_data))

  
  res <- stringr::str_glue(
"Bonjour le monde. This is a baby barcode Object.
----------
It contains: 
    {subjects}
----------
$meta_data: {meta_data_n} meta data field(s) available
    {meta_data_info}
----------
$messyBc: {messyBc_n} Samples for raw barcodes
{messyBc_info}")
  
  if (!is.null(summary_res$cleanBc_barcode_n)) {

    cleanBc_n <- length(summary_res$cleanBc_barcode_n)
    cleanBc_info <- lapply(seq_along(summary_res$cleanBc_barcode_n), function(i) { 
      sample_name <- names(summary_res$cleanBc_barcode_n)[i]
      n <- summary_res$cleanBc_barcode_n[i]
      stringr::str_glue("    In sample ${sample_name} there are: {n} barcodes")
    }) %>% unlist %>% paste(collapse = "\n")
    res_cleanBc <- stringr::str_glue(
"----------
$cleanBc: {cleanBc_n} Samples for cleaned barcodes
{cleanBc_info}")
    res <- stringr::str_c(res, "\n", res_cleanBc)
  }
  res
}

#' Print the summary of BarcodeObj
#'
#' @param x A BarcodeObj
#' @param ... Additional arguments
#'
#' @examples
#' data(bc_obj)
#'
#' format(bc_obj)
#' print(bc_obj)
#' @export
print.BarcodeObj <- function(x, ...) {
  cat(format(x), "\n")
  #   cat("Bonjour le monde. This is a baby barcode Object.\n----------\n")
  #   cat("It contains:\n")
  ## The items in the list
  #   cat(names(barcodeObj), collapse = "  ", "\n----------\n")
  #   if (!is.null(barcodeObj$messyBc)) {
  ## How many samples in messyBc item
  #     cat("$messyBc:", length(barcodeObj$messyBc), "Samples for uncleaned barcodes\n")
  #     for (i in names(barcodeObj$messyBc)) {
  ## The number of barcode in each sample
  #       cat("    In sample", paste0("$", i), "there are:", nrow(barcodeObj$messyBc[[i]]), "Tags\n")
  #     }
  #     cat("\n")
  #   }
  #   if (!is.null(barcodeObj$cleanBc)) {
  ## How many samples in cleanBc iterm
  #     cat("$cleanBc: ", length(barcodeObj$cleanBc), "Samples for cleaned barcodes\n")
  #     for (i in names(barcodeObj$cleanBc)) {
  ## The number of barcode in each sample
  #       cat("    In sample", paste0("$", i), "there are:", nrow(barcodeObj$cleanBc[[i]]), "barcodes\n")
  #     }
  #     cat("\n")
  #   }
}
