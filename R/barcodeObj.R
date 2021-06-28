check_sample_name <- function(barcodeObj) {

  meta_name <- rownames(barcodeObj$metadata)

  # metadata v.s. messyBc
  if (all(meta_name != names(barcodeObj$messyBc))) 
    stop("The messyBc sample names are not consistent with metadta")

  # metadata v.s. cleanBc
  if (!is.null(barcodeObj$cleanBc)) {
    if (all(meta_name != names(barcodeObj$cleanBc))) 
      stop("The cleanBc sample names are not consistent with metadta")
  }
}

#' Subset operation of BarcodeObj object
#'
#' The subset operators perform the sample and barcode management with
#' BarcodeObj.
#'
#' @param barcodeObj A BarcodeObj
#' @param barcode A vector of integer or string, select the barcode subset from
#' BarcodeObj
#' @param sample A vector or an expression, select the sample subset from
#' BarcodeObj
#' @param black_list A character vector. Remove the listed barcodes.
#' @param white_list A character vector. Only keep the listed barcodes.
#' @param barcodeObj_x A BarcodeObj
#' @param barcodeObj_y A BarcodeObj
#' @param is_sample_quoted_exp A bool value. If the sample parameter is quoted
#' expression or not
#' @param ... Additional arguments
#' @return A BarcodeObj
#'
#' @details
#' `bc_subset` and `[]` used for selecting samples barcodes.
#'
#' `+` used for combining two BarcodesObj
#'
#' `-` used for removing barcodes in a black_list from BarcodesObj
#'
#' `*` used for selecting barcodes in a white_list from BarcodesObj
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
#' bc_subset(
#'   bc_obj, 
#'   sample = c("sample1_rep1", "sample1_rep2"), 
#'   barcode = c("AACCTT"))
#'
#' # Join two samples with different barcode sets
#' bc_obj["AGAG", "test1"] + bc_obj["AAAG", "test1"]
#'
#' # Join two samples with overlap barcodes
#' bc_obj_join <- bc_obj["AGAG", "test1"] + bc_obj["AGAG", "test1"]
#' bc_obj_join
#' # The same barcode will merged after applying bc_cure_depth()
#' bc_cure_depth(bc_obj_join)
#'
#' # Remove barcodes
#' bc_obj
#' bc_obj - "AAAG"
#' 
#' # Select barcode in white list
#' bc_obj
#' bc_obj * "AAAG"
#' ###
#' @export bc_subset
#' @rdname bc_subset
bc_subset <- function(
  barcodeObj, 
  sample = NULL, 
  barcode = NULL, 
  black_list = NULL, 
  is_sample_quoted_exp = FALSE) {

  reads_eq <- NULL # due to NOTE in check

  # check sample names consistency
  check_sample_name(barcodeObj)


  metadata <- barcodeObj$metadata

  # preprocess sample selecting expression
  if (is_sample_quoted_exp) {
    sample_call <- sample
  } else {
    sample_call <- substitute(sample)
  }

  # TODO: The funciton only can apply the operation to the `messyBc` `cleanBc`
  # and `metadata`. We need to make it capable to apply the selection to all
  # information in the object.

  ## TODO: How do handle messyBc

  # select barcodes
  if (!is.null(barcode)) {
    # select barcodes in messyBc
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <- lapply(barcodeObj$messyBc, function(d) { 
        d[barcode_seq %in% barcode] 
      })
    }
    # select barcodes in cleanBc
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- lapply(barcodeObj$cleanBc, function(d) { 
        d[barcode_seq %in% barcode] 
      })
    }
  }

  # remove barcodes
  if (!is.null(black_list)) {

    # remove barcodes in messyBc
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <- lapply(barcodeObj$messyBc, function(d) { 
        d[!(barcode_seq %in% black_list)] 
      })
    }

    # remove barcodes in cleanBc
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- lapply(barcodeObj$cleanBc, function(d) { 
        d[!(barcode_seq %in% black_list)] 
      })
    }
  }

  # select samples
  if (!is.null(sample_call)) {

    # evaluate the sample argument
    sample_i <- eval(sample_call, metadata, parent.frame())

    # subset metadata
    barcodeObj$metadata <- metadata[sample_i, , drop=FALSE]

    # subset messyBc
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc <-  barcodeObj$messyBc[rownames(barcodeObj$metadata)]
    }

    # subset cleanBc
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc <- barcodeObj$cleanBc[rownames(barcodeObj$metadata)]
    }
  }

  return(barcodeObj)
}

#' @export
#' @rdname bc_subset
"[.BarcodeObj" <- function(barcodeObj, barcode = NULL, sample = NULL, ...) {

  # do not evaluate the expression
  y_call <- substitute(sample)
  # invoke bc_subset to done the job
  return(bc_subset(barcodeObj, sample = y_call, barcode = barcode, is_sample_quoted_exp = TRUE))
}

#' @rdname bc_subset
#' @export
"+.BarcodeObj" <- function(barcodeObj_x, barcodeObj_y) {
  # TODO: Apply the merge to all parts of the data
  #       How to deal when two BarcodeObj have the same samples, the same samples will merged

  # merge metadata
  metadata_x <- barcodeObj_x$metadata
  metadata_y <- barcodeObj_y$metadata
  metadata_xy <- merge(metadata_x, metadata_y, by = "sample_name", all = TRUE)
  rownames(metadata_xy) <- metadata_xy$sample_name

  # merge messyBc
  barcodeObj_x$messyBc <- lapply(rownames(metadata_xy), function(sample_name) {
    rbind(barcodeObj_x$messyBc[[sample_name]], barcodeObj_y$messyBc[[sample_name]])
  })
  names(barcodeObj_x$messyBc) <- metadata_xy$sample_name

  # merge cleanBc
  if (!is.null(barcodeObj_x$cleanBc) & !is.null(barcodeObj_y$cleanBc)) {
    barcodeObj_x$cleanBc <- lapply(rownames(metadata_xy), function(sample_name) {
      rbind(barcodeObj_x$cleanBc[[sample_name]], barcodeObj_y$cleanBc[[sample_name]])
    })
    names(barcodeObj_x$cleanBc) <- metadata_xy$sample_name
  }

  barcodeObj_x$metadata <- metadata_xy
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

#' Get barcode sequences
#' 
#' This function retrieve the barcode sequences from cleanBc element in
#' BarcodeObj.
#'
#' @param barcodeObj A BarcodeObj.
#' @param unlist A bool value, if TRUE, then the function will return a vector
#' of unique barcodes from all samples; otherwise a list will be returned, each
#' element is corresponding to a sample.
#' @return A character vector or a list
#' @examples
#' data(bc_obj)
#'
#' # get unique barcode vector from all samples
#' bc_barcodes(bc_obj)
#' 
#' # get barcode for each sample in a list
#' bc_barcodes(bc_obj, unlist = FALSE)
#' ###
#' @export
bc_barcodes <- function(barcodeObj, unlist = TRUE) {

  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found, please run bc_cure_* first.")
  }

  d <- lapply(barcodeObj$cleanBc, function(x) { x$barcode_seq }) 
  if (unlist) {
    unlist(d) %>% as.character %>% unique
  } else {
    names(d) = names(barcodeObj$cleanBc)
    d
  }
}

#' Access sample names in BarcodeObj
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
  rownames(barcodeObj$metadata)
}

#' @rdname bc_names
#' @export
"bc_names<-" <- function(barcodeObj, value) {

  # if the sample names are not consistent stop
  check_sample_name(barcodeObj)

  # check if the new names fit the sample number
  if (length(value) != nrow(barcodeObj$metadata)) 
    stop("The sample names do not have correct length.")

  # If exists messyBc renew the name
  if (!is.null(barcodeObj$messyBc)) {
    names(barcodeObj$messyBc) <- value
  } 

  # If exists cleanBc renew the name
  if (!is.null(barcodeObj$cleanBc)) {
    names(barcodeObj$cleanBc) <- value
  }

  # renew sample name in metadata
  rownames(barcodeObj$metadata) <- value
  barcodeObj$metadata$sample_name <- value

  barcodeObj
}

#' Access metadata in BarcodeObj
#'
#' @param barcodeObj A BarcodeObj
#' @param key A string, the meta data record name
#' @param value A string vector, the meta data record value
#' @return A BarcodeObj
#' @examples
#' data(bc_obj)
#' 
#' # get the metadata data.frame
#' bc_meta(bc_obj)
#'
#' # assign value to metadata by $ operation
#' bc_meta(bc_obj)$phenotype <- c("l", "b")
#'
#' # assign value to metasta by key argument
#' bc_meta(bc_obj, key = "sample_type") <- c("l", "b")
#' 
#' # show the updated metadata
#' bc_meta(bc_obj)
#'
#' # assign a new data.frame as metadata
#' metadata <- data.frame(
#'   sample_name <- c("test1", "test2"),
#'   phenotype <- c("l", "b")
#'   )
#' rownames(metadata) = bc_names(bc_obj)
#' bc_meta(bc_obj) <- metadata
#' ###
#' @export
bc_meta <- function(barcodeObj) {
  check_sample_name(barcodeObj)
  barcodeObj$metadata
}

#' @rdname bc_meta
#' @export
"bc_meta<-" <- function(barcodeObj, key = NULL, value) {

  # check sample names consistency
  check_sample_name(barcodeObj)

  # if no key is given, update the metadata
  if (is.null(key)) {
    if (!is(value, "data.frame"))
      stop("The input data is not data.frame")

    # if new value matches the sample number
    if (nrow(value) != nrow(barcodeObj$metadata)) 
      stop("The given meta data does not have correct length.")

    barcodeObj$metadata <- value
  } else {
    # if new value matches the sample number
    if (length(value) != nrow(barcodeObj$metadata)) 
      stop("The given meta data does not have correct length.")

    barcodeObj$metadata[[key]] <- value
  }

  # check sample names consistency
  check_sample_name(barcodeObj)

  barcodeObj
}


#' Transform BarcodeObj into other data type
#' 
#' @param barcodeObj A BarcodeObj
#' @return data.frame
#' @examples
#' data(bc_obj)
#'
#' bc_obj <- bc_cure_depth(bc_obj)
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
  # TODO: enable messyBc?

  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found, please run bc_cure_* first.")
  }

  d <- barcodeObj$cleanBc %>% rbindlist(idcol = TRUE)
  names(d)[1] <- "sample_name"
  return(d)
}


#' @rdname bc_2df
#' @return matrix
#' @export
bc_2matrix <- function(barcodeObj) {
  # TODO: enable messyBc?

  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found, please run bc_cure_* first.")
  }

  # TODO: use sparse matrix?
  # long shape data to wide shape data
  # with barcode_seq by row and sample by column
  d <- bc_2dt(barcodeObj) %>% 
    data.table::dcast(barcode_seq ~ sample_name, value.var = "count", fill = 0) 
  m <- data.frame(d[, -1]) %>% data.matrix
  rownames(m) = d[, barcode_seq]
  m
} 

count_BarcodeObj <- function(barcodeObj) {

  summary_res <- list()

  # number of barcodes for each samples in messyBc
  summary_res$messyBc_barcode_n <- lapply(barcodeObj$messyBc, nrow)
  names(summary_res$messyBc_barcode_n) <- rownames(barcodeObj$metadata)

  # number of barcodes for each samples in cleanBc
  if (!is.null(barcodeObj$cleanBc)) {
    summary_res$cleanBc_barcode_n <- lapply(barcodeObj$cleanBc, nrow)
    names(summary_res$cleanBc_barcode_n) <- rownames(barcodeObj$metadata)
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
#' # format BarcodeObj for pretty print
#' format(bc_obj)
#' ###
#' @seealso [Bc::print.BarcodeObj]
#' @export
format.BarcodeObj <- function(x, ...) {

  summary_res <- count_BarcodeObj(x)

  # elements list in BarcodeObj
  subjects <- paste(names(x), collapse = "  ")

  # number of samples
  messyBc_n <- length(summary_res$messyBc_barcode_n)

  messyBc_info <- lapply(seq_along(summary_res$messyBc_barcode_n), function(i) { 
    # sample name
    sample_name <- names(summary_res$messyBc_barcode_n)[i]
    # barcode number
    n <- summary_res$messyBc_barcode_n[i]
    stringr::str_glue("    In sample ${sample_name} there are: {n} Tags")
  }) %>% unlist %>% paste(collapse = "\n")

  # elements in metadata
  metadata_info <- paste(names(x$metadata), collapse = "  ")
  metadata_n <- length(names(x$metadata))

  
  res <- stringr::str_glue(
"Bonjour le monde, This is a BarcodeObj.
----------
It contains: 
    {subjects}
----------
$metadata: {metadata_n} meta data field(s) available
    {metadata_info}
----------
$messyBc: {messyBc_n} Samples for raw barcodes
{messyBc_info}")
  

  # cleanBc
  if (!is.null(summary_res$cleanBc_barcode_n)) {

    # sample number
    cleanBc_n <- length(summary_res$cleanBc_barcode_n)
    cleanBc_info <- lapply(seq_along(summary_res$cleanBc_barcode_n), 
      function(i) { 
        # sample name
        sample_name <- names(summary_res$cleanBc_barcode_n)[i]
        # sample number
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
#' # print BarcodeObj
#' print(bc_obj)
#' ###
#' @seealso [Bc::format.BarcodeObj] format BarcodeObj for pretty print
#' @export
print.BarcodeObj <- function(x, ...) {
  cat(format(x), "\n")
}
