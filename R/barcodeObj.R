check_sample_name = function(barcodeObj) {
  meta_name = rownames(barcodeObj$meta_data)
  if (all(meta_name != names(barcodeObj$messyBc))) stop("The messyBc sample names are not consistent with metadta")
  if (!is.null(barcodeObj$cleanBc)) {
    if (all(meta_name != names(barcodeObj$cleanBc))) stop("The cleanBc sample names are not consistent with metadta")
  }
}

#' Subset operation of BarcodeObj object
#'
#' @param barcodeObj A BarcodeObj
#' @param x A vector of integer or string, select the barcode subset from BarcodeObj
#' @param y A vector or an expression, select the sample subset from BarcodeObj
#' @return A BarcodeObj
#'
#' @export
#' @examples
#' data(bc_obj)
#' bc_obj[c("AGAG", "AAAG"), ]
#' bc_obj[, "test1"]
#' bc_obj[, sample_name == "test1"]
#' bc_obj[, c("test1", "test2")]
#'
"[.BarcodeObj" = function(barcodeObj, x = NULL, y = NULL) {

  y_call = substitute(y)
  return(subset.BarcodeObj(barcodeObj, sample = y_call, barcode = x, sample_is_expression = T))
}

#' Subset operation of BarcodeObj object
#'
#' @param barcodeObj A BarcodeObj
#' @param sample A vector of integer or string, select the barcode subset from BarcodeObj
#' @param barcode A vector or an expression, select the sample subset from BarcodeObj
#' @param black_list A vector or string. Remove the the barcode in the list.
#' @return A BarcodeObj
#'
#' @export
#' @examples
#' # Not run
#' # subset(bc_obj, barcode = c("AACCTT", "AACCTT"))
#' # subset(bc_obj, sample = "sample1_rep1", barcode = c("AACCTT", "AACCTT"))
#' # subset(bc_obj, sample = c("sample1_rep1", "sample1_rep2"), barcode = c("AACCTT", "AACCTT"))
subset.BarcodeObj = function(barcodeObj, sample = NULL, barcode = NULL, black_list = NULL, sample_is_expression = F) {

  reads_eq = NULL # due to NOTE in check

  meta_data = barcodeObj$meta_data
  if (sample_is_expression) {
    sample_call = sample
  } else {
    sample_call = substitute(sample)
  }

  # TODO: The funciton only can apply the operation to the `messyBc` and `cleanBc`. We need to make it
  # capable to apply the selection to all information in the object.

  ## TODO: How do handle messyBc
  if (!is.null(barcode)) {
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc = lapply(barcodeObj$messyBc, function(d) { d[barcode_seq %in% barcode] })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc = lapply(barcodeObj$cleanBc, function(d) { d[barcode_seq %in% barcode] })
    }
  }
  if (!is.null(black_list)) {
    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc = lapply(barcodeObj$messyBc, function(d) { d[!(barcode_seq %in% black_list)] })
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc = lapply(barcodeObj$cleanBc, function(d) { d[!(barcode_seq %in% black_list)] })
    }
  }
  if (!is.null(sample_call)) {
    check_sample_name(barcodeObj)

    sample_i = eval(sample_call, meta_data, parent.frame())

    barcodeObj$meta_data = meta_data[sample_i, , drop=F]

    if (!is.null(barcodeObj$messyBc)) {
      barcodeObj$messyBc =  barcodeObj$messyBc[rownames(barcodeObj$meta_data)]
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc =  barcodeObj$cleanBc[rownames(barcodeObj$meta_data)]
    }
  }
  return(barcodeObj)
}

#' Combined the BarcodeObj
#'
#' @param x A BarcodeObj
#' @param y A BarcodeObj
#' @return A BarcodeObj
#'
#' @export
"+.BarcodeObj" = function(x, y) {
  # TODO: Apply the merge to all parts of the data
  #       How to deal when two BarcodeObj have the same samples, the same samples will merged

  meta_data_x = x$meta_data
  meta_data_y = y$meta_data

  meta_data_xy = merge(meta_data_x, meta_data_y, by = "sample_name", all = T)
  rownames(meta_data_xy) = meta_data_xy$sample_name

  x$messyBc = lapply(rownames(meta_data_xy), function(sample_name) {
    rbind(x$messyBc[[sample_name]], y$messyBc[[sample_name]])
  })
  names(x$messyBc) = meta_data_xy$sample_name

  if (!is.null(x$cleanBc) & !is.null(y$cleanBc)) {

    x$cleanBc = lapply(rownames(meta_data_xy), function(sample_name) {
      rbind(x$cleanBc[[sample_name]], y$cleanBc[[sample_name]])
    })
    names(x$cleanBc) = meta_data_xy$sample_name
  }
  x$meta_data = meta_data_xy
  x
}

#' Remove the barcode in black list
#'
#' @param x A BarcodeObj
#' @param y A string or vector, the barcode black list to revmoe
#' @return A BarcodeObj
#'
#' @export
"-.BarcodeObj" = function(x, y) {
  subset(x, black_list = y)
}

#' Select the barcode by white list
#'
#' @param x A BarcodeObj
#' @param y A string or vector, the barcode white list
#' @return A BarcodeObj
#'
#' @export
"*.BarcodeObj" = function(x, y) {
  subset(x, barcode = y)
}

#' Get barcode sequences
#'
#' @export
bc_barcodes = function(barcodeObj, unlist = T) {
  d = lapply(barcodeObj$cleanBc, function(x) { x$barcode_seq }) 
  if (unlist) {
    unlist(d) %>% as.character %>% unique
  } else {
    names(d) = names(barcodeObj$cleanBc)
    d
  }
}

#' Get sample names
#'
#' @export
bc_names = function(barcodeObj) {
  check_sample_name(barcodeObj)
  rownames(barcodeObj$meta_data)
}


#' Assign sample names
#'
#' @export
"bc_names<-" = function(barcodeObj, value) {
  check_sample_name(barcodeObj)
  if (length(value) != nrow(barcodeObj$meta_data)) stop("The sample names do not have correct length.")
  if (!is.null(barcodeObj$messyBc)) {
    names(barcodeObj$messyBc) = value
  } 
  if (!is.null(barcodeObj$cleanBc)) {
    names(barcodeObj$cleanBc) = value
  }
  rownames(barcodeObj$meta_data) = value
  barcodeObj$meta_data$sample_name = value
  barcodeObj
}

#' Get meta data
#'
#' @export
bc_meta = function(barcodeObj) {
  check_sample_name(barcodeObj)
  barcodeObj$meta_data
}


#' Assign meta data
#'
#' @export
"bc_meta<-" = function(barcodeObj, key = NULL, value) {
  check_sample_name(barcodeObj)
  if (is.null(key)) {
    if (nrow(value) != nrow(barcodeObj$meta_data)) stop("The given meta data does not hvae correct length.")
    barcodeObj$meta_data = value
  } else {
    if (length(value) != nrow(barcodeObj$meta_data)) stop("The given meta data does not hvae correct length.")
    barcodeObj$meta_data[[key]] = value
  }
  barcodeObj
}


#' Get barcode reslt data.table format
#'
#' @export
bc_2dt = function(barcodeObj) {
  # TODO: enable messyBc
  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found.")
  }
  d = barcodeObj$cleanBc %>% rbindlist(idcol = T)
  names(d)[1] = "sample_name"
  return(d)
}

#' Get barcode reslt data.frame format
#'
#' @export
bc_2df = function(barcodeObj) {
  bc_2dt(barcodeObj) %>% as.data.frame
}

#' Get barcode result in matrix format
#'
#' @export
bc_2matrix = function(barcodeObj) {
  # TODO: enable messyBc
  if (is.null(barcodeObj$cleanBc)) {
    stop("No cleanBc found.")
  }
  # TODO: use sparse matrix
  d = bc_2dt(barcodeObj) %>% data.table::dcast(barcode_seq ~ sample_name, value.var = "count", fill = 0) 
  m = data.frame(d[, -1]) %>% data.matrix
  rownames(m) = d[, barcode_seq]
  m
} 

# TODO: How to update the sample name?
# sampleNames = function(...) UseMethod("rename")
# 
# sampleNames.BarcodeObj = function(barcodeObj, rename) {
#   if (!is.null(barcodeObj$messyBc)) {
#     old_name = names(barcodeObj$messyBc)
#     old_name[names(rename)] = rename
#     names(barcodeObj$messyBc) = old_name
#   }
# 
#   if (!is.null(barcodeObj$cleanBc)) {
#     old_name = names(barcodeObj$cleanBc)
#     old_name[names(rename)] = rename
#     names(barcodeObj$cleanBc) = old_name
#   }
#   barcodeObj
# }


summary_BarcodeObj = function(barcodeObj) {

  summary_res = list()
  summary_res$messyBc_barcode_n = lapply(barcodeObj$messyBc, nrow)
  names(summary_res$messyBc_barcode_n) = rownames(barcodeObj$meta_data)
  if (!is.null(barcodeObj$cleanBc)) {
    summary_res$cleanBc_barcode_n = lapply(barcodeObj$cleanBc, nrow)
    names(summary_res$cleanBc_barcode_n) = rownames(barcodeObj$meta_data)
  }

  summary_res
}

#' Format the summary of BarcodeObj
#'
#' @export
format.BarcodeObj = function(barcodeObj) {

  summary_res = summary_BarcodeObj(barcodeObj)

  subjects = paste(names(barcodeObj), collapse = "  ")
  messyBc_n = length(summary_res$messyBc_barcode_n)
  messyBc_info = lapply(seq_along(summary_res$messyBc_barcode_n), function(i) { 
    sample_name = names(summary_res$messyBc_barcode_n)[i]
    n = summary_res$messyBc_barcode_n[i]
    stringr::str_glue("    In sample ${sample_name} there are: {n} Tags")
  }) %>% unlist %>% paste(collapse = "\n")

  meta_data_info = paste(names(barcodeObj$meta_data), collapse = "  ")
  meta_data_n = length(names(barcodeObj$meta_data))

  
  res = stringr::str_glue(
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

    cleanBc_n = length(summary_res$cleanBc_barcode_n)
    cleanBc_info = lapply(seq_along(summary_res$cleanBc_barcode_n), function(i) { 
      sample_name = names(summary_res$cleanBc_barcode_n)[i]
      n = summary_res$cleanBc_barcode_n[i]
      stringr::str_glue("    In sample ${sample_name} there are: {n} barcodes")
    }) %>% unlist %>% paste(collapse = "\n")
    res_cleanBc = stringr::str_glue(
"----------
$cleanBc: {cleanBc_n} Samples for cleaned barcodes
{cleanBc_info}")
    res = stringr::str_c(res, "\n", res_cleanBc)
  }
  res
}

#' Print out the summary of BarcodeObj
#'
#' @param barcodeObj A BarcodeObj
#' @return The general information about the BarcodeObj
#'
#' @export
print.BarcodeObj = function(barcodeObj) {
  cat(format(barcodeObj), "\n")
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

# view = function(x, ...) UseMethod("view", x)
# view.BarcodeObj = function(barcodeObj) {
#   print.default(barcodeObj)
# }
