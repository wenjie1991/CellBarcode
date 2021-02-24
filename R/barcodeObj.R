#' Subset operation of BarcodeObj object
#'
#' @param x A vector or string. Select the barcode subset from BarcodeObj
#' @param y A vector or string. Select the sample subset from BarcodeObj
#' @return A BarcodeObj
#' @export
#' @examples
#' # Not run
#' #bc_obj[c("AACCTT", "AACCTT"), ]
#' #bc_obj[, "sample1_rep1"]
#' #bc_obj[, c("sample1_rep1", "sample1_rep2")]
"[.BarcodeObj" = function(barcodeObj, x = NULL, y = NULL) {
  if (is.null(x) & is.null(y)) {
    error("Error message. TODO.")
  }
  if (!is.null(y)) {
    return(subset.BarcodeObj(barcodeObj, sample = y))
  }
  if (!is.null(x)) {
    return(subset.BarcodeObj(barcodeObj, barcode = x))
  }
  return(subset.BarcodeObj(barcodeObj, sample = y, barcode = x))
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
#' # subset(bc_obj, barcode = c("AACCTT", "AACCTT"))
#' # subset(bc_obj, sample = "sample1_rep1", barcode = c("AACCTT", "AACCTT"))
#' # subset(bc_obj, sample = c("sample1_rep1", "sample1_rep2"), barcode = c("AACCTT", "AACCTT"))
subset.BarcodeObj = function(barcodeObj, sample = NULL, barcode = NULL, black_list = NULL) {
  # TODO: The funciton only can apply the operation to the `messyBc` and `cleanBc`. We need to make it
  # capable to apply the selection to all information in the object.
  if (is.null(sample) & is.null(barcode)) {
    error("Error message. TODO.")
  }
  if (!is.null(barcode)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc = lapply(barcodeObj$messyBc, function(d) {
        d[reads_seq %in% barcode]
})
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc = lapply(barcodeObj$cleanBc, function(d) {
        d[reads_seq %in% barcode]
})
    }
  }
  if (!is.null(black_list)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc = lapply(barcodeObj$messyBc, function(d) {
        d[!(reads_seq %in% black_list)]
})
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc = lapply(barcodeObj$cleanBc, function(d) {
        d[!(reads_seq %in% black_list)]
})
    }
  }
  if (!is.null(sample)) {
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$messyBc =  barcodeObj$messyBc[sample]
    }
    if (!is.null(barcodeObj$cleanBc)) {
      barcodeObj$cleanBc =  barcodeObj$cleanBc[sample]
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
  x$messyBc = append(x$messyBc, y$messyBc)
  x$cleanBc = append(x$cleanBc, y$cleanBc)
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


sampleNames = function(...) UseMethod("rename")

sampleNames.BarcodeObj = function(barcodeObj, rename) {
  if (!is.null(barcodeObj$messyBc)) {
    old_name = names(barcodeObj$messyBc)
    old_name[names(rename)] = rename
    names(barcodeObj$messyBc) = old_name
  }

  if (!is.null(barcodeObj$cleanBc)) {
    old_name = names(barcodeObj$cleanBc)
    old_name[names(rename)] = rename
    names(barcodeObj$cleanBc) = old_name
  }
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
  paste(cat(names(barcodeObj), collapse = "  "), "\n----------\n")
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
    cat("$cleanBc: ", length(barcodeObj$cleanBc), "Samples for cleaned barcodes\n")
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
