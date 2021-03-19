#' @title Generate genBarCode object from Bc BarcodeObj
#'
#' @description Use the uncorrected barcode to generate the genBarCode data object
#' @param barcodeObj A BarcodeObj
#' @param BC_backbone A string, the backbone of the barcode
#' @param resDir A string, the location to store the data of genBarCode object
#' @return A list, contains genBarCode object
#' @export
bc_2genBarCode = function(barcodeObj, BC_backbone = "", resDir = getwd(), type = "messyBc") {
  count = barcode_seq = read_count = NULL # due to NOTE in check

  # TODO: include the genBaRcode barcode object, so that the user no need to load the genBaRcode package
  from_obj = barcodeObj[[type]]
  lapply(1:length(from_obj), function(i) {
    d = data.table::data.table(from_obj[[i]])[, .(read_count = sum(count)), by = barcode_seq][order(read_count, decreasing = T)][, .(read_count = read_count, barcode = factor(barcode_seq))]
    d = as.data.frame(d)
    return(methods::new(Class = "BCdat", reads = d, results_dir = resDir, label = names(from_obj)[i], BC_backbone = BC_backbone))
  })
}
