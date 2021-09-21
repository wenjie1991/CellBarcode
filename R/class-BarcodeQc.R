## barcodeQc class
###########################

#' @rdname bc_seq_qc
#' @exportClass BarcodeQc
BarcodeQc <- setClass("BarcodeQc",
    slots=list(
        top="integer",
        distribution="data.frame",
        base_quality_per_cycle="ANY",
        base_freq_per_cycle="data.frame",
        summary="numeric"
        ),
    )

## barcodeQcSet class
###########################
###

#' @rdname bc_seq_qc
#' @exportClass BarcodeQcSet
BarcodeQcSet <- setClass("BarcodeQcSet",
    slots=list(
        qc_list="list"
        )
    )


