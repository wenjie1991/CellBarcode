use extendr_api::prelude::*;

mod lib_10x_barcode;
mod lib_clustering;
mod lib_read_seq;

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod CellBarcode;
    use lib_10x_barcode;
    use lib_clustering;
    use lib_read_seq;
}

