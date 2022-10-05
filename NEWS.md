
## Changes in version 1.2.2

### Enhancements

- Reduce memory usage by removing the original sequence in barcodeObj.

- Impromve fastq reading efficiency by using Rcpp.

### New features

- `bc_splitVDJ` split VDJ barcode (experimental).

- `bc_extract_10X_scSeq` function used to extract transcripted barcode in 10X
Genomics scRNASeq data.

### Bug fix

- fix the error of `bc_cure_cluster`.
