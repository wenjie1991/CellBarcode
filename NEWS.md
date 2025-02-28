## Changes in Version 1.13.2

8a38560 [doc] Update package document
f37a920 [bug] bc_cure_cluster() error when there is no @cleanBc slot
dcee0d8 [doc] Update install guide
910ce93 [feat] transform the Rcpp to rust
d1d354f [misc] add error message when no barcode detected in bc_extract_sc_fastq()
63f7c95 [misc] update the error message when no barcode detected
f32bb33 [mis] Update action cache version
5c89289 [fix] zlibbioc Deprecation
4509346 [mis] Update Makevars.win
9589498 [mis] Document update
9f427dc [doc] Update license to Artistic-2.0
99af8f4 [doc] update scSeq vignettes

## Changes in Version 1.9.1

### New features

#### Adding function `bc_create_BarcodeObj` enabling import processed barcode data

a16a9e7 feat: Adding function bc_create_BarcodeObj enabling import process barcode data
45d886a test: adding test for processed data import
1d5f4de doc: add bc_create_BarcodeObj document

### Fix:

#### Make sure object is `data.frame` instead of `tibble` or `data.table`

895e9a9 fix: make sure the data.frame is pure data.frame instead of tibble or data.table
31e6989 fix: error while checking data.frame

### Chore:

#### Error messages and input checks

f84b2a1 chore: message if no barcodes found in bulk sequencings
e0cded9 chore: Error message when no barcodes are found
3beef87 chore: Check input for sam and bam file

## Changes in Version 1.7.2

### New features

- `bc_extract_sc_bam`: function to extract single-cell barcodes from BAM file.

## Changes in Version 1.7.0

### New features

- `bc_extract_sc_sam`: function to extract single-cell barcodes from SAM file
- `bc_extract_sc_fastq`: function to extract single-cell barcodes from fastq
  file
- `bc_plot_count`: function to plot the count feature of the barcodes such as
  barcode reads versus total reads, reads per UMI distribution.

### Incompatible changes

`bc_extract` now only accepts a vector or list, and returns `BarcodeObj`,
instead of returning `data.frame` when the input is a single sample.

## Changes in Version 1.4.0

Version bumping

## Changes in version 1.3.2

### Enhancements

- Reduce memory usage by removing the original sequence in barcodeObj.

- Impromve fastq reading efficiency by using Rcpp.

### New features

- `bc_splitVDJ` split VDJ barcode (experimental).

- `bc_extract_10XscSeq` function used to extract transcripted barcode in 10X
Genomics scRNASeq data.

### Bug fix

- fix the error of `bc_cure_cluster`.

## Changes in 1.2.0

Just version bumping

## Changes in 1.1.0

Just version bumping

## Start with 1.0.0
