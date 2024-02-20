# CellBarcode

[![check](https://github.com/wenjie1991/CellBarcode/actions/workflows/r.yml/badge.svg)](https://github.com/wenjie1991/CellBarcode/actions)
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/CellBarcode.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/CellBarcode)

**CellBarcode** is an R package for dealing with **Cellular DNA barcoding** sequencing data.

The R package was created by Wenjie SUN, Anne-Marie Lyne, and Leïla Perié at Institut Curie.

## Types of barcodes

**CellBarcode** can handle all types of DNA barcodes, provided that:

- The barcodes have a pattern that can be matched by a regular expression.
- Each barcode is within a single sequencing read.

## What you can do with **CellBarcode**

- Perform quality control for the DNA sequence results, and filter the sequences according
  to their quality metrics.

- Identify barcode (and UMI) information in sequencing results.

- Performs quality control and deal with the spurious sequences that come from potential PCR & sequence errors.

- Provide toolkits to make it easier to manage samples and barcodes with metadata.

## Installing

### Install the development version from GitHub

```r
if(!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("wenjie1991/CellBarcode")
```

## Getting Started

Here is an example of a basic workflow:

```r
library(CellBarcode)
library(magrittr)

# The example data is a mix of MEF lines with known barcodes
# 2000 reads for each file have been sampled for this test dataset
# Data can be accessed here: https://zenodo.org/records/10027002
example_data <- system.file("extdata", "mef_test_data", package = "CellBarcode")
fq_files <- dir(example_data, "gz", full=TRUE)

# prepare metadata
metadata <- stringr::str_split_fixed(basename(fq_files), "_", 10)[, c(4, 6)]
metadata <- data.frame(metadata)
sample_name <- apply(metadata, 1, paste, collapse = "_")
colnames(metadata) = c("cell_number", "replication")
rownames(metadata) = sample_name
metadata

# extract UMI barcode with regular expression
bc_obj <- bc_extract(
  fq_files,
  pattern = "(.{12})CTCGAGGTCATCGAAGTATCAAG(.+)TAGCAAGCTCGAGAGTAGACCTACT", 
  pattern_type = c("UMI" = 1, "barcode" = 2),
  sample_name = sample_name,
  metadata = metadata
)
bc_obj

# sample subset operation, select 'mixa'
bc_sub <- bc_subset(bc_obj, sample=replication == "mixa")
bc_sub 

# filter the barcode, UMI barcode amplicon >= 2 & UMI counts >= 2
bc_sub <- bc_cure_umi(bc_sub, depth = 2) %>% bc_cure_depth(depth = 2)

# select barcodes with a white list
bc_sub[c("AAGTCCAGTACTATCGTACTA", "AAGTCCAGTACTGTAGCTACTA"), ]

# export the barcode counts to data.frame
head(bc_2df(bc_sub))

# export the barcode counts to matrix
head(bc_2matrix(bc_sub))
```

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Citation

If you use **CellBarcode** in your research, please cite the following paper:
[Sun, W. et al. Extracting, filtering and simulating cellular barcodes using CellBarcode tools. Nat Comput Sci 1–16 (2024)](https://www.nature.com/articles/s43588-024-00595-7)

