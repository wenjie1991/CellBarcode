# CellBarcode

**CellBarcode** is an R package for dealing with **Cellular DNA barcoding** sequencing data.

## Kinds of barcodes

**CellBarcode** handle all kinds of DNA barcodes, as long as:

- The barcode have a pattern which be matched by a regular expression.
- Each barcode is within a single sequencing read.

## What you can do with **CellBarcode**

- Performs quality control the DNA sequence results and filters the sequences according
  to the quality metrics.

- Parses sequences, extracts barcode (and UMI) information.

- Performs quality control and filters the barcode.

- Provides toolkits make it easier to manipulate samples and barcodes with metadata.

## Installing

1. Use the `devtools` install package from GitHub

```
library(devtools)
install_github("wenjie1991/CellBarcode")
```

2. Bioconductor
TBD (I hope the package can be accepted in Bioconductor).

## Get start

Here is an example of a basic workflow:

```{r}
library(CellBarcode)
library(magrittr)

# The example data is the mix of MEF lines with known barcodes
# 2000 reads for each file have been sampled for this test dataset
# TODO: Citation of the paper:
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

# sample subset operation, select technical repeats 'mixa'
bc_sub = bc_obj[, replication == "mixa"]
bc_sub 

# filter the barcode, UMI barcode amplicon > 1 & UMI counts > 1
bc_sub <- bc_cure_umi(bc_sub, depth = 1) %>% bc_cure_depth(depth = 1)

# select barcodes with a white list
bc_sub[c("AAGTCCAGTACTATCGTACTA", "AAGTCCAGTACTGTAGCTACTA"), ]

# export the barcode counts to data.frame
head(bc_2df(bc_sub))

# export the barcode counts to matrix
head(bc_2matrix(bc_sub))
```

## Depends

R (>= 4.0.0)
methods,
Rcpp (>= 1.0.5),
data.table (>= 1.14.0),
plyr (>= 1.8.6),
ggplot2 (>= 3.3.5),
stringr (>= 1.4.0),
magrittr (>= 2.0.1),
ShortRead (>= 1.48.0),
Biostrings (>= 2.58.0),
egg (>= 0.4.5),
utils,
S4Vectors

## License

[MIT](https://choosealicense.com/licenses/mit/)
