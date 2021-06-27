library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="Bc")
sr <- readFastq(fq_file)
ds <- sr@sread
iv <- tables(ds, n = Inf)$top
df <- data.frame(seq = names(iv), freq = as.character(iv))
l <- list(sample1 = ds, sample2 = ds)

test_that("read in data", {

  expect_equal(
    bc_extract(sr, pattern = "AAAAA(.*)CCCCC")$barcode_seq, "GG")

  expect_equal(
    bc_extract(ds, pattern = "AAAAA(.*)CCCCC")$barcode_seq, "GG")

  expect_equal(
    bc_extract(iv, pattern = "AAAAA(.*)CCCCC")$barcode_seq, "GG")

  expect_equal(
    bc_extract(df, pattern = "AAAAA(.*)CCCCC")$barcode_seq, "GG")

  expect_equal(
    bc_extract(fq_file, pattern = "AAAAA(.*)CCCCC")$barcode_seq, "GG")

  expect_equal(
    bc_barcodes(bc_cure_depth(bc_extract(c(fq_file, fq_file), sample_name = c("test1", "test2"), pattern = "AAAAA(.*)CCCCC"))),
    "GG")

})





