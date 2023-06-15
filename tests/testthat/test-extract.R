library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
sr <- list(readFastq(fq_file))
ds <- list(sr[[1]]@sread)
iv <- list(tables(ds[[1]], n = Inf)$top)
df <- list(data.frame(seq = names(iv[[1]]), freq = as.integer(iv[[1]])))
l <- list(sample1 = ds[[1]], sample2 = ds[[1]])

test_that("read in data", {

  expect_equal(
    bc_extract(sr, pattern = "AAAAA(.*)CCCCC")@messyBc[[1]]$barcode_seq, "GG")

  expect_equal(
    bc_extract(ds, pattern = "AAAAA(.*)CCCCC")@messyBc[[1]]$barcode_seq, "GG")

  expect_equal(
    bc_extract(iv, pattern = "AAAAA(.*)CCCCC")@messyBc[[1]]$barcode_seq, "GG")

  expect_equal(
    bc_extract(df, pattern = "AAAAA(.*)CCCCC")@messyBc[[1]]$barcode_seq , "GG")

  expect_equal(
    bc_extract(fq_file, pattern = "AAAAA(.*)CCCCC")@messyBc[[1]]$barcode_seq, "GG")

  expect_equal(
    bc_barcodes(bc_cure_depth(bc_extract(c(fq_file, fq_file), sample_name = c("test1", "test2"), pattern = "AAAAA(.*)CCCCC"), depth = 0)), "GG")

})





