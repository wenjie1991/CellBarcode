# library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
sr <- readFastq(fq_file)
ds <- sr@sread
iv <- tables(ds, n = Inf)$top
df <- data.frame(seq = names(iv), freq = as.character(iv))
l <- list(sample1 = ds, sample2 = ds)
l_sr <- list(sample1 = sr, sample2 = sr)

bc_filterSeq(l_sr)

test_that("Sequence filter", {

  expect_equal(ShortRead::tables(bc_filterSeq(fq_file)[[1]])$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(sr))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(ds))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(iv))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(df))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(l)[[1]])$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_filterSeq(l_sr)[[1]])$top, c("AAAAAGGCCCCC" = 1))

})


