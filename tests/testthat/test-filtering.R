# library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
sr <- readFastq(fq_file)
ds <- sr@sread
iv <- tables(ds, n = Inf)$top
df <- data.frame(seq = names(iv), freq = as.character(iv))
l <- list(sample1 = ds, sample2 = ds)
l_sr <- list(sample1 = sr, sample2 = sr)

bc_seq_filter(l_sr)

test_that("Sequence filter", {

  expect_equal(ShortRead::tables(bc_seq_filter(fq_file)[[1]])$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(sr))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(ds))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(iv))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(df))$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(l)[[1]])$top, c("AAAAAGGCCCCC" = 1))

  expect_equal(ShortRead::tables(bc_seq_filter(l_sr)[[1]])$top, c("AAAAAGGCCCCC" = 1))

})


