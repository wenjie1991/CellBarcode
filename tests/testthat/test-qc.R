# library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="Bc")
sr <- readFastq(fq_file)
ds <- sr@sread
iv <- tables(ds, n = Inf)$top
df <- data.frame(seq = names(iv), freq = as.character(iv))
l <- list(sample1 = ds, sample2 = ds)
l_sr <- list(sample1 = sr, sample2 = sr)

test_that("Quality Control", {

  #expect_s3_class(plot(bc_seqQC(l_sr)), "egg")
  #expect_s3_class(plot(bc_seqQC(fq_file[[1]])), "egg")

  expect_snapshot(bc_seqQC(sr))

  expect_snapshot(bc_seqQC(ds))

  expect_snapshot(bc_seqQC(iv))

  expect_snapshot(bc_seqQC(df))

  expect_snapshot(bc_seqQC(l))

  expect_snapshot(bc_seqQC(l_sr))
})

