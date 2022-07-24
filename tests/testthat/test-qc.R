# library(ShortRead)

fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
sr <- readFastq(fq_file)
ds <- sr@sread
iv <- tables(ds, n = Inf)$top
df <- data.frame(seq = names(iv), freq = iv)

l_ds <- list(sample1 = ds, sample2 = ds)
l_sr <- list(sample1 = sr, sample2 = sr)

test_that("Quality Control", {

  #expect_s3_class(plot(bc_seq_qc(l_sr)), "egg")
  #expect_s3_class(plot(bc_seq_qc(fq_file[[1]])), "egg")

  expect_snapshot(
      bc_seq_qc(sr)
  )

  expect_snapshot(
      bc_seq_qc(ds)
  )

  expect_snapshot(
      bc_seq_qc(iv)
  )

  expect_snapshot(
      bc_seq_qc(df)
  )

  expect_snapshot(
      bc_seq_qc(l_ds)
  )

  expect_snapshot(
      bc_seq_qc(l_sr)
  )

})

