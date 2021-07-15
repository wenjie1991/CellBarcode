## set up a dummy data set
d1 <- data.frame(
  seq = c(
    "ACTTCGATCGATCGAAAAGATCGATCGATC",
    "AATTCGATCGATCGAAGAGATCGATCGATC",
    "CCTTCGATCGATCGAAGAAGATCGATCGATC",
    "TTTTCGATCGATCGAAAAGATCGATCGATC",
    "AAATCGATCGATCGAAGAGATCGATCGATC",
    "CCCTCGATCGATCGAAGAAGATCGATCGATC",
    "GGGTCGATCGATCGAAAAGATCGATCGATC",
    "GGATCGATCGATCGAAGAGATCGATCGATC",
    "ACTTCGATCGATCGAACAAGATCGATCGATC",
    "GGTTCGATCGATCGACGAGATCGATCGATC",
    "GCGTCCATCGATCGAAGAAGATCGATCGATC"
    ),
  freq = c(
    30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
    )
  )


test_that("Senerio1: Backbone no error, Depth cutoff > 5", {
  pattern <- "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"))
  bc_obj <- bc_cure_depth(bc_obj, depth=5)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name = "test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6)))
})


test_that("Senerio1.1: Backbone no error, Depth cutoff > 5, hamming dist 1", {
  pattern <- "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"))
  bc_obj <- bc_cure_cluster(bc_cure_depth(bc_obj, depth=5), distance = 1)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name = "test", barcode_seq = c("AGAG", "AGAAG"), count = c(104 + 50, 14 + 6)))
})

test_that("Senerio1.2: Backbone no error, Depth cutoff > 5, levenshtein dist 1", {
  pattern <- "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"))
  bc_obj <- bc_cure_cluster(bc_cure_depth(bc_obj, depth=5), distance = 1, dist_method = "leven")
  expect_equal(bc_2df(bc_obj), data.frame(sample_name = "test", barcode_seq = c("AGAG"), count = c(104 + 50 + 14 + 6)))
})


test_that("Senerio2: Backbone 1 error, Depth cutoff > 5", {
  pattern <- "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern,  maxLDist=1)
  bc_obj <- bc_cure_depth(bc_obj, depth=5)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count=c(104, 50, 20, 6)))
})

test_that("Senerio3: Backbone no error, Depth cutoff > 5, No fishing, UMI is unique", {
  pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), pattern_type=c(UMI=1, barcode=2))
  bc_obj <- bc_cure_depth(bc_cure_umi(bc_obj, depth=5, doFish=F, isUniqueUMI=T), depth=0)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG"), count=c(3, 3, 1)))
})


test_that("Senerio4: Backbone 1 error, Depth cutoff > 5, Fishing, UMI is not unique", {
  pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), pattern_type=c(UMI=1, barcode=2), maxLDist=1)
  bc_obj <- bc_cure_depth(bc_cure_umi(bc_obj, depth=5, doFish=T, isUniqueUMI=F), depth=0)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count=c(3, 3, 3, 1)))
})

test_that("Senerio5: Backbone 1 error, Depth cutoff > 5, Fishing, UMI is not unique", {
  pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), pattern_type=c(UMI=1, barcode=2), maxLDist=1)
  bc_obj <- bc_cure_depth(bc_cure_umi(bc_obj, depth=5, doFish=T, isUniqueUMI=F), depth=0)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count=c(3, 3, 3, 1)))
})


