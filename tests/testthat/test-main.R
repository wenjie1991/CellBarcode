## set up a dummy data set
d = data.frame(
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

test_that("Senerio: Backbone no error, Depth cutoff > 5", {
  pattern = "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj = bc_extract(list(test = d), pattern, sample_name=c("test"))
  bc_obj = bc_cure(bc_obj, depth=5)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name = "test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6)))
})

test_that("Senerio2: Backbone 1 error, Depth cutoff > 5", {
  pattern = "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj = bc_extract(list(test = d), pattern, sample_name=c("test"), maxLDist=1)
  bc_obj = bc_cure(bc_obj, depth=5)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count=c(104, 50, 20, 6)))
})

test_that("Senerio3: Backbone no error, Depth cutoff > 5, No fishing, UMI is unique", {
  pattern = "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj = bc_extract(list(test = d), pattern, sample_name=c("test"), pattern_type=c(UMI=1, barcode=2))
  bc_obj = bc_cure(bc_obj, depth=5, doFish=F, with_umi=T, isUniqueUMI=T)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG"), count=c(3, 3, 1)))
})


test_that("Senerio4: Backbone 1 error, Depth cutoff > 5, Fishing, UMI is not unique", {
  pattern = "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
  bc_obj = bc_extract(list(test = d), pattern, sample_name=c("test"), pattern_type=c(UMI=1, barcode=2), maxLDist=1)
  bc_obj = bc_cure(bc_obj, depth=5, doFish=T, with_umi=T, isUniqueUMI=F)
  expect_equal(bc_2df(bc_obj), data.frame(sample_name="test", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count=c(3, 3, 3, 1)))
})
