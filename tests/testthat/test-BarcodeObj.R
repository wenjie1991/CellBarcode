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


d2 <- data.frame(
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

d_l <- list(test1 = d1, test2 = d2)

pattern <- "TCGATCGATCGA([ACTG]+)ATCGATCGATC"
bc_obj <- bc_extract(d_l, pattern, sample_name=c("test1", "test2"))
bc_obj <- bc_cure_depth(bc_obj, depth=6, isUpdate = FALSE)


test_that("data type transformation", {

  expect_equal(
    bc_2df(bc_obj[, "test1"]), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6), stringsAsFactors=FALSE))

  expect_equal(
    bc_2dt(bc_obj[, "test1"]), 
    data.table(sample_name = "test1", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6), stringsAsFactors=FALSE))

  out_test <- bc_2matrix(bc_obj) 
  out_test <- out_test[order(out_test[, 1]), ]
  out_truth <- matrix(c(104, 50, 14, 6, 104, 50, 14, 6), byrow=FALSE, ncol=2)
  dimnames(out_truth) <- list(c("AGAG", "AAAG", "AGAAG", "ACAAG"), c("test1", "test2"))
  out_truth <- out_truth[order(out_truth[, 1]), ]
  expect_equal(out_test, out_truth)
})


test_that("subset operation", {

  expect_equal(
    bc_2df(bc_obj[, "test1"]), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_obj[, "test1"]), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG", "AGAAG", "ACAAG"), count = c(104, 50, 14, 6), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_obj["AGAG", "test1"]), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG"), count = c(104), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_subset(bc_obj, barcode = "AGAG", sample = "test1")), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG"), count = c(104), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_obj["AGAG", "test1"] + bc_obj["AAAG", "test1"]), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG"), count = c(104, 50), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_cure_depth(bc_obj["AGAG", "test1"] + bc_obj["AGAG", "test2"], depth = 0)),
    data.frame(sample_name = c("test1", "test2"), barcode_seq = c("AGAG", "AGAG"), count = c(104, 104), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_obj[, "test1"] - "ACAAG")
    , 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG", "AGAAG"), count = c(104, 50, 14), stringsAsFactors=FALSE))

  expect_equal(
    bc_2df(bc_obj[, "test1"] * c("AGAG", "AAAG")), 
    data.frame(sample_name = "test1", barcode_seq = c("AGAG", "AAAG"), count = c(104, 50), stringsAsFactors=FALSE))

})

test_that("meta data", {

  expect_equal(bc_barcodes(bc_obj[, "test1"]), c("AGAG", "AAAG", "AGAAG", "ACAAG"))

  expect_equal(
    bc_names(bc_obj), c("test1", "test2"))

  
  bc_meta(bc_obj)$id <- c(1, 2) 
  expect_equal(bc_meta(bc_obj)$id, c(1, 2))

  bc_names(bc_obj) <- c("test11", "test12")
  expect_equal(
    bc_names(bc_obj), c("test11", "test12"))

  expect_snapshot(print(bc_obj))

})

## import processed data as a BarcodeObj
test_that("import", {
    x = bc_2matrix(bc_obj)
    expect_equal(
        bc_2matrix(
            bc_cure_depth(
                bc_create_BarcodeObj(x)
            )
        ), 
    x)

    x = data.frame(x)
    expect_equal(
        data.frame(
            bc_2matrix(
                bc_cure_depth(
                    bc_create_BarcodeObj(x)
                )
            )
        ),
    x)
})


