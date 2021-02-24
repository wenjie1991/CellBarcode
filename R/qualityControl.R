get_base_quality_per_cycle = function(bstringset) {
  # Output: data.table Cycle, Mean, Median, P5, P95
  quality_m = bstringset
  quality_encode = encoding(quality_m)
  m = alphabetByCycle(quality_m)[-1, ]
  quality_v = quality_encode[rownames(m)]

  y = lapply(1:ncol(m), function(i) {
    Rle(quality_v, m[, i])
  })

  quality_stat = sapply(y, function(x) {
    c(Mean = mean(x), Median = stats::median(x), P5 = stats::quantile(x, 0.5)[[1]], P95 = stats::quantile(x, 0.95)[[1]])
  }) 
  quality_stat = t(quality_stat)
  data.table(Cycle = seq(1, nrow(quality_stat)), quality_stat)
}

get_base_freq_per_cycle = function(dnastringset) {
  # Output: data.table Cycle, Count, Base
  base_m = alphabetByCycle(dnastringset)
  base_m = t(base_m[c("A", "C", "G", "T"), ])
  base_m = data.table(Cycle = seq(1:nrow(base_m)), base_m)
  melt(base_m, id.vars = "Cycle", measure.vars = c("A", "C", "G", "T"), value.name = "Count", variable.name = "Base")
}

#' get the QC information of the raw sequence
#' fastq file name, ShortReadQ Obj, DNAStringSet Obj, data.frame, vector
runQC = function(...) UseMethod("runQC")

runQC.ShortReadQ = function(x, n = 50, plot = F) {
  # output: top, distribution (nOccurences, nReads), base_quality_per_cycle, base_freq_per_cycle
  output = ShortRead::tables(x)
  output$base_quality_per_cycle = get_base_quality_per_cycle(x@quality)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x@sread)
  output
}

runQC.DNAStringSet = function(x, plot = F) {
  # output: top, distribution (nOccurences, nReads), base_freq_per_cycle
  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  output
}

runQC.data.frame = function(x) {
  ## Input: data.frame(seq, freq)
  ## convert data.frame to DNAStringSet
  sequences = x$seq
  freq = x$freq
  x = DNAStringSet(rep(sequences, freq))

  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  output
}

runQC.character = function(x) {
  x = DNAStringSet(x)

  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  output
}

runQC.integer = function(x) {
  x = DNAStringSet(rep(names(x), x))

  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  output
}

runQC.list = function(x) {
  lapply(x, runQC)
}

runFqQC = function(files) {
  # TODO: use qa to do fast check for the quality
  ## generate reads quality information
  # TODO: Create a new function to check the quality of fastq file.
  lapply(files, function(file) {
    runQC(ShortRead::readFastq(file))
  })
}

