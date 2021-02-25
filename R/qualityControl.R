get_base_quality_per_cycle = function(bstringset) {
  # Output: data.table Cycle, Mean, Median, P5, P95
  quality_m = bstringset
  quality_encode = encoding(quality_m)
  m = alphabetByCycle(quality_m)[-1, ]
  quality_v = quality_encode[rownames(m)]

  y = lapply(1:ncol(m), function(i) {
    S4Vectors::Rle(quality_v, m[, i])
  })

  quality_stat = sapply(y, function(x) {
    quantile_result = stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(quantile_result) = c("P5", "P25", "Median", "P75", "P95")
    c(Mean = mean(x), quantile_result)
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
#' 
#' fastq file name, ShortReadQ Obj, DNAStringSet Obj, data.frame, vector
#' 
runQC = function(...) UseMethod("runQC")

runQC.ShortReadQ = function(x, n = 50, plot = F) {
  # output: top, distribution (nOccurrences, nReads), base_quality_per_cycle, base_freq_per_cycle
  output = ShortRead::tables(x)
  output$base_quality_per_cycle = get_base_quality_per_cycle(x@quality)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x@sread)
  class(output) = append(class(output), "barcodeQc")
  output
}

runQC.DNAStringSet = function(x, plot = F) {
  # output: top, distribution (nOccurrences, nReads), base_freq_per_cycle
  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  class(output) = append(class(output), "barcodeQc")
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
  class(output) = append(class(output), "barcodeQc")
  output
}

runQC.integer = function(x) {
  x = DNAStringSet(rep(names(x), x))

  output = ShortRead::tables(x)
  output$base_freq_per_cycle = get_base_freq_per_cycle(x)
  class(output) = append(class(output), "barcodeQc")
  output
}

runQC.character = function(file) {
  # TODO: use qa to do fast check for the quality
  ## generate reads quality information
  # TODO: Create a new function to check the quality of fastq file.
  if (length(file) == 1) {
    runQC(ShortRead::readFastq(file))
  } else {
    qc_list = lapply(file, function(f) {
      runQC(ShortRead::readFastq(f))
    })
    names(qc_list) = basename(file)
    class(qc_list) = append(class(qc_list), "barcodeQcSet")
    qc_list
  }
}


runQC.list = function(x) {
  qc_list = lapply(x, runQC)
  class(qc_list) = append(class(output), "barcodeQcSet")
  qc_list
}

plot_reads_depth_distribution = function(distribution) {

  nReads = nOccurrences = NULL

  # distribution is data.frame with two columns "nOccurrences" and "nReads".
  d = distribution
  g = ggplot(d) + aes(x = nReads, y = nOccurrences) + geom_line() + theme_bw()
  g
}

plot_base_percentage_distribution = function(base_freq_per_cycle) {

  Cycle = Count = Base = NULL

  # base_freq_per_cycle is data.frame with three columns "Cycle", "Count" and "Base"
  d = data.table(base_freq_per_cycle)
  d[, all_base_per_cycle := sum(Count), by = .(Cycle)]
  d[, base_percent := Count / all_base_per_cycle, by = .(Cycle, Base)]
  g = ggplot(d) + aes(x = Cycle, y = base_percent, color = Base) + geom_point() + geom_line() + theme_bw()
  g
}

plot_base_quality_distribution = function(base_quality_per_cycle) {

  Cycle = Median = P95 = P75 = P25 = P5 = NULL

  # base_quality_per_cycle is data.table with columns Cycle, Mean, Median, P5, P95
  d = base_quality_per_cycle
  g = ggplot(d) + aes(x = Cycle, middle  = Median, upper = P75, lower = P25, ymax = P95, ymin = P5, group = Cycle) + geom_boxplot(stat = "identity") + theme_bw()
  g
}

plot.barcodeQc = function(x, opt = "all") {
  # barcodeQc has: top, distribution (nOccurrences, nReads), base_quality_per_cycle, base_freq_per_cycle
  # TODO: user can select which figure to draw
  g_list = list()
  if ("distribution" %in% names(x)) {
    g_list = append(g_list, list(plot_reads_depth_distribution(x$distribution)))
  }
  if ("base_freq_per_cycle" %in% names(x)) {
    g_list = base::append(g_list, list(plot_base_percentage_distribution(x$base_freq_per_cycle)))
  }
  if ("base_quality_per_cycle" %in% names(x)) {
    g_list = base::append(g_list, list(plot_base_quality_distribution(x$base_quality_per_cycle)))
  }
  g_num = length(g_list)
  col_n = ceiling(g_num / 2)
  egg::ggarrange(plots = g_list, nrow = 2, ncol=col_n)
}

plot.barcodeQcSet = function(x) {

  Count = Cycle = fileName = base_num = Base = base_percent = Median = NULL

  # TODO:
  # Raw data:
  #   Most frequency reads
  #   Reads length distribution
  #   nucleic acide per base
  #   quality per base
  g_list = list()

  d = lapply(1:length(x), function(i) { data.table(x[[i]]$base_freq_per_cycle, fileName = names(x)[i]) }) %>% rbindlist()
  d[, base_num := sum(Count), by = .(Cycle, fileName)]
  d[, base_percent := Count / base_num, by = .(Base, Cycle, fileName)]
  p1 = ggplot(d, aes(Cycle, fileName, fill = Base, alpha = base_percent)) + geom_tile(color = "black") + labs(y = "Sample Name") + theme_bw()
  g_list = append(g_list, list(p1))

  if ("base_quality_per_cycle" %in% names(x[[1]])) {
    d = lapply(1:length(x), function(i) { data.table(x[[i]]$base_quality_per_cycle, sample_name = names(x)[i]) }) %>% rbindlist()
    p2 = ggplot(d, aes(Cycle, sample_name, fill = Median)) + geom_tile(color = "black") + labs(y = "Sample Name", fill = "Median Base Quality") + theme_bw()
    g_list = append(g_list, list(p2))
  }

  row_n = length(g_list)
  egg::ggarrange(plots = g_list, nrow = row_n, ncol=1)
}
