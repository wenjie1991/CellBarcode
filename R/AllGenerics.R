#' Manages barcodes and samples in a BarcodeObj object
#'
#' A set of functions and operators for subset or join of
#' BarcodeObj object(s). 
#' The \code{bc_subset}, \code{*} and \code{-} are used to select barcodes or
#' samples in a \code{BarcodeObj} object.
#' Two BarcodeObj objects can be joined by \code{+}.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param barcode A vector of integer or string, indicating the selected
#' barcode.
#' @param sample A character vector or integer vector or an expression
#' (expressio not applicable for \code{[]} operator), specifying the samples in
#' the subsets. When the value is an expression, the columns in the metadata can
#' be used as variable.
#' @param black_list A character vector, specifying the black list with excluded
#' barcodes.
#' @param white_list A character vector, giving the barcode white list. 
#' @param barcodeObj_x A BarcodeObj object.
#' @param barcodeObj_y A BarcodeObj object.
#' @param is_sample_quoted_exp A logical value. If TRUE, the expression in
#' \code{sample} argument will not be evaluated before executing the function. 
#' @return A BarcodeObj object.
#'
#' @details
#' \code{bc_subset} and \code{[]}: Gets samples and barcodes subset from a
#' \code{BarcodeObj} object.
#'
#' \code{+}: Combines two \code{BarcodeObj} objects. The \code{metadata},
#' \code{cleanBc} and
#' \code{messyBc} slot in the BarcodeObj objects will be joined. 
#' For the \code{metadata} slot, the \code{sample_name} column, and the
#' \emph{Full outer join} (the record in either BarcodeObj object) will be
#' performed with rownames as the key.
#' The \code{messyBc} and \code{cleanBc} from two objects are combined by rows
#' for the same sample from two \code{BarcodeObj} objects.
#' 
#' \code{-}: removes barcodes in the black_list.
#'
#' \code{*}: selects barcodes in the white_list.
#'
#' @examples
#' data(bc_obj)
#'
#' bc_obj
#'
#' # Select barcodes
#' bc_subset(bc_obj, barcode = c("AACCTT", "AACCTT"))
#' bc_obj[c("AGAG", "AAAG"), ]
#'
#' # Select samples by meta data
#' bc_meta(bc_obj)$phenotype <- c("l", "b")
#' bc_meta(bc_obj)
#' bc_subset(bc_obj, phenotype == "l")
#'
#' # Select samples by sample name
#' bc_obj[, "test1"]
#' bc_obj[, c("test1", "test2")]
#' bc_subset(bc_obj, sample = "test1", barcode = c("AACCTT", "AACCTT"))
#'
#' # Apply barcodes black list
#' bc_subset(
#' bc_obj,
#'     sample = c("test1", "test2"),
#'     barcode = c("AACCTT"))
#'
#' # Join two samples with different barcode sets
#' bc_obj["AGAG", "test1"] + bc_obj["AAAG", "test2"]
#'
#' # Join two samples with overlap barcodes
#' bc_obj_join <- bc_obj["AGAG", "test1"] + bc_obj["AGAG", "test2"]
#' bc_obj_join
#' # The same barcode will removed after applying bc_cure_depth()
#' bc_cure_depth(bc_obj_join)
#'
#' # Remove barcodes
#' bc_obj
#' bc_obj - "AAAG"
#'
#' # Select barcodes in white list
#' bc_obj
#' bc_obj * "AAAG"
#' ###
#' @rdname bc_subset
#' @export
setGeneric("bc_subset", 
    function(
        barcodeObj, 
        sample=NULL, 
        barcode=NULL, 
        black_list=NULL,
        is_sample_quoted_exp=FALSE) {
        standardGeneric("bc_subset") 
    }
)

# setGeneric("[", function(barcodeObj, barcode, sample) { standardGeneric("[") 
# })
# setGeneric("+", function(barcodeObj_x, barcodeObj_y) { standardGeneric("+") 
# })
# setGeneric("-", function(barcodeObj, black_list) { standardGeneric("-") })
# setGeneric("*", function(barcodeObj, white_list) { standardGeneric("*") })

#' @rdname bc_subset
#' @export
setGeneric("bc_merge", function(barcodeObj_x, barcodeObj_y) {
    standardGeneric("bc_merge") })

#' Gets barcode sequences
#'
#' \code{bc_barcodes} used to get the barcode sequences in \code{BarcodeObj}
#' object. The input 
#' \code{BarcodesObj} object should be pre-processed by \code{bc_cure_*}
#' functions, such as \code{bc_cure_depth}, \code{bc_cure_umi}.
#'
#' @param barcodeObj A \code{BarcodeObj} object.
#' @param unlist A logical value. If TRUE, the function returns a vector of
#' unique barcode list from all samples; otherwise a list will be returned. In
#' the later case, each element of the list contains the barcodes of a sample.
#' @return A character vector or a list.
#' @examples
#' data(bc_obj)
#'
#' # get unique barcode vector of all samples
#' bc_barcodes(bc_obj)
#'
#' # get a list with each element containing barcodes from one sample
#' bc_barcodes(bc_obj, unlist = FALSE)
#'
#' ###
#' @rdname bc_barcodes
#' @export
setGeneric("bc_barcodes", function(barcodeObj, unlist=TRUE) {
    standardGeneric("bc_barcodes") })

#' Access & update sample names in BarcodeObj & and BarcodeQcSet
#'
#' Get or update sample names in BarcodeObj object and BarcodeQcSet.
#'
#' @param x A \code{BarcodeObj} object or a \code{BarcodeQcSet} object.
#' @param value A character vector setting the new sample names, with the length
#' of the samples number in \code{BarcodeObj} or \code{BarcodeQcSet} object.
#' @return A character vector
#' @examples
#' data(bc_obj)
#'
#' bc_names(bc_obj)
#' bc_names(bc_obj) <- c("new1", "new2")
#' @rdname bc_names
#' @export
setGeneric("bc_names", function(x) { standardGeneric("bc_names") })


#' @rdname bc_names
#' @export
setGeneric("bc_names<-", function(x, value) { standardGeneric("bc_names<-") })

#' Accesses messyBc slot in the BarcodeObj object
#'
#' \code{messyBc} slot of BarcodeObj object contains the raw barcode reads
#' frequency data. For more detail about the \code{messyBc} slot, see
#' \code{\link[CellBarcode]{BarcodeObj}}. \code{bc_messyBc} is used to access
#' the `messyBc` slot in the \code{BarcodeObj}.
#'
#' @param barcodeObj A \code{BarcodeObj} objects.
#' @param isList A logical value, if TRUE (default), the return is a list with each
#' sample as an element. Otherwise, the function will return a \code{data.frame}
#' contains the data from all the samples with a column named \code{sample_name}
#' to keep the sample information.
#' @return
#' If a \code{list} is requested, in the \code{list} each element is a
#'  \code{data.frame} corresponding to the successive samples. Each
#'  \code{data.frame} has at most 3 columns: 1. \code{umi_seq} (optional): UMI
#' sequence. 2. \code{barcode_seq}: barcode sequence. 3. \code{count}: how many
#' reads a full sequence has. 
#'
#' If a \code{data.frame} is requested, the \code{data.frame} in the list
#' described above are combined into one \code{data.frame} by row, with an extra
#' column named \code{sample_name} for identifying sample.
#'
#' @examples
#'
#'  data(bc_obj)
#' # get the data in messyBc slot
#' # default the return value is a list
#' bc_messyBc(bc_obj)
#'
#' # the return value can be a data.frame
#' bc_messyBc(bc_obj, isList=FALSE)
#' ###
setGeneric("bc_messyBc", function(barcodeObj, isList=TRUE){ standardGeneric("bc_messyBc") })

#' Accesses cleanBc slot in the BarcodeObj object
#'
#' \code{cleanBc} slot of BarcodeObj object contains the processed barcode reads
#' frequency data. For more detail about the \code{cleanBc} slot, see
#' \code{\link[CellBarcode]{BarcodeObj}}. \code{bc_cleanBc} is used to access
#' the `cleanBc` slot in the \code{BarcodeObj}.
#'
#' @param barcodeObj A \code{BarcodeObj} objects.
#' @param isList A logical value, if TRUE (default), the return is a list with each sample
#' as an element. Otherwise, the function will return a \code{data.frame}
#' contains the data from all the samples with a column named \code{sample_name}
#' to keep the sample information.
#' @return
#' If a \code{list} is requested, each \code{list} element is a \code{data.frame}
#' for each sample. In a \code{data.frame}, there are 2 columns 1.
#' \code{barcode_seq}: barcode sequence 2. \code{counts}: reads count, or UMI
#' count if the \code{cleanBc} was created by \code{bc_cure_umi}.
#'
#' If a \code{data.frame} is requested, the \code{data.frame} in the list
#' described above are combined into one \code{data.frame} by row, with an extra
#' column named \code{sample_name} for identifying sample.
#'
#' @examples
#'
#' data(bc_obj)
#' # get the data in cleanBc slot
#' # default the return value is a list
#' bc_cleanBc(bc_obj)
#'
#' # the return value can be a data.frame
#' bc_cleanBc(bc_obj, isList=FALSE)
#' ###
setGeneric("bc_cleanBc", function(barcodeObj, isList=TRUE){ standardGeneric("bc_cleanBc") })


#' Accesses and sets metadata in BarcodeObj object
#'
#' Sample information is kept in metadata. \code{bc_meta} is for accessing and
#' updating metadata in \code{BarcodeObj} object
#'
#' @param barcodeObj A \code{BarcodeObj} object.
#' @param key A string, identifying the metadata record name to be modified.
#' @param value A string vector or a data.frame. If the \code{value} is a
#' vector, it should have the same length of sample number in the BarcodeObj
#' object.  Otherwise, if the \code{value} is \code{data.frame}, the row name of
#' the \code{data.frame} should be the sample name, and each column as a
#' metadata variable. 
#' @return A data.frame
#' @examples
#' data(bc_obj)
#'
#' # get the metadata data.frame
#' bc_meta(bc_obj)
#'
#' # assign value to a metadata by $ operation
#' bc_meta(bc_obj)$phenotype <- c("l", "b")
#'
#' # assign value to a metasta by "key" argument
#' bc_meta(bc_obj, key = "sample_type") <- c("l", "b")
#'
#' # show the updated metadata
#' bc_meta(bc_obj)
#'
#' # assign a new data.frame to metadata
#' metadata <- data.frame(
#'     sample_name <- c("test1", "test2"),
#'     phenotype <- c("l", "b")
#'     )
#' rownames(metadata) = bc_names(bc_obj)
#' bc_meta(bc_obj) <- metadata
#' ###
#' @rdname bc_meta
#' @export
setGeneric("bc_meta", function(barcodeObj) { standardGeneric("bc_meta") })

#' @rdname bc_meta
#' @export
setGeneric("bc_meta<-", function(barcodeObj, key=NULL, value) {
    standardGeneric("bc_meta<-") })

#' Transforms BarcodeObj object into other data type
#'
#' Transforms BarcodeObj object into \code{data.frame}, \code{data.table} or
#' \code{matrix}.
#'
#' @param barcodeObj A \code{BarcodeObj} object.
#' @return A \code{data.frame}, with two columns: \code{barcode_seq} and
#' \code{count}.
#' @examples
#' data(bc_obj)
#'
#' bc_obj <- bc_cure_depth(bc_obj)
#'
#' # BarcodeObj to data.frame
#' bc_2df(bc_obj)
#'
#' # BarcodeObj to data.table
#' bc_2dt(bc_obj)
#'
#' # BarcodeObj to matrix
#' bc_2matrix(bc_obj)
#'
#' ###
 
#' @rdname bc_2df
#' @export
setGeneric("bc_2df", function(barcodeObj) { standardGeneric("bc_2df") })

#' @rdname bc_2df
#' @export
setGeneric("bc_2dt", function(barcodeObj) { standardGeneric("bc_2dt") })


#' @rdname bc_2df
#' @export
setGeneric("bc_2matrix", function(barcodeObj) { standardGeneric("bc_2matrix") }) 


# setGeneric("format", function(x) { standardGeneric("format") })
# setGeneric("show", function(object) { standardGeneric("show")})

#' Extract barcode from sequences
#' 
#' \code{bc_extract} identifies the barcodes (and UMI) from the sequences using
#' regular expressions.  \code{pattern} and \code{pattern_type} arguments are
#' necessary, which provide the barcode (and UMI) pattern and their location
#' within the sequences.
#'
#' @param x A single or a list of fastq file, ShortReadQ, DNAStringSet,
#' data.frame, or named integer.
#' @param pattern A string or a string vector with the same number of files,
#' specifying the regular expression with capture. It matchs the barcode (and
#' UMI) with capture pattern.
#' @param sample_name A string vector, applicable when \code{x} is a list or
#' fastq file vector. This argument specifies the sample names. If not provided,
#' the function will look for sample name in the rownames of metadata,
#' the fastqfile name or the \code{list} names.
#' @param metadata A \code{data.frame} with sample names as the row names, with 
#' each metadata record by column, specifying the sample characteristics. 
#' @param maxLDist An integer. The minimun mismatch threshold for barcode
#' matching, when maxLDist is 0, the \code{\link[stringr]{str_match}}  is
#' invoked for barcode matching which is faster, otherwise
#' \code{\link[utils]{aregexec}} is invoked and the \code{costs} parameters can
#' be used to specifying the weight of the distance calculation.
#' @param pattern_type A vector. It defines the barcode (and UMI) capture
#' group. See Details.
#' @param costs A named list, applicable when maxLDist > 0, specifying the
#' weight of each mismatch events while extracting the barcodes.  The list
#' element name have to be \code{sub} (substitution), \code{ins} (insertion) and
#' \code{del} (deletion). The default value is \code{list(sub = 1, ins = 99, del
#' = 99)}.  See \code{\link[utils]{aregexec}} for more detail information.
#' @param ordered A logical value. If the value is true, the return barcodes
#' (UMI-barcode tags) are sorted by the reads counts.
#' @details
#' The \code{pattern} argument is a regular expression, the capture operation
#' \code{()} identifying the barcode or UMI. \code{pattern_type} argument
#' annotates capture, denoting the UMI or the barcode captured pattern. In the
#' example:
#' \preformatted{
#' ([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC
#' |---------| starts with 3 base pairs UMI.
#'            |----------| constant sequence in the backbone.
#'                        |-------| flexible barcode sequences.
#'                                 |---------| 3' constant sequence.
#' }
#'
#' In UMI part \code{[ACGT]{3}}, \code{[ACGT]} means it can be one of
#' the "A", "C", "G" and "T", and \code{{3}} means it repeats 3 times. 
#' In the barcode pattern \code{[ACGT]+}, the \code{+} denotes
#' that there is at least one of the \code{A} or \code{C} or \code{G} or
#' \code{T.}
#' 
#' @return 
#' This function returns a BarcodeObj object if the input is a \code{list} or a
#' \code{vector} of Fastq files, otherwise it returns a \code{data.frame.} In
#' the later case
#' the \code{data.frame} has columns:
#' \enumerate{
#'   \item \code{umi_seq} (optional): UMI sequence, applicable when there is UMI
#'     in `pattern` and `pattern_type` argument.
#'   \item \code{barcode_seq}: barcode sequence.
#'   \item \code{count}: reads number.
#' }
#' 
#'
#' @examples
#' fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
#'
#' library(ShortRead)
#'
#' # barcode from fastq file
#' bc_extract(fq_file, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from ShortReadQ object
#' sr <- readFastq(fq_file)  # ShortReadQ
#' bc_extract(sr, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from DNAStringSet object
#' ds <- sread(sr)  # DNAStringSet
#' bc_extract(ds, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from integer vector
#' iv <- tables(ds, n = Inf)$top # integer vector
#' bc_extract(iv, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from data.frame 
#' df <- data.frame(seq = names(iv), freq = as.integer(iv)) # data.frame
#' bc_extract(df, pattern = "AAAAA(.*)CCCCC")
#'
#' # barcode from list of DNAStringSet
#' l <- list(sample1 = ds, sample2 = ds) # list
#' bc_extract(l, pattern = "AAAAA(.*)CCCCC")
#'
#' # Extract UMI and barcode
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'     )
#'   ) 
#' # barcode backbone with UMI and barcode
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_extract(
#'     list(test = d1), 
#'     pattern, 
#'     sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' ###
#' @rdname bc_extract
#' @export
setGeneric("bc_extract", 
    function(
        x,
        pattern="",
        sample_name=NULL,
        metadata=NULL,
        maxLDist=0,
        pattern_type=c(barcode=1),
        costs=list(sub=1, ins=99, del=99),
        ordered=TRUE
        ) { standardGeneric("bc_extract") })

#' Finds barcode count cutoff point
#'
#' Finds the cutoff point for the barcode count filtering based on the barcode
#' count distribution.
#' 
#' @param barcodeObj A \code{BarcodeObj} object.
#' @param useCleanBc A logical value, if \code{TRUE}, the \code{cleanBc} slot
#' in the \code{BarcodeObj} object will be used, otherwise the \code{messyBc}
#' slot will be used.
#' @return a numeric \code{vector} of the cutoff point.
#'
#' @details The one dimension kmeans clustering is applied for identify the 
#' "true barcode" based on read count. The the algorithm detail is:
#' 1. Remove the barcodes with count below the median of counts.
#' 2. Transform the count by log2(x+1).
#' 3. Apply the 1 dimension clustering to the logarized count, with
#' the cluster number of 2 and weights of the logarized count.
#' 4. Choose the minimum count value in the cluster with more count as
#' cutoff point.
#'
#' For more info about 1 dimension kmeans used here please refer to
#' \code{\link[Ckmeans.1d.dp]{Ckmeans.1d.dp}}.
#' @examples
#' 
#' data(bc_obj)
#' 
#' bc_auto_cutoff(bc_obj)
#' 
#' @rdname bc_auto_cutoff
#' @export
setGeneric("bc_auto_cutoff", 
    function(barcodeObj, useCleanBc=TRUE) { standardGeneric("bc_auto_cutoff") })

#' Filters barcodes by counts
#'
#' bc_cure_depth filters barcodes by the read counts or the UMI counts.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A numeric or a vector of numeric, specifying the threshold of
#' minimum count for a barcode to kept. If the input is a vector, if the vector
#' length is not the same to the sample number the element will be repeatedly
#' used. And when the depth argument is a number with negative value, automatic
#' cutoff point will be chosen by \code{bc_auto_cutoff} function for each
#' samples. See \code{\link[CellBarcode]{bc_auto_cutoff}} for details.
#' @param isUpdate A logical value. If TRUE, the \code{cleanBc} slot in
#' \code{BarcodeObj} will be used preferentially, otherwise the \code{messyBc}
#' slot will be used. If no cleanBc is available, \code{messyBc} will be used.
#' @return A \code{BarcodeObj} object with \code{cleanBc} slot updated or
#' created.
#'
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#'
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"),
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth < 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' bc_2matrix(bc_cured)
#'
#' # Use UMI information, filter the barcode < 5 UMI
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#'
#' ###
#' @rdname bc_cure_depth
#' @export
setGeneric("bc_cure_depth", 
    function(
        barcodeObj,
        depth=0,
        isUpdate=TRUE
        ) { standardGeneric("bc_cure_depth") })

#' Clean barcodes by editing distance
#'
#' \code{bc_cure_cluster} performs clustering of barcodes by editing distance,
#' and remove the minority barcodes with similar sequence. This function is only
#' applicable for the BarcodeObj object with a \code{cleanBc} slot. The barcodes
#' with smaller reads count will be removed.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param dist_threshold A single integer, or vector of integers with the length of
#' sample number, specifying the editing distance threshold for defining two
#' similar barcode sequences. If the input is a vector, each value in the vector
#' relates to one sample according to its order in \code{BarcodeObj} object.
#' The sequences with editing distance equal or less than the threshold will be
#' considered similar barcodes.
#' @param depth_fold_threshold A single numeric or vector of numeric with the
#' length of sample number, specifying the depth fold change threshold of
#' removing the similar minority barcode. The majority barcode should have at
#' least \code{depth_fold_threshold} times of reads of the similar minotiry
#' one, in order to remove the minority similar barcode. (TODO: more preciouse
#' description)
#' @param dist_method A  character string, specifying the editing distance 
#' used for evaluating barcodes similarity. It can be "hamm" for Hamming
#' distance or "leven" for Levenshtein distance.
#' @param cluster_method A character string specifying the algorithm used to
#' perform the clustering of barcodes. Currently only "greedy" is
#' available, in this case, The most and the least abundant barcode will
#' be used for comparing, the least abundant barcode is preferentially removed. 
#' @param count_threshold An integer, read depth threshold to consider a barcode
#' as a true barcode. If a barcode with count higher than this threshold
#' it will not be removeg, even the barcode is similar to more abundant one.
#' Default is 1e9.
#' @param dist_costs A list, the cost of the events of distance algorithm, 
#' applicable when Levenshtein distance is applied. The
#' names of vector have to be \code{insert}, \code{delete} and \code{replace},
#' specifying the weight of insertion, deletion, replacement events
#' respectively. The default cost for each event is 1.
#' @return A BarcodeObj object with cleanBc slot updated.
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'     seq = c(
#'         "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AATTCGATCGATCGAAGAGATCGATCGATC",
#'         "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'         "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'         "AAATCGATCGATCGAAGAGATCGATCGATC",
#'         "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'         "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'         "GGATCGATCGATCGAAGAGATCGATCGATC",
#'         "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'         "GGTTCGATCGATCGACGAGATCGATCGATC",
#'         "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'         ),
#'     freq = c(
#'         30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'         )
#'     )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Remove barcodes with depth < 5
#' (bc_cured <- bc_cure_depth(bc_obj, depth=5))
#' 
#' # Do the clustering, remove the less abundent barcodes
#' # one by hamming distance <= 1 
#' bc_cure_cluster(bc_cured, dist_threshold = 1)
#' 
#' # Levenshtein distance <= 1
#' bc_cure_cluster(bc_cured, dist_threshold = 2, dist_method = "leven",
#'     dist_costs = list("insert" = 2, "replace" = 1, "delete" = 2))
#' 
#' ###
#' @rdname bc_cure_cluster
#' @export
setGeneric("bc_cure_cluster", 
    function(
        barcodeObj,
        dist_threshold=1,
        depth_fold_threshold=1,
        dist_method="hamm",
        cluster_method="greedy",
        count_threshold=1e9,
        dist_costs=list("replace"=1, "insert"=1, "delete"=1)
        ) { standardGeneric("bc_cure_cluster") })

#' Filters UMI-barcode tag by counts 
#'
#' When the UMI is applied, \code{bc_cure_umi} can filter the UMI-barcode tags
#' by counts. 
#'
#' @param barcodeObj A BarcodeObj object.
#' @param depth A numeric or a vector of numeric, specifying the UMI-barcode
#' tag count threshold. Only the barcodes with UMI-barcode tag count equal or larger than
#' the threshold are kept. 
#' @param doFish A logical value, if true, for barcodes with UMI read depth
#' above the threshold, “fish” for identical barcodes with UMI read depth below
#' the threshold. The consequence of \code{doFish} will not increase the number
#' of identified barcodes, but the UMI counts will increase due to including the
#' low depth UMI barcodes. 
#' @param isUniqueUMI A logical value. In the case that a UMI
#' relates to several barcodes, if you believe that the UMI is absolute unique,
#' then only the UMI-barcodes tags with highest count are kept for each UMI.
#' @return A \code{BarcodeObj} object with \code{cleanBc} slot updated (or
#' created).
#' @details When invoke this function, it processes the data with following
#' steps:
#' \enumerate{
#'   \item (if isUniqueUMI is TRUE) Find dominant UMI-barcode tag with highest
#'   reads count in each UMI.
#'   \item UMI-barcode depth filtering.
#'   \item (if doFish is TRUE) Fishing the UMI-barcode tags with low reads
#'   count.
#' }
#'
#' @examples
#' data(bc_obj)
#'
#' d1 <- data.frame(
#'    seq = c(
#'        "ACTTCGATCGATCGAAAAGATCGATCGATC",
#'        "AATTCGATCGATCGAAGAGATCGATCGATC",
#'        "CCTTCGATCGATCGAAGAAGATCGATCGATC",
#'        "TTTTCGATCGATCGAAAAGATCGATCGATC",
#'        "AAATCGATCGATCGAAGAGATCGATCGATC",
#'        "CCCTCGATCGATCGAAGAAGATCGATCGATC",
#'        "GGGTCGATCGATCGAAAAGATCGATCGATC",
#'        "GGATCGATCGATCGAAGAGATCGATCGATC",
#'        "ACTTCGATCGATCGAACAAGATCGATCGATC",
#'        "GGTTCGATCGATCGACGAGATCGATCGATC",
#'        "GCGTCCATCGATCGAAGAAGATCGATCGATC"
#'        ),
#'    freq = c(
#'        30, 60, 9, 10, 14, 5, 10, 30, 6, 4 , 6
#'        )
#'    )
#' 
#' pattern <- "([ACTG]{3})TCGATCGATCGA([ACTG]+)ATCGATCGATC"
#' bc_obj <- bc_extract(list(test = d1), pattern, sample_name=c("test"), 
#'     pattern_type=c(UMI=1, barcode=2))
#'
#' # Use UMI information to remove the barcode <= 5 UMI-barcode tags
#' bc_umi_cured <- bc_cure_umi(bc_obj, depth =0, doFish=TRUE, isUniqueUMI=TRUE)
#' bc_cure_depth(bc_umi_cured, depth = 5)
#'
#' @rdname bc_cure_umi
#' @export
setGeneric("bc_cure_umi", function(
        barcodeObj,
        depth=2,
        doFish=FALSE,
        isUniqueUMI=FALSE
        ) { standardGeneric("bc_cure_umi") })

#' Summary and evaluate barcode diversity
#' 
#' \code{bc_summary_barcode} evaluates sequence diversity metrics using the
#' barcodes data in the \code{cleanBc} slot of \code{BarcodeObj} object. It
#' also generates Lorenz curve and barcode frequency distribution graphs.
#'
#' @param barcodeObj A BarcodeObj object.
#' @param plot A logical value, if TRUE, draw the Lorenz curve and barcode
#' distribution graphs. 
#' @param log_x A logical value, if TRUE, the \code{x} axis is logarized.
#' @return A data.frame with following columns:
#' \itemize{
#'   \item \code{total_reads}: total read number.
#'   \item \code{uniq_barcode}: how many barcodes in the dataset.
#'   \item \code{shannon_index}: Shannon's diversity index or Shannon–Wiener
#'     index.
#'   \item \code{equitability_index}: Shannon's equitability.
#'   \item \code{bit_index}: Shannon bit information.
#' }
#'
#' @details
#' Followings are the metrics used for evaluating the barcode diversity:
#'
#' \emph{Richness}: The unique barcodes number \eqn{R}, it evaluates the
#' richness of the barcodes.
#'
#' \emph{Shannon index}: Shannon diversity index is weighted geometric
#' average of the proportion \eqn{p} of barcodes.
#' \deqn{ H' = - \sum_{i=1}^{R}p_ilnp_i }
#'
#' \emph{Equitability index}: Shannon equitability \eqn{E_H} characterize the
#' evenness of the barcodes, it is a value between 0 and 1, with 1 being
#' complete evenness.
#' \deqn{ E_H = H' / H'_{max} = H / ln(R) }
#'
#' \emph{Bit}:
#' Shannon entropy \eqn{H}, with a units of bit, 
#' \deqn{ H = - \sum_{i=1}^{R}p_ilog_2p_i }
#'
#' @examples
#' data(bc_obj)
#' 
#' # filter barcode by depth
#' bc_obj <- bc_cure_depth(bc_obj)
#'
#' # Output the summary of the barcodes
#' bc_summary_barcode(bc_obj)
#' @rdname bc_summary_barcode
#' @export
setGeneric("bc_summary_barcode", function(
        barcodeObj,
        plot=TRUE,
        log_x=TRUE
        ) { standardGeneric("bc_summary_barcode") })

#' Barcode read count 2D scatter plot of sample combination
#'
#' Draw barcode count scatter plot for all pairwise combination of samples
#' within a \code{BarcodeObj} object. It uses \code{cleanBc} slot in the
#' \code{BarcodeObj} object is used to draw the figure. If the \code{BarcodeObj}
#' object does not have a cleanBc slot, you have to run the \code{bc_cure*}
#' functions in ahead, such as \code{\link[CellBarcode]{bc_cure_depth}},
#' \code{\link[CellBarcode]{bc_cure_umi}}. 
#'
#' @param barcodeObj A \code{BarcodeObj} object, which has a \code{cleanBc} slot
#' @param count_marks A numeric or numeric vector, specifying the read count
#' cutoff in the scatter plot for each sample.
#' @param highlight A character vector, specifying the barcodes to be
#' highlighted.
#' @param log_coord A logical value, if TRUE (default), the \code{x} and
#' \code{y} coordinates of the scatter plot will be logarized by \code{log10.}
#' @param alpha A numeric between 0 and 1, specifies the transparency of the
#' dots in the scatter plot.
#' @return A scatter plot matrix.
#'
#' @examples
#'
#' data(bc_obj)
#'
#' bc_plot_mutual(barcodeObj=bc_obj, count_marks=c(30, 20))
#' ###
#' @rdname bc_plot_mutual
#' @export
setGeneric("bc_plot_mutual", function(
        barcodeObj,
        count_marks=NULL,
        highlight=NULL,
        log_coord=TRUE,
        alpha=0.7
        ) { standardGeneric("bc_plot_mutual") })


#' Scatter plot of barcode count distribution per sample
#'
#' Draws barcode count distribution for each sample in a
#' BarcodeObj object.
#'
#' @param barcodeObj A \code{BarcodeObj} object has a cleanBc slot
#' @param sample_names A character vector or integer vector, specifying the
#' samples used for plot.
#' @param count_marks A numeric or numeric vector, specifying the read count
#' cutoff in the scatter plot for each sample.
#' @param highlight A character vector, specifying the barcodes need to be
#' highlighted.
#' @param log_coord A logical value, if TRUE (default), the \code{x} and
#' \code{y} coordinates of the scatter plot will be logarized by log10.
#' @param alpha A numeric between 0 and 1, specifies the transparency of the
#' dots in the scatter plot.
#' @return 1D distribution graph matrix.
#'
#' @examples
#' data(bc_obj) 
#'
#' bc_plot_single(bc_obj, count_marks=c(10, 11))
#' ###
#' @rdname bc_plot_single
#' @export
setGeneric("bc_plot_single", function(
        barcodeObj,
        sample_names=NULL,
        count_marks=NULL,
        highlight=NULL,
        log_coord=TRUE,
        alpha=0.7
        ) { standardGeneric("bc_plot_single") })


#' Barcode read count 2D scatter plot for given pairs
#'
#' Draws scatter plot for barcode read count between given pairs of samples with
#' a \code{BarcodeObj} object. This function will return scatter plot matrix
#' contains the scatter plots for all given sample pairs.
#'
#' @param barcodeObj A \code{BarcodeObj} object.
#' @param sample_x A character vector or a integer vector, specifying the sample
#' in \code{x} axis of each scatter plot. It can be the sample names in
#' BarcodeObj or the sample index value.
#' @param sample_y A character vector or a integer vector, similar to
#' \code{sample_x}, specifying the samples used for \code{y} axis. It can be the
#' sample names or the
#' sample index value.  
#' @param count_marks_x A numeric vector used to mark the cutoff
#' point for samples in x
#' axis
#' @param count_marks_y A number vector used to mark the cutoff point for
#' samples in y axis.
#' @param highlight A character vector, specifying the barcodes need to be
#' highlighted.
#' @param log_coord A logical value, if TRUE (default), the \code{x} and
#' \code{y} coordinates
#' of the scatter will be logarized by log10.
#' @param alpha A numeric between 0 and 1, specifies the transparency of the
#' dots in the scatter plot.
#' @return Scatter plot matrix.
#' @examples
#' 
#' data(bc_obj)
#'
#' bc_names(bc_obj)
#'
#' bc_plot_pair(barcodeObj=bc_obj, sample_x="test1", sample_y="test2",
#'     count_marks_x=30, count_marks_y=20)
#' ###
#' @rdname bc_plot_pair
#' @export
setGeneric("bc_plot_pair", function(
        barcodeObj,
        sample_x,
        sample_y,
        count_marks_x=NULL,
        count_marks_y=NULL,
        highlight=NULL,
        log_coord=TRUE,
        alpha=0.7
        ) { standardGeneric("bc_plot_pair") })



#' Evaluates sequences quality
#'
#' \code{bc_seq_qc} evaluates sequences quality. See the return value for detail.
#' 
#' @param x A single or list of Fastq file, ShortReadQ object, DNAStringSet
#' object, data.frame or named integer vector.
#' @param sample_name A character vector with the length of sample number, used
#' to set the sample name.
#' @param reads_sample_size A integer value define the sample size of the
#' sequences for quality control analysis. If the there are less sequences comparing
#' to this value, all the sequences will be used. The default is 1e5.
#' @return A barcodeQc or a barcodeQcSet class. 
#' The barcodeQc is a list with four slots, 
#' \itemize{
#'   \item \code{top}: a \code{data.frame} with top 50 most frequency sequence, 
#'   \item \code{distribution}: a \code{data.frame} with the distribution of
#'     read depth. It contains \code{nOccurrences} (depth), and \code{nReads}
#'     (unique sequence) columns.
#'   \item \code{base_quality_per_cycle}: \code{data.frame} with base-pair
#'     location (NGS sequencing cycle) by row, and the base-pair quality summary
#'     by column, including Mean, P5 (5% quantile), P25 (25% quantile), Median,
#'     P75 (75% quantile) and P95 (95% quantile).
#'   \item \code{base_freq_per_cycle}: \code{data.frame} with three columns: 1.
#'     \code{Cycle}, the sequence base-pair location (NGS sequencing cycle); 2.
#'     \code{Base}, DNA base;
#'     \code{Count}: reads count.
#'   \item{summary}: a numeric vector with following elements:
#'     \code{total_read}, \code{median_read_length},
#'     \code{p5_read_length}, \code{p95_read_length}.
#' }
#' The barcodeQcSet is a list of barcodeQc.
#' 
#' @examples
#' library(ShortRead)
#' # fastq file
#' fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
#' bc_seq_qc(fq_file)
#'
#' # ShortReadQ
#' sr <- readFastq(fq_file[1])
#' bc_seq_qc(sr)
#'
#' # DNAStringSet
#' ds <- sread(sr)
#' bc_seq_qc(ds)
#'
#' # List of DNAStringSet
#' l <- list(sample1 = ds, sample2 = ds)
#' bc_plot_seqQc(bc_seq_qc(l))
#'
#' # List of ShortRead
#' l_sr <- list(sample1 = sr, sample2 = sr)
#' bc_plot_seqQc(bc_seq_qc(l_sr))
#'
#' ###
#' @rdname bc_seq_qc
#' @export
setGeneric("bc_seq_qc", function(x, sample_name=NULL, reads_sample_size = 1e5) {
    standardGeneric("bc_seq_qc") })

#' @rdname bc_seq_qc
#' @export
setGeneric("bc_plot_seqQc", function(x) { standardGeneric("bc_plot_seqQc") })

#' Summary barcodeQcSet
#'
#' Summary the "total read count" and "read length" of each samples within
#' a \code{BarcodeQcSet} object, and output a \code{data.frame} with sample by
#' row and different metrics by column.
#'
#' @param x a barcodeQcSet object.
#' @return A \code{data.frame} with 5 columns: \code{sample_name},
#' \code{total_read}, \code{median_read_length}, \code{p5_read_length} and
#' \code{p95_read_length.}
#' 
#' @examples
#'
#' fq_file <- dir(
#'     system.file("extdata", "mef_test_data", package = "CellBarcode"),
#'     full=TRUE)
#'
#' bc_summary_seqQc(bc_seq_qc(fq_file))
#' ###
#' @rdname bc_summary_seqQc
#' @export
setGeneric("bc_summary_seqQc", function(x) { standardGeneric("bc_summary_seqQc") })


#' Remove low quality sequence
#' 
#' Remove low quality sequences by base-pair quality, sequence length or unknown
#' base "N".
#'
#' @param x A single or a list of Fastq file, \code{ShortReadQ},
#' \code{DNAStringSet}, \code{data.frame}, integer vector.
#' @param min_average_quality A numeric or a vector of numeric, specifying the
#' threshold of the minimum average base quality of a sequence to be kept. 
#' @param min_read_length A single or a vector of integer, specifying the
#' sequence length threshold.
#' @param N_threshold A integer or a vector of integer, specifying the maximum
#' \code{N} can be in a sequence.
#' @param sample_name A string vector, specifying the sample name in the output.
#' @return A ShortReadQ or DNAStringSet object with sequences passed the
#' filters.
#' @examples
#' library(ShortRead)
#' 
#' fq_file <- system.file("extdata", "simple.fq", package="CellBarcode")
#'
#' # apply filter to fastq files
#' bc_seq_filter(fq_file)
#'
#' # read in fastq files to get ShortReadQ object
#' sr <- readFastq(fq_file[1])
#' # apply sequencing quality filter to ShortReadQ
#' bc_seq_filter(sr)
#'
#' # get DNAStringSet object
#' ds <- sread(sr)
#' # apply sequencing quality filter to DNAStringSet
#' bc_seq_filter(ds)
#'
#' ###
#' @rdname bc_seq_filter
#' @export
setGeneric("bc_seq_filter", function(
        x,
        min_average_quality=30,
        min_read_length=0,
        N_threshold=0,
        sample_name="") { standardGeneric("bc_seq_filter") })





