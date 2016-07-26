#' Example datasets
#'
#' Simulated example datasets to illustrate the format of data and test the
#' basic functions.
#'
#' These are the datasets of cases and controls, each have two replicates. Case
#' 1 and case 2 have 4714 and 5197 observations while control 1 and control 2
#' have 1817 and 5249 observations respectively.
#'
#' @format Both of the data sets are data frame with 7 variables:
#' \describe{
#'   \item{chrBase}{Factor with 4714 levels, unique ID.}
#'   \item{chr}{Factor with 1 level "chr1", Chromosome name.}
#'   \item{base}{Integer of base position.}
#'   \item{strand}{Factor with 2 levels "F" and "R" for strand.}
#'   \item{coverage}{Integer of read coverage, i.e, total number of reads.}
#'   \item{freqC}{Numeric percentage (0 - 100) of C bases on that location,
#'                i.e. percentage of reads showing methylation.}
#'   \item{freqT}{Numeric percentage of T bases on that location, i.e.
#'                percentage of reads showing un-methylation.}
#' }
#'
#' @rdname demo_data
"case1"

#' @rdname demo_data
"case2"

#' @rdname demo_data
"control1"

#' @rdname demo_data
"control2"
