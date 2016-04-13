#' Example datasets
#'
#' These are the datasets of case and controls, each have two replicates. Case
#' 1 and case 2 have 4714 and 5197 observations and control 1 and control 2
#' have 1817 and 5249 observations respectively.
#'
#' @format Both of the data sets are data frame with 7 variables:
#' \describe{
#'   \item{chrBase}{Factor}
#'   \item{chr}{Factor with 1 level "chr1".}
#'   \item{base}{integer ...}
#'   \item{strand}{Factor with 2 levels "F" and "R".}
#'   \item{coverage}{integer ...}
#'   \item{freqC}{numeric ...}
#'   \item{freqT}{numeric ...}
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
