#' This is a function to produce methylation identification result.
#'
#' This function returns the methylation identification result for all sites
#' commonly mapped by all replicates under two experimental conditions.
#'
#' @param data A list of read counts, methylated read counts and
#'        sites commonly mapped by all replicates.
#' @param n1 Number of replicates in treatment experiment.
#' @param n2 Number of replicates in control experiment.
#' @return Methylation identification result matrix, first column is the list
#' of sites, second column is the list of identification result.
#'
#' @references Wang, H., He, C., Kushwaha, G., Xu, D., \& Qiu, J. (2016).
#' A full Bayesian partition model for identifying hypo-and hyper-methylated
#' loci from single nucleotide resolution sequencing data. BMC bioinformatics,
#' 17(1), 71.
#'
#' @examples
#' \dontrun{
#'  data("case1", "case2", "control1", "control2")
#'  read_counts <- seqdata(list(case1, case2), list(control1, control2))
#'  ident <- MethyBayes(read_counts, 2, 2)
#' }
#' @export
#' @useDynLib MethyBayes
MethyBayes<-function(data, n1, n2){
    C<-data[[1]]
    M<-data[[2]]
    ones<-which(colSums(M==C)==4)
    zeros<-which(colSums(M)==0)
    index<-sort(union(ones, zeros))
    p<-process(data, n1, n2)
    pnew<-p[,-index]
    l<-length(data[[3]])-length(index)
    g0<-round(l/10)
    g1<-round(l*2/10)
    g<-max(l*150L, 1000000L)
    #g<-1000000L
    # g<-l*1000L
    res <- matrix(0L, l, 3)
    sum_counts<- .Fortran("mcmc_froutine", l=as.integer(l), g0=as.integer(g0),
                    g1=as.integer(g1), g=g, p=pnew, res)[[6]]
    Iprop <- apply(sum_counts, MARGIN = 1, which.max) - 1
    result<-rep(0, length(data[[3]]))
    result[-index]<-Iprop
    final<-matrix(NA, length(data[[3]]), 2)
    final[,1]<-data[[3]]
    final[,2]<-result

    return(final)

}
