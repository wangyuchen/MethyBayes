#' This is a function to produce methylation identification result
#' This is a function
#'
#' Returns the methylation identification result for all sites
#' commonly mapped by all replicates under two experimental conditions
#'
#' @param data a list of read counts , methylated read counts and
#'        sites commonly mapped by all replicates
#' @param n1 number of replicates in treatment experiment
#' @param n2 number of replicates in control experiment
#' @return methylation identification result matrix, first
#'         column is the list of sites, second column is the
#'         list of identification result
#'
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
  sum_counts<- .Fortran("mcmc", l=as.integer(l), g0=as.integer(g0), g1=as.integer(g1), g=g, p=pnew, res)[[6]]
  Iprop <- apply(sum_counts, MARGIN = 1, which.max) - 1
  result<-rep(0, length(data[[3]]))
  result[-index]<-Iprop
  final<-matrix(NA, length(data[[3]]), 2)
  final[,1]<-data[[3]]
  final[,2]<-result

  return(final)

}
