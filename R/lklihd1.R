#' This is a likelihood function for group 0 (euqal-methylated group)
#'
#' This function calculates the likelihood of data given
#' the site belonging to group 0
#' @param Ccol a column of read counts for all replicates under two conditions
#' @param Mcol a column of methylated read counts for all replicates
#'        under two conditions
#' @param alpha1 parameter value in likelihood for group 0
#' @param beta1 parameter value in likelihood for group 0
#' @return the likelihood of data given the site belonging to group 0
lklihd1<-function(Ccol,Mcol,alpha1,beta1){
  Ncol<-Ccol-Mcol
  l1<-lbeta(sum(Mcol)+alpha1,sum(Ncol)+beta1)-lbeta(alpha1,beta1)
  l2<-sum(lgamma(Ccol+1))-sum(lgamma(Ncol+1))-sum(lgamma(Mcol+1))
  #l1<-lbeta(sum(Mcol)+alpha1,sum(Ncol)+beta1)
  #l2<--sum(lgamma(Ncol+1))-sum(lgamma(Mcol+1))
  l1+l2
}
