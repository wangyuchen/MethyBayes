#' This is a likelihood function for group 1 (hyper-methylated group)
#'
#' This function calculates the likelihood of the data
#' given the site belonging to group 1
#' @param CDcol a column of read counts under treatment experiment
#' @param CUcol a column of read counts under contrl experiment
#' @param MDcol a column of methylated read counts under treatment experiment
#' @param MUcol a column of methylated read counts under control experiment
#' @param integ1 probability of x less than y (x and y follow beta distribtuion
#'        with parameters calculated from the data in treatment and control respectively
#'        please see paper for detail information)
#' @param alpha2 parameter value for likelihood in group 1
#' @param beta2 parameter value for likelihood in group 1
#' @return the likelihood of the data given the site belonging to group 1
#'


lklihd2<-function(CDcol,CUcol,MDcol,MUcol,integ1, alpha2, beta2){
  NDcol<-CDcol-MDcol
  NUcol<-CUcol-MUcol
  l1<-lbeta(sum(MDcol)+alpha2,sum(NDcol)+beta2)-lbeta(alpha2,beta2)
  l2<-lbeta(sum(MUcol)+alpha2,sum(NUcol)+beta2)-lbeta(alpha2,beta2)
  #l1<-lbeta(sum(MDcol)+alpha2,sum(NDcol)+beta2)
  #l2<-lbeta(sum(MUcol)+alpha2,sum(NUcol)+beta2)

  Ccol<-c(CDcol,CUcol)
  Mcol<-c(MDcol,MUcol)
  Ncol<-c(NDcol,NUcol)
  l3<-sum(lgamma(Ccol+1))-sum(lgamma(Ncol+1))-sum(lgamma(Mcol+1))
  #l3<--sum(lgamma(Ncol+1))-sum(lgamma(Mcol+1))
  prob1<-integ1
  prob1[prob1<0.00001]<-0.00001
  l4<-log(2)+log(prob1)
  l1+l2+l3+l4
}
