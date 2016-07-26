#' This is a function to merge BS-seq data from multiple replicates.
#'
#' Return read counts and methylated reads counts for sites commonly
#' mapped by all replicates.
#' @param cases A list of data frames for multiple biological replicates
#'       under treatment experiment.
#' @param controls A list of data frames for multiple biological replicates
#'       under control experiment.
#' @return A list of read counts, methylated reads and commonly mapped sites.
#' @examples
#' data("case1", "case2", "control1", "control2")
#' read_counts <- seqdata(list(case1, case2), list(control1, control2))
#' @export
seqdata<-function(cases, controls){
    n1<-length(cases)
    n2<-length(controls)

    #if only one sample in each condition
    if(n1==1 & n2==1){
        case_1<-cases[[1]]
        control_1<-controls[[1]]
        rownames(case_1)<-case_1[,1]
        rownames(control_1)<-control_1[,1]
        base_case<-case_1[,1]
        base_control<-control_1[,1]
        base<-intersect(base_case, base_control)
        case1<-case_1[base,]
        control1<-control_1[base, ]
        L<-length(base)
        C<-matrix(0, nrow=sum(n1,n2), ncol=L)
        C[1,]<-case1[, 5]
        C[2,]<-control1[,5]
        p<-matrix(0, nrow=sum(n1,n2), ncol=L)
        p[1,]<-case1[,6]
        p[2,]<-control1[,6]
        M<-round(p*C/100)

    } else if(n1>1 | n2>1){
        case_1<-cases[[1]]
        base_case<-case_1[,1]
        if(n1>1){
            for (i in 2:n1){
                case_oth<-cases[[i]]
                base_case_oth<-case_oth[,1]
                base_case<-intersect(base_case,base_case_oth)
            }
        }
        control_1<-controls[[1]]
        base_control<-control_1[,1]
        if(n2>1){
            for(i in 2:n2){
                control_oth<-controls[[i]]
                base_control_oth<-control_oth[,1]
                base_control<-intersect(base_control,base_control_oth)
            }
        }
        base<-intersect(base_case, base_control)
        L<-length(base)
        C<-matrix(0,nrow=sum(n1,n2), ncol=L)
        p<-matrix(0,nrow=sum(n1,n2), ncol=L)
        for (i in 1:n1){
            rownames(cases[[i]])<-cases[[i]][,1]
            C[i,]<-cases[[i]][base,5]
            p[i,]<-cases[[i]][base,6]
        }
        for (i in (n1+1):(n1+n2)){
            rownames(controls[[i-n1]])<-controls[[i-n1]][,1]
            C[i,]<-controls[[i-n1]][base,5]
            p[i,]<-controls[[i-n1]][base,6]
        }
        M<-round(C*p/100)
    }

    return(list(C, M, base))

}









