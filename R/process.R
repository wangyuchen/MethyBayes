process<-function(data, n1, n2){
    # a utility function to process the read counts and
    # methylated read counts to a new data matrix used for MCMC s
    # @param data a list of read counts, methylated read counts and commonly
    # mapped sites
    # @param n1 number of replicates in treatment experiment
    # @param n2 number of replicates in control experiments
    # @return p matrix used for MCMC
    # set.seed(123123123)
    C<-data[[1]]
    M<-data[[2]]
    #reorganized data
    CD<-matrix(C[1:n1, ],nrow=n1)
    CU<-matrix(C[(n1+1):(n1+n2),], nrow=n2)
    MD<-matrix(M[1:n1, ],nrow=n1)
    MU<-matrix(M[(n1+1):(n1+n2),],nrow=n2)
    ND<-CD-MD
    NU<-CU-MU
    # hyper-parameters set up for beta priors
    alpha1<-1
    beta1<-1
    alpha2<-1
    beta2<-1
    alpha3<-1
    beta3<-1

    #############integration(prob(xi<yi)&prob(xi>yi))###########
    ###xi ~ beta(alp1,bta1);yi~beta(alp2,bta2)########################
    MDclsum<-colSums(MD)
    MUclsum<-colSums(MU)
    NDclsum<-colSums(ND)
    NUclsum<-colSums(NU)

    alp1<-alpha2+MDclsum
    bta1<-beta2+NDclsum
    alp2<-alpha2+MUclsum
    bta2<-beta2+NUclsum
    L<-ncol(C)
    integrat1<-rep(0,L)
    s<-1000
    for(i in 1:L){
        x<- stats::rbeta(s,alp1[i],bta1[i])
        y<- stats::rbeta(s,alp2[i],bta2[i])
        integrat1[i]<-sum(x<y)/s
    }
    integrat2<-1-integrat1

    p<-matrix(NA,nrow=3,ncol=L)
    for(i in 1:L){
        Ccol<-C[,i]
        Mcol<-M[,i]
        CDcol<-CD[,i]
        CUcol<-CU[,i]
        MDcol<-MD[,i]
        MUcol<-MU[,i]
        integ1<-integrat1[i]
        integ2<-integrat2[i]
        p[1,i]<-lklihd1(Ccol,Mcol, alpha1,beta1)
        p[2,i]<-as.numeric(lklihd2(CDcol,CUcol,MDcol,MUcol,integ1,
                                   alpha2,beta2))
        p[3,i]<-as.numeric(lklihd3(CDcol,CUcol,MDcol,MUcol,integ2,
                                   alpha3,beta3))
    }

    return(p)
}
