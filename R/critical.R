#' Compute the critical value of generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the critical value of generalized Kolmogorov-Smirnov test statistics,
#' see details on how to use the function.
#' 
#' @param alpha numeric, the type I error rate for the critical value. Please do 
#' not be confused with `alpha0`.
#' 
#' @examples 
#' ## Compute the critical value of the KS test
#' ## of sample size 10
#' GKSCritical(alpha = 0.05, n = 10, statName = "KS")
#' 
#' ## The critical value for the test that
#' ## only considers the first 3 ordered samples
#' ## All gives the same result.
#' GKSCritical(alpha = 0.05, n = 10, alpha0 = 0.3, statName = "KS")
#' GKSCritical(alpha = 0.05, n = 10, index = 1:3, statName = "KS")
#' GKSCritical(alpha = 0.05, n = 10, indexL = 1:3, indexU = 1:3, statName = "KS")
#' 
#' 
#' @return A critical value
#' @inheritParams GKSStat
#' @inherit GKSStat details
#' @rdname critical
#' @export
GKSCritical <-function(alpha,n=NULL,alpha0=1,
                       index=NULL,indexL=NULL,indexU= NULL,
                       statName = c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-")){
    statName <- match.arg(statName)
    if(statName=="KS"){
        stat <- KSCritical(alpha=alpha,n=n,alpha0=alpha0,
                         index=index,indexL=indexL,indexU=indexU)
        return(stat)
    }
    if(statName=="KS+"){
        stat <- KSPlusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="KS-"){
        stat <- KSMinusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="BJ"){
        stat <- BJCritical(alpha=alpha,n=n,alpha0=alpha0,
                         index=index,indexL=indexL,indexU=indexU)
        return(stat)
    }
    if(statName=="BJ+"){
        stat <- BJPlusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="BJ-"){
        stat <- BJMinusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC"){
        stat <- HCCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC+"){
        stat <- HCPlusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC-"){
        stat <- HCMinusCritical(alpha=alpha,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    
    
    
    
}

genericCritical<-function(pvalueFunc, searchRange,
                          alpha,n=NULL,alpha0=NULL,index=NULL,
                          indexL=NULL,indexU= NULL){
    rootFunc=function(stat) 
        sapply(stat, function(stat)
            pvalueFunc(stat=stat,n=n,alpha0=alpha0,
                   index=index,indexL=indexL,indexU=indexU)-alpha)
    res=uniroot(rootFunc,searchRange,extendInt = "yes")
    res$root
}


#' @rdname critical
#' @export
HCCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= HCPvalue,searchRange=c(0,100),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}

#' @rdname critical
#' @export
BJCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= BJPvalue,searchRange=c(0,1),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}

#' @rdname critical
#' @export
KSCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= KSPvalue,searchRange=c(0,1),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}
#' @rdname critical
#' @export
HCPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    HCCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @rdname critical
#' @export
HCMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    HCCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}
#' @rdname critical
#' @export
BJPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    BJCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @rdname critical
#' @export
BJMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    BJCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}
#' @rdname critical
#' @export
KSPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    KSCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @rdname critical
#' @export
KSMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    KSCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}