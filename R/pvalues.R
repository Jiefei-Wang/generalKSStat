#' Compute the pvalue of generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the pvalue of generalized Kolmogorov-Smirnov test statistics,
#' see details on how to use the function.
#' 
#' @param x Numeric, the sample data.
#' @param stat a `generalKSStat` object or numeric, the statistic that the p-value
#' is computed for. This parameter will be ignored if the parameter `x` is not null.
#' If `stat` is a `generalKSStat` object, there is no need to 
#' specify the other parameters unless you want to overwrite them.
#' If a numeric value is provided to the parameter `stat`, you must at least 
#' specify the sample size `n`.
#' @param n Integer, the sample size of the data.
#' 
#' @examples 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' 
#' ## Perform KS test
#' ks_res <- GKSStat(x = x, statName = "KS")
#' 
#' ## Compute the pvalue for the KS test
#' GKSPvalue(stat = ks_res)
#' 
#' ## For any observed statistic
#' GKSPvalue(stat = 0.2, n = 10, statName = "KS")
#' 
#' ## Change the detection range of the KS test
#' ## to test only the first 3 ordered samples
#' ## All gives the same result
#' GKSPvalue(stat = 0.2, n = 10, alpha0 = 0.3, statName = "KS")
#' GKSPvalue(stat = 0.2, n = 10, index = 1:3, statName = "KS")
#' GKSPvalue(stat = 0.2, n = 10, indexL = 1:3, indexU = 1:3, statName = "KS")
#' 
#' 
#' 
#' @return A numeric value representing the pvalue
#' @inheritParams GKSStat
#' @inherit GKSStat details
#' @rdname pvalue
#' @export
GKSPvalue<-function(stat=NULL , n =NULL, alpha0 = NULL, 
                     index=NULL,indexL=NULL,indexU=NULL,
                     x=NULL, statName = NULL){
    
    if(is.generalKSStat(stat)){
        statName <- stat$statName
        statValue <- stat$statValue
        n <- stat$n
        indexL <- stat$indexL
        indexU <- stat$indexU
    }else{
        statName <- match.arg(statName,
                              c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-"))
        statValue <- stat
        if(!is.null(x)){
            stat <- GKSStat(x=x,alpha0=alpha0,index=index,
                            indexL=indexL,indexU=indexU,
                            statName = statName,pvalue=TRUE)
            return(getPvalue(stat))
        }
        
    }
    stopifnot(!is.null(n))
    sideIndex <- getTwoSideIndex(statName=statName,
                                 n=n,
                                 alpha0=alpha0,
                                 index=index,
                                 indexL=indexL,
                                 indexU=indexU)
    indexL <- sideIndex$indexL
    indexU <- sideIndex$indexU
    pvalue <- call_func(root = "Pvalue",prefix = sideIndex$statName,
                        stat=statValue,n=n,indexL=indexL,indexU=indexU)
    
    return(pvalue)
}




genericPvalue<-function(statName,statValue ,n=n,indexL,indexU){
    localCritical <- call_func(root = statName,postfix = "LocalCritical",
              stat=statValue,n=n)
    l=localCritical$l
    h=localCritical$h
    if(length(indexL)!=0){
        l[-indexL] <- 0
    }else{
        l=rep(0,length(l))
    }
    if(length(indexU)!=0){
        h[-indexU] <- 1
    }else{
        h=rep(1,length(h))
    }
    orderedProb(l,h)
}



#' @rdname pvalue
#' @export
HCPvalue<-function(stat,n,indexL=NULL,indexU=NULL){
    res=1-genericPvalue("HC",statValue=stat,n=n,indexL=indexL,indexU =indexU)
}


#' @rdname pvalue
#' @export
BJPvalue<-function(stat,n,indexL=NULL,indexU=NULL){
    res=1-genericPvalue("BJ",statValue=stat,n=n,indexL=indexL,indexU =indexU)
}

#' @rdname pvalue
#' @export
KSPvalue<-function(stat,n,indexL=NULL,indexU=NULL){
    res=1-genericPvalue("KS",statValue=stat,n=n,indexL=indexL,indexU =indexU)
}
