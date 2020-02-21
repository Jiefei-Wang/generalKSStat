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
    if(!is.generalKSStat(stat)){
        if(all.null(alpha0,index,indexL,indexU)){
            alpha0 <- 1
        }
    }
    if(is.null(statName)){
        if(is.generalKSStat(stat)){
            statName=stat$statName
        }else{
            statName="KS"
        }
    }
    if(!is.null(x)){
        if(is.null(alpha0)){
            alpha0 <- 1
        }
        stat <- GKSStat(x=x,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU,
                        statName = statName,pvalue=TRUE)
        return(getPvalue(stat))
    }
    
    if(statName=="KS"){
        stat <- KSPvalue(stat=stat,n=n,alpha0=alpha0,
                       index=index,indexL=indexL,indexU=indexU)
        return(stat)
    }
    if(statName=="KS+"){
        stat <- KSPlusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="KS-"){
        stat <- KSMinusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="BJ"){
        stat <- BJPvalue(stat=stat,n=n,alpha0=alpha0,
                       index=index,indexL=indexL,indexU=indexU)
        return(stat)
    }
    if(statName=="BJ+"){
        stat <- BJPlusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="BJ-"){
        stat <- BJMinusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC"){
        stat <- HCPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC+"){
        stat <- HCPlusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    if(statName=="HC-"){
        stat <- HCMinusPvalue(stat=stat,n=n,alpha0=alpha0,index=index)
        return(stat)
    }
    stop("Undefined stat name:",statName)
}




genericPvalue<-function(l,h,indexL,indexU){
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


getC = function(x, a){
    out=(x+(a^2-a*(a^2+4*(1-x)*x)^0.5)/2)/(1+a^2);
    return(out)
}

getHCLower<-function(stat,n){
    stat <- stat/sqrt(n)
    const<-seq_len(n)/n
    a<-1+stat^2
    b<--2*const-stat^2
    (-b-sqrt(b^2-4*a*const^2))/2/a
}
getHCUpper<-function(stat,n){
    stat <- stat/sqrt(n)
    const<-(seq_len(n)-1)/n
    a<-1+stat^2
    b<--2*const-stat^2
    (-b+sqrt(b^2-4*a*const^2))/2/a
}
#' @rdname pvalue
#' @export
HCPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,indexL=NULL,indexU=NULL){
    args <- getArgs(stat=stat,n=n,alpha0=alpha0,
                    index=index,indexL=indexL,indexU=indexU)
    
    n <- args$n
    stat <- args$statValue
    l=getHCLower(stat,n)
    h=getHCUpper(stat,n)
    res=1-genericPvalue(l,h,indexL=args$indexL,indexU =args$indexU)
    res
}


#' @rdname pvalue
#' @export
BJPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,indexL=NULL,indexU=NULL){
    args <- getArgs(stat=stat,n=n,alpha0=alpha0,
                    index=index,indexL=indexL,indexU=indexU)
    n <- args$n
    stat <- args$statValue
    l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
    h=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
    res=1-genericPvalue(l,h,indexL=args$indexL,indexU =args$indexU)
    res
}

#' @rdname pvalue
#' @export
KSPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL,indexL=NULL,indexU=NULL){
    args <- getArgs(stat=stat,n=n,alpha0=alpha0,
                    index=index,indexL=indexL,indexU=indexU)
    n <- args$n
    stat <- args$statValue
    l <- 1:n/n - stat
    h <- stat + 1:n/n-1/n
    l[l<0]=0
    h[h>1]=1
    res=1-genericPvalue(l,h,indexL=args$indexL,indexU =args$indexU)
    res
}


#' @rdname pvalue
#' @export
HCPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    HCPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @rdname pvalue
#' @export
HCMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    HCPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}
#' @rdname pvalue
#' @export
BJPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    BJPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @rdname pvalue
#' @export
BJMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    BJPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}
#' @rdname pvalue
#' @export
KSPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    KSPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @rdname pvalue
#' @export
KSMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    index <- getIndexSimple(n,alpha0,index)
    KSPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}