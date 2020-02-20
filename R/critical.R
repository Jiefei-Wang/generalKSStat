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


#' @export
HCCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= HCPvalue,searchRange=c(0,100),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}

#' @export
BJCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= BJPvalue,searchRange=c(0,1),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}

#' @export
KSCritical<-function(alpha,n=NULL,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU= NULL){
    genericCritical(
        pvalueFunc= KSPvalue,searchRange=c(0,1),
        alpha=alpha,n=n,alpha0=alpha0,index=index,
        indexL=indexL,indexU= indexU
    )
}
#' @export
HCPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    HCCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @export
HCMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    HCCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}
#' @export
BJPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    BJCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @export
BJMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    BJCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}
#' @export
KSPlusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    KSCritical(alpha = alpha, n = n, alpha0 = alpha0, indexL = index)
}
#' @export
KSMinusCritical<-function(alpha,n=NULL,alpha0=NULL,index=NULL){
    KSCritical(alpha = alpha, n = n, alpha0 = alpha0, indexU = index)
}