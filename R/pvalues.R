GKSPvalue<-function(stat , n =NULL, alpha0 = NULL, 
                     index=NULL,indexL=NULL,indexU=NULL,
                     statName = NULL){
    if(is.null(statName)){
        if(is(stat,"generalKSStat")){
            statName=stat$statName
        }else{
            statName="KS"
        }
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


#' @export
HCPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    HCPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @export
HCMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    HCPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}
#' @export
BJPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    BJPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @export
BJMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    BJPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}
#' @export
KSPlusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    KSPvalue(stat = stat, n=n, alpha0=alpha0,
             indexL= index)
}
#' @export
KSMinusPvalue<-function(stat,n=NULL,alpha0=NULL,index=NULL){
    KSPvalue(stat = stat, n=n, alpha0=alpha0,
             indexU= index)
}