## The function returns the index of non-null values
## if noEmpty = TRUE, return 0L when the index is empty  
which.null<-function(..., noEmpty = FALSE){
    args<-list(...)
    res <- vapply(args,is.null,logical(1))
    index <- which(res)
    if(noEmpty&&length(index)==0){
        index = 0L
    }
    index
}
which.NonNull<- function(...,noEmpty = FALSE){
    args<-list(...)
    res <- vapply(args,is.null,logical(1))
    index <- which(!res)
    if(noEmpty&&length(index)==0){
        index = 0L
    }
    index
}

which.firstNull <- function(...,noEmpty = FALSE){
    res <- which.null(...,noEmpty = TRUE)
    if(noEmpty||res[1]!=0){
        res[1]
    }else{
        integer(0)
    }
}

which.firstNonNull<- function(...,noEmpty = FALSE){
    res <- which.NonNull(...,noEmpty = TRUE)
    if(noEmpty||res[1]!=0){
        res[1]
    }else{
        integer(0)
    }
}

all.null<-function(...){
    args<-list(...)
    res <- vapply(args,is.null,logical(1))
    all(res)
}

getIndexSimple<-function(n,alpha0,index){
    if(!is.null(n)&&!is.null(alpha0)&&is.null(index)){
        nRegion <- max(floor(alpha0 * n), 1)
        index <- seq_len(nRegion)
    }
    index
}

getIndex <- function(n,alpha0,index=NULL,indexL=NULL,indexU=NULL){
    if(!is.null(index)&&!all.null(indexL,indexU)){
        stop("The argument `index` and arguments `indexL`, `indexU`
             cannot be both non-null")
    }
    firstNonNull <- which.firstNonNull(indexL,indexU,index,alpha0,noEmpty=TRUE)
    if(firstNonNull==0L){
        stop("The range argument `alpha0` and `index` are both null")
    }else if(firstNonNull==4L){
        if(is.null(n))stop("The sample size is missing")
        nRegion <- max(floor(alpha0 * n), 1)
        index <- seq(1, nRegion)
        indexL <- seq(1, nRegion)
        indexU <- seq(1, nRegion)
    }else if(firstNonNull==3L){
        indexL <- index
        indexU <- index
    }
    list(n = n, alpha0=alpha0,index =index,indexL=indexL,indexU=indexU)
}


getArgs<-function(stat,n=NULL,alpha0=NULL,index=NULL,indexL=NULL,indexU=NULL){
    if(!is.generalKSStat(stat)){
        args <- getIndex(n=n,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU)
        args[["statValue"]] <- stat
        return(args)
    }
    if(is.null(n)){
        n <- stat[["n"]]
    }
    if(all.null(alpha0,index,indexL,indexU)){
        alpha0<-stat[["alpha0"]]
        index<-stat[["index"]]
        indexU<-stat[["indexU"]]
        indexL<-stat[["indexL"]]
    }
    if(length(grep("+",stat$statName,fixed=TRUE))==1&&is.null(indexL)){
        nRegion <- max(floor(alpha0 * n), 1)
        indexL <- seq_len(nRegion)
    }
    if(length(grep("-",stat$statName,fixed=TRUE))==1&&is.null(indexU)){
        nRegion <- max(floor(alpha0 * n), 1)
        indexU <- seq_len(nRegion)
    }
    args <- getIndex(n=n,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU)
    args$statValue <- stat$statValue
    args
}

getStatFullName <-function(statName,indexL, indexU){
    if(!is.null(indexL)&&is.null(indexU))
        statSign <- "+"
    else if(is.null(indexL)&&!is.null(indexU))
        statSign <- "-"
    else 
        statSign <- ""
    paste0(statName,statSign)
}

is.generalKSStat<-function(x){
    is(x,"generalKSStat")
}
