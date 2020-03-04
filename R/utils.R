## The cached result
cache <- list()
cache$criticals <- new.env()

call_func <- function(root, prefix=NULL, postfix=NULL, ...){
    func_name <- paste0(c(prefix,root,postfix),collapse="")
    do.call(func_name,args = list(...))
}

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

getIndexOneSide<-function(n,alpha0,index,indexOneSide){
    if(!is.null(indexOneSide)){
        return(indexOneSide)
    }else{
        if(!is.null(index)){
            return(index)
        }else{
            if(!is.null(n)&&!is.null(alpha0)){
                nRegion <- max(floor(alpha0 * n), 1)
                index <- seq_len(nRegion)
                return(index)
            }
        }
    }
    NULL
}

getTwoSideIndex <- function(statName,n,alpha0,index,indexL,indexU){
    side <- substring(statName,nchar(statName))
    oneSideInd = which(side==c("+","-"))
    if(length(oneSideInd)!=0){
        statName <- substr(statName,1,nchar(statName)-1)
        if(oneSideInd==1L){
            indexL <- getIndexOneSide(n,alpha0,index,indexL)
            if(is.null(indexL))
                indexL <- seq_len(n)
            indexU <- NULL
        }else{
            indexU <- getIndexOneSide(n,alpha0,index,indexU)
            if(is.null(indexU))
                indexU <- seq_len(n)
            indexL <- NULL
        }
    }else{
        indexL <- getIndexOneSide(n,alpha0,index,indexL)
        indexU <- getIndexOneSide(n,alpha0,index,indexU)
        if(all.null(indexL,indexU)){
            indexL <- seq_len(n)
            indexU <- seq_len(n)
        }
    }
    if(!is.null(n)){
        stopifnot(length(indexL)<=n)
        stopifnot(length(indexU)<=n)
    }
    list(statName = statName, side =side, indexL= indexL,indexU=indexU)
}

# 
# getStatFullName <-function(statName,indexL, indexU){
#     if(!is.null(indexL)&&is.null(indexU))
#         statSign <- "+"
#     else if(is.null(indexL)&&!is.null(indexU))
#         statSign <- "-"
#     else 
#         statSign <- ""
#     paste0(statName,statSign)
# }

is.generalKSStat<-function(x){
    is(x,"generalKSStat")
}
