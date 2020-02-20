#' @export
GKSStat <- function(
    x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL,
    statName = c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-"),
    pvalue = TRUE){
    statName <- match.arg(statName)
    if(statName=="KS"){
        stat <- KSStat(x=x,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU)
    }
    if(statName=="KS+"){
        stat <- KSPlusStat(x=x,alpha0=alpha0,index=index)
    }
    if(statName=="KS-"){
        stat <- KSMinusStat(x=x,alpha0=alpha0,index=index)
    }
    if(statName=="BJ"){
        stat <- BJStat(x=x,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU)
    }
    if(statName=="BJ+"){
        stat <- BJPlusStat(x=x,alpha0=alpha0,index=index)
    }
    if(statName=="BJ-"){
        stat <- BJMinusStat(x=x,alpha0=alpha0,index=index)
    }
    if(statName=="HC"){
        stat <- HCStat(x=x,alpha0=alpha0,index=index,indexL=indexL,indexU=indexU)
    }
    if(statName=="HC+"){
        stat <- HCPlusStat(x=x,alpha0=alpha0,index=index)
    }
    if(statName=="HC-"){
        stat <- HCMinusStat(x=x,alpha0=alpha0,index=index)
    }
    if(pvalue)
        stat$pvalue <- GKSPvalue(stat)
    stat
}



## Generic statistical function
genericStat<-function(statName,isMin , levelFunc,
                      x,alpha0,
                      index=NULL,indexL=NULL,indexU=NULL){
    localLevels <- partialLevelStat(statFunc = levelFunc,
                                 x = x,
                                 alpha0 = alpha0,
                                 index= index,
                                 indexL = indexL,
                                 indexU = indexU
    )
    if(isMin){
        stat <- min(localLevels)
    }else{
        stat <- max(localLevels)
    }
    .generalKSStat(statName = statName,
               statValue= stat,
               n=length(x),
               alpha0=alpha0,
               index=index,
               indexL=indexL,
               indexU=indexU)
}






#' @export
HCStat<-function(x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="HC",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = FALSE,levelFunc = HCLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @export
BJStat<-function(x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="BJ",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = TRUE,levelFunc = BJLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @export
KSStat<-function(x,alpha0=1,index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="KS",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = FALSE,levelFunc = KSLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @export
HCPlusStat<-function(x,alpha0 = 1, index=NULL){
    HCStat(x=x,alpha0=alpha0,indexL=index)
}
#' @export
HCMinusStat<-function(x,alpha0 = 1, index=NULL){
    HCStat(x=x,alpha0=alpha0,indexU=index)
}

#' @export
BJPlusStat<-function(x,alpha0 = 1, index=NULL){
    BJStat(x=x,alpha0=alpha0,indexL=index)
}
#' @export
BJMinusStat<-function(x,alpha0 = 1, index=NULL){
    BJStat(x=x,alpha0=alpha0,indexU=index)
}

#' @export
KSPlusStat <- function(x,alpha0=1,index=NULL){
    KSStat(x=x,alpha0=alpha0,indexL=index)
}

#' @export
KSMinusStat <- function(x,alpha0=1,index=NULL){
    
    KSStat(x=x,alpha0=alpha0,indexU=index)
}
