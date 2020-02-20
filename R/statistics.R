#' Compute generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the Kolmogorov-Smirnov, Berk-Jones or the higher criticism statistics
#' to test for uniform distributions. The function `GKSStat` provides an uniform 
#' way to computes different test statistics and the functions `HCStat`, `BJStat`
#' and `KSStat` are for a specific statistic. The function `XPlusStat` 
#' and `XMinusStat`provides one-sided version of the corresponding `X` test statistics.
#' To be consistent with the other statistics, the traditional higher criticism 
#' statistic is in the function `HCPlusStat` and the function `HCStat` computes 
#' two-sided higher criticism statistic.
#' 
#' @param x Numeric, the samples that the test statistics will be based on.
#' @param alpha0 Numeric, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param index Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param indexL Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param indexU Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param statName Character, the name of the statistic that will be computed.
#' The default is "KS".
#' @param pvalue Logical, whether to compute the p-value of the statistic. 
#' The default is `TRUE`
#' 
#' @details 
#' \bold{statistics definitions}
#' 
#' The function compute the test statistics which defines the geometrical 
#' distance between the null distribution and empirical distribution based on
#' the ordered samples, that is, if `T` is a statistic and `X_1`,`X_2`,...,`X_n` 
#' are the samples, the value of `T` is purely based on the value of 
#' `X_(1)`,`X_(2)`,...,`X_(n)`, where `X_(i)` is the ith ordered samples 
#' of `X1`,`X2`,...,`Xn`. Moreover, the rejection region of the 
#' statistic `T` can be written as a set of rejection regions of the 
#' ordered samples `X_(1)`,`X_(2)`,...,`X_(n)`. In other words, 
#' there exist two sequences `{l_i}` and `{u_i}` for `i=1,...,n` and
#' the statistic `T` is in the rejection region if and only if there exist
#' one `i` such that `X_(i) < l_i` or `X_(i) > u_i`.
#' 
#' The most well-known statistic which takes this form is the Kolmogorov-Smirnov 
#' statistic. Other statistics like Berk-Jones or the higher criticism also have
#' the same form but define different sets of `{l_i}` and `{u_i}`. 
#' 
#' \bold{alpha0, index, indexL and indexU}
#' 
#' As mentioned previouly, the rejection of a test can be determined by the 
#' sequences of `{l_i}` and `{u_i}`. Therefore, the parameter `alpha0`, `index`
#' `indexL` and `indexU`. provide a fine control of which `l_i` and `u_i` 
#' will be considered in the test procedure. If no argument is provided, all `l_i`s
#' and `u_i`s will be compared with their corresponding ordered sample `X_(i)`. 
#' This yields the traditional test statistics. If `alpha0` is used, the 
#' test statistics is based on `X_(1),...X_(k)` where k is the nearest integer of 
#' `alpha0*n`. If `index` is provided, only `X_(i)` for `i` in `index` will be considered
#' in a test. If `indexL` and/or `indexU` is not `NULL`, only `l_i` for `i` in 
#' `indexL` and `u_i` for `i` in `indexU` will be used in a test. These can be used to
#' generate an one-sided version of the test statistic.
#' 
#' @examples 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' ## Perform KS test
#' GKSStat(x = x, statName = "KS")
#' 
#' ## Perform one-sided KS test
#' GKSStat(x = x, statName = "KS+")
#' GKSStat(x = x, statName = "KS-")
#' 
#' @return a `generalKSStat` S3 object
#' @rdname statistics
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






#' @rdname statistics
#' @export
HCStat<-function(x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="HC",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = FALSE,levelFunc = HCLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @rdname statistics
#' @export
BJStat<-function(x,alpha0 = 1, index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="BJ",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = TRUE,levelFunc = BJLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @rdname statistics
#' @export
KSStat<-function(x,alpha0=1,index=NULL,indexL=NULL,indexU=NULL){
    statName <- getStatFullName(statName="KS",indexL=indexL,indexU=indexU)
    genericStat(statName=statName, isMin = FALSE,levelFunc = KSLevel,
                x=x,alpha0=alpha0,
                index=index,indexL=indexL,indexU=indexU)
}

#' @rdname statistics
#' @export
HCPlusStat<-function(x,alpha0 = 1, index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    HCStat(x=x,alpha0=alpha0,indexL=index)
}
#' @rdname statistics
#' @export
HCMinusStat<-function(x,alpha0 = 1, index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    HCStat(x=x,alpha0=alpha0,indexU=index)
}

#' @rdname statistics
#' @export
BJPlusStat<-function(x,alpha0 = 1, index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    BJStat(x=x,alpha0=alpha0,indexL=index)
}
#' @rdname statistics
#' @export
BJMinusStat<-function(x,alpha0 = 1, index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    BJStat(x=x,alpha0=alpha0,indexU=index)
}

#' @rdname statistics
#' @export
KSPlusStat <- function(x,alpha0=1,index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    KSStat(x=x,alpha0=alpha0,indexL=index)
}

#' @rdname statistics
#' @export
KSMinusStat <- function(x,alpha0=1,index=NULL){
    index <- getIndexSimple(length(x),alpha0,index)
    KSStat(x=x,alpha0=alpha0,indexU=index)
}
