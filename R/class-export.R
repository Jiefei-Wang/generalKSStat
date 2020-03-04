.generalKSStat<-function(statName,statValue,n,indexL=NULL,indexU=NULL){
    stat <- list(statName = statName,
                 statValue = statValue,
                 n=n)
    stat[["indexL"]]=indexL
    stat[["indexU"]]=indexU
    structure(stat,class = "generalKSStat")
}

#' Class methods for the generalKSStat S3 object
#' 
#' 
#' @param x the generalKSStat S3 object
#' @param ... Ignored.
#' @examples 
#' 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' ## Perform KS test
#' GKSStat(x = x, statName = "KS")
#' 
#' @return  
#' print: invisible `x`
#' other fucntions: a numeric value
#' 
#' @rdname classMethod
#' @export
print.generalKSStat <- function(x,...){
    #print(x$statValue)
    # class(x)=NULL
    # print(x)
    cat("The", x$statName, "test statistics\n")
    cat("Sample size:",x$n,"\n")
    cat("Stat value:", x$statValue,"\n")
    if(!is.null(x$pvalue))
        cat("P-value:", x$pvalue)
    invisible(x)
}

#' @rdname classMethod
#' @export
getStatName<-function(x){
    x$statName
}
#' @rdname classMethod
#' @export
getStatValue<-function(x){
    x$statValue
}
#' @rdname classMethod
#' @export
getPvalue<-function(x){
    x$pvalue
}
