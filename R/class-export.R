.generalKSStat<-function(statName,statValue,n,alpha0=NULL,
                     index=NULL,indexL=NULL,indexU=NULL){
    stat <- list(statName = statName,
                 statValue = statValue,
                 n=n)
    if(all.null(index,indexL,indexU)){
    stat[["alpha0"]]=alpha0
    }else{
        stat[["index"]]=index
        stat[["indexL"]]=indexL
        stat[["indexU"]]=indexU
    }
    structure(stat,class = "generalKSStat")
}

#' print function for the generalKSStat S3 object
#' 
#' @param x the generalKSStat S3 object
#' @examples 
#' 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' ## Perform KS test
#' GKSStat(x = x, statName = "KS")
#' 
#' @return invisible x
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