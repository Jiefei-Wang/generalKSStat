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


#' @export
print.generalKSStat <- function(x,...){
    #print(x$statValue)
    class(x)=NULL
    print(x)
    invisible(x)
}
#' @export
as.numeric.generalKSStat <- function(x,...){
    x[["statValue"]]
}