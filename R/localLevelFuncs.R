## Indivial level function
HCPlusLevel <- function(n, x, sx, index){
  if (length(index) == 0)
    return(numeric(0))
  sqrt(n) * (seq(1, n)[index] / n - sx[index]) / sqrt(sx[index] * (1 - sx[index]))
}
HCMinusLevel<- function(n, x, sx, index){
  if (length(index) == 0)
    return(numeric(0))
  sqrt(n) * (sx[index]-(seq(1, n)[index]-1) / n) / sqrt(sx[index] * (1 - sx[index]))
}
HCLevel <- function(n, x, sx, index, indexU, indexL) {
  HCPlus <- HCPlusLevel(n, x, sx,indexL)
  HCMinus <- HCMinusLevel(n, x, sx,indexU)
  c(HCPlus, HCMinus)
}
BJPlusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  sapply(seq_along(x)[index], function(x)
    pbeta(sx[x], x, n - x + 1))
}
BJMinusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  1 - BJPlusLevel(n, x, sx, index)
}
BJLevel <- function(n, x, sx, index, indexL, indexU) {
  BJPlus <- BJPlusLevel(n, x, sx,indexL)
  BJMinus <- BJMinusLevel(n, x, sx,indexU)
  c(BJPlus, BJMinus)
}
KSPlusLevel <- function(n, x, sx, index) {
  seq(1, n)[index] / n - sx[index]
}
KSMinusLevel <- function(n, x, sx, index) {
  sx[index] - (seq(1, n)[index] - 1) / n
}
KSLevel <- function(n, x, sx, index, indexL, indexU) {
  KSPlus <- KSPlusLevel(n, x, sx,indexL)
  KSMinus <- KSMinusLevel(n, x, sx,indexU)
  c(KSPlus, KSMinus)
}

## These functions return a set of level stat
partialLevelStat <- function(statFunc, x, alpha0, index, indexL, indexU) {
  n <- length(x)
  sx <- sort(x)
  sx[sx == 0] <- min(10 ^ -6, sx[sx != 0])
  sx[sx == 1] <- max(1 - 10 ^ -6, sx[sx != 1])
  rangeIndex <- getIndex(n=n,alpha0=alpha0, 
                         index=index, indexL=indexL, indexU=indexU)
  statFunc(
    n = n,
    x = x,
    sx = sx,
    index = rangeIndex$index,
    indexU = rangeIndex$indexU ,
    indexL = rangeIndex$indexL
  )
}





## get local critical value

HCLocalCritical<-function(stat,n){
  stat <- stat/sqrt(n)
  a<-1+stat^2
  ## lower
  const<-seq_len(n)/n
  b<--2*const-stat^2
  l <- (-b-sqrt(b^2-4*a*const^2))/2/a
  ## upper
  const<-(seq_len(n)-1)/n
  b<--2*const-stat^2
  h <- (-b+sqrt(b^2-4*a*const^2))/2/a
  list(l =l,h= h)
}

BJLocalCritical<-function(stat,n){
  l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
  h=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))
  list(l =l,h= h)
}

KSLocalCritical<-function(stat,n){
  l <- 1:n/n - stat
  h <- stat + 1:n/n-1/n
  l[l<0]=0
  h[h>1]=1
  list(l =l,h= h)
}

