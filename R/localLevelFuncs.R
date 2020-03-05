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
HCLevel <- function(n, x, sx, indexU, indexL) {
  HCPlus <- max(HCPlusLevel(n, x, sx,indexL),0)
  HCMinus <- max(HCMinusLevel(n, x, sx,indexU),0)
  max(c(HCPlus, HCMinus))
}
BJPlusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  vapply(seq_along(x)[index], function(x)
    pbeta(sx[x], x, n - x + 1),numeric(1))
}
BJMinusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  1 - BJPlusLevel(n, x, sx, index)
}
BJLevel <- function(n, x, sx, indexL, indexU) {
  BJPlus <- min(BJPlusLevel(n, x, sx,indexL),1)
  BJMinus <- min(BJMinusLevel(n, x, sx,indexU),1)
  min(BJPlus, BJMinus)
}
KSPlusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  index / n - sx[index]
}
KSMinusLevel <- function(n, x, sx, index) {
  if (length(index) == 0)
    return(numeric(0))
  sx[index] - (index - 1) / n
}
KSLevel <- function(n, x, sx, indexL, indexU) {
  KSPlus <- max(KSPlusLevel(n, x, sx,indexL),0)
  KSMinus <- max(KSMinusLevel(n, x, sx,indexU),0)
  max(KSPlus, KSMinus)
}

## These functions return a set of level stat
partialLevelStat <- function(statFunc, x, indexL, indexU) {
  n <- length(x)
  sx <- sort(x)
  sx[sx == 0] <- min(10 ^ -6, sx[sx != 0])
  sx[sx == 1] <- max(1 - 10 ^ -6, sx[sx != 1])
  statFunc(
    n = n,
    x = x,
    sx = sx,
    indexU = indexU ,
    indexL = indexL
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
  l=vapply(seq_len(n),function(x)qbeta(stat,x,n-x+1),numeric(1))
  h=vapply(seq_len(n),function(x)qbeta(1 - stat,x,n-x+1),numeric(1))
  list(l =l,h= h)
}

KSLocalCritical<-function(stat,n){
  l <- seq_len(n)/n - stat
  h <- stat + seq_len(n)/n-1/n
  l[l<0]=0
  h[h>1]=1
  list(l =l,h= h)
}

