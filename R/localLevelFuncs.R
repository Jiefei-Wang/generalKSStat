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
