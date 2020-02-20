context("Test statistics")
x <- 1:9/10
index <- c(2,3)

test_that("KS",{
    stat <- GKSStat(x=x,index=index,statName = "KS",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.089)
    expect_equal(round(stat$pvalue,3),0.946)
})
test_that("KS+",{
    stat <- GKSStat(x=x,index=index,statName = "KS+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.033)
    expect_equal(round(stat$pvalue,3),0.643)
})
test_that("KS-",{
    stat <- GKSStat(x=x,index=index,statName = "KS-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.089)
    expect_equal(round(stat$pvalue,3),0.542)
})


test_that("BJ",{
    stat <- GKSStat(x=x,index=index,statName = "BJ",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.436)
    expect_equal(round(stat$pvalue,3),0.978)
})

test_that("BJ+",{
    stat <- GKSStat(x=x,index=index,statName = "BJ+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.537)
    expect_equal(round(stat$pvalue,3),0.646)
})

test_that("BJ-",{
    stat <- GKSStat(x=x,index=index,statName = "BJ-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.436)
    expect_equal(round(stat$pvalue,3),0.543)
})

test_that("HC",{
    stat <- GKSStat(x=x,index=index,statName = "HC",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.667)
    expect_equal(round(stat$pvalue,3),0.940)
})
test_that("HC+",{
    stat <- GKSStat(x=x,index=index,statName = "HC+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.218)
    expect_equal(round(stat$pvalue,3),0.648)
})

test_that("HC-",{
    stat <- GKSStat(x=x,index=index,statName = "HC-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.667)
    expect_equal(round(stat$pvalue,3),0.527)
})
