 qfun.basic <- function(gamma)
{
        qarr <- array(0, dim = c(2, 4, 4))
        theta <- exp(gamma)
        qarr[1, 1, 1] <-  - theta[1]
        qarr[1, 1, 2] <- theta[1]
        qarr[1, 2, 1] <- theta[2]
        qarr[1, 2, 2] <-  - theta[2] - theta[3] - theta[6]
        qarr[1, 2, 3] <- theta[3]
        qarr[1, 2, 4] <- theta[6]
        qarr[1, 3, 2] <- theta[4]
        qarr[1, 3, 3] <-  - theta[4] - theta[5]
        qarr[1, 3, 4] <- theta[5]
        qarr[2, 1, 1] <-  - theta[7]
        qarr[2, 1, 2] <- theta[7]
        qarr[2, 2, 1] <- theta[8]
        qarr[2, 2, 2] <-  - theta[8] - theta[9] - theta[12]
        qarr[2, 2, 3] <- theta[9]
        qarr[2, 2, 4] <- theta[12]
        qarr[2, 3, 2] <- theta[10]
        qarr[2, 3, 3] <-  - theta[10] - theta[11]
        qarr[2, 3, 4] <- theta[11]
        return(qarr)
}
qderivf <- function(gamma)
{
        rmat <- array(0, c(12, 2, 4, 4))
        theta <- exp(gamma)
        rmat[1, 1, 1, 1] <- ( - theta[1])
        rmat[1, 1, 1, 2] <- theta[1]
        rmat[7, 2, 1, 1] <- ( - theta[7])
        rmat[7, 2, 1, 2] <- theta[7]
        rmat[2, 1, 2, 1] <- theta[2]
        rmat[2, 1, 2, 2] <- ( - theta[2])
        rmat[8, 2, 2, 1] <- theta[8]
        rmat[8, 2, 2, 2] <- ( - theta[8])
        rmat[3, 1, 2, 2] <- ( - theta[3])
        rmat[3, 1, 2, 3] <- theta[3]
        rmat[9, 2, 2, 2] <- ( - theta[9])
        rmat[9, 2, 2, 3] <- theta[9]
        rmat[4, 1, 3, 2] <- theta[4]
        rmat[4, 1, 3, 3] <- ( - theta[4])
        rmat[10, 2, 3, 2] <- theta[10]
        rmat[10, 2, 3, 3] <- ( - theta[10])
        rmat[5, 1, 3, 3] <- ( - theta[5])
        rmat[5, 1, 3, 4] <- theta[5]
        rmat[11, 2, 3, 3] <- ( - theta[11])
        rmat[11, 2, 3, 4] <- theta[11]
        rmat[6, 1, 2, 4] <- theta[6]
        rmat[6, 1, 2, 2] <- ( - theta[6])
        rmat[12, 2, 2, 4] <- theta[12]
        rmat[12, 2, 2, 2] <- ( - theta[12])
        return(rmat)
}
qfun.m2 <- function(gamma)
{
        qarr <- array(0, dim = c(2, 4, 4))
        theta <- exp(gamma)
        qarr[1, 1, 1] <-  - theta[1]
        qarr[1, 1, 2] <- theta[1]
        qarr[1, 2, 1] <- theta[2]
        qarr[1, 2, 2] <-  - theta[2] - theta[3] - theta[4]
        qarr[1, 2, 3] <- theta[3]
        qarr[1, 2, 4] <- theta[4]
        qarr[1, 3, 2] <- theta[5]
        qarr[1, 3, 3] <-  - theta[5] - theta[6]
        qarr[1, 3, 4] <- theta[6]
        qarr[2, 1, 1] <-  - theta[7]
        qarr[2, 1, 2] <- theta[7]
        qarr[2, 2, 1] <- theta[2]
        qarr[2, 2, 2] <-  - theta[2] - theta[8] - theta[10]
        qarr[2, 2, 3] <- theta[8]
        qarr[2, 2, 4] <- theta[10]
        qarr[2, 3, 2] <- theta[4]
        qarr[2, 3, 3] <-  - theta[4] - theta[9]
        qarr[2, 3, 4] <- theta[9]
        return(qarr)
}
qderiv.m2 <- function(gamma)
{
        rmat <- array(0, c(10, 2, 4, 4))
        theta <- exp(gamma)
        rmat[1, 1, 1, 1] <- ( - theta[1])
        rmat[1, 1, 1, 2] <- theta[1]
        rmat[2, 1, 2, 1] <- theta[2]
        rmat[2, 1, 2, 2] <- ( - theta[2])
        rmat[2, 2, 2, 1] <- theta[2]
        rmat[2, 2, 2, 2] <- ( - theta[2])
        rmat[3, 1, 2, 2] <- ( - theta[3])
        rmat[3, 1, 2, 3] <- theta[3]
        rmat[4, 1, 3, 2] <- theta[4]
        rmat[4, 1, 3, 3] <- ( - theta[4])
        rmat[4, 2, 3, 2] <- theta[4]
        rmat[4, 2, 3, 3] <- ( - theta[4])
        rmat[5, 1, 3, 3] <- ( - theta[5])
        rmat[5, 1, 3, 4] <- theta[5]
        rmat[6, 1, 2, 4] <- theta[6]
        rmat[6, 1, 2, 2] <- ( - theta[6])
        rmat[7, 2, 1, 1] <- ( - theta[7])
        rmat[7, 2, 1, 2] <- theta[7]
        rmat[8, 2, 2, 2] <- ( - theta[8])
        rmat[8, 2, 2, 3] <- theta[8]
        rmat[9, 2, 3, 3] <- ( - theta[9])
        rmat[9, 2, 3, 4] <- theta[9]
        rmat[10, 2, 2, 4] <- theta[10]
        rmat[10, 2, 2, 2] <- ( - theta[10])
        return(rmat)
}
