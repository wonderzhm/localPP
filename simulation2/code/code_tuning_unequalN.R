# This R code is used to select the tuning parameters for 
# local-PP-PEB, local-PP-GEB and JSD in Section 3.4 of main paper. 
# NOTE: the grid points of tuning parameters were reduced to save the running time
#       for real applications, please consider using the full range suggested in the paper. 

# Load required packages
library(parallel)
library(doParallel)
library(R2jags)
library(foreach)
library(bhmbasket) 
library(basket) 
library(partitions)
library(plyr)
library(cluster)
library(coda)
library(rjags)
library(BasketTrial)
source("../functions/functions.R")
source("../functions/functions_parallel.R")

nperclust <- 500 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 5000 MC replicates
N <- rbind(c(10, 26),
           c(10, 16),
           c(4,   8), # the first number does not matter since no interim for this basket
           c(10, 17),
           c(10, 22)
) # interim sample size and total sample size for each indication
B <- nrow(N) # total number of baskets
B <- nrow(N) # total number of baskets
BFs <- rep(NA, B)
for(i in 1:B) BFs[i] <- sum(N[-i,ncol(N)])/N[i,ncol(N)]
a.max <- max(BFs) ## maximum borrowing factor
p0 <- 0.15
p1 <- 0.30
p2 <- 0.45
pnull <- c(p0, p0, p0, p0, p0) # null response rate for each indication
ptarget <- c(p1, p1, p1, p1, p1) # target response rate for each indication

## BOP2 for each indication with error rate sig.level=0.1
sig.level <- 0.1 # type I error
stopbounds <- rbind(1, 
                    1, 
                    -1, # no interim for this basket due to the small maximum sample size
                    1, 
                    1) # obtained from BOP2 app
beta.a0 <- pnull # default beta prior 
beta.b0 <- 1-pnull # default beta prior 
ndigits = 3 ## number of digits for Q

## scenarios 
scenarios <- rbind( c(p0, p0, p0, p0, p0),
                    c(p0, p0, p0, p1, p1),
                    c(p0, p1, p1, p1, p1),
                    c(p0, p1, p1, p2, p2),
                    c(p0, p2, p2, p2, p2),
                    c(p1, p1, p1, p1, p1)
)

## generate data which will be analyzed by all methods
seed <- 2024
data.object<-generate.data(N, scenarios, ntrial = nperclust*nclust, seed = seed)

####################################################################################
## Selecting hyperparameter a for local PP with pairwise EB weights
####################################################################################
# ranges of a and delta have been reduced to save time; please use the following in practice
# as a starting point and then refine it using finer grid points by fixing delta at one or two values.
# aseq <- seq(0.1, a.max, length.out = 10)
# deltaseq <- seq(0.1, 0.4, 0.1)
aseq <- c(seq(0.5, 0.6, 0.05)) 
deltaseq <- seq(0.4, 0.4, 0.1) 
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
results <- matrix(NA, length(aseq)*length(deltaseq), 7)
colnames(results) <- c("a","delta","BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                    a = aseq[i], delta = deltaseq[j], method = "PEB")
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits))
    res <- get.weighted.power(res.post, Q = Q0)
    results[indx, ]<- sprintf(c(aseq[i], deltaseq[j], res$error.tw, 
                                mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                                res$power.cdr, res$power.ccr), fmt = '%#.3f')
  }
}
(res.PEB <- results[order(as.numeric(results[,5])),])

####################################################################################
## Selecting hyperparameter a for local PP with global EB weights
####################################################################################
# ranges of a and delta have been reduced to save time; please use the following in practice
# as a starting point and then refine it using finer grid points by fixing delta at one or two values.
# aseq <- seq(0.1, a.max, length.out = 10)
# deltaseq <- seq(0.1, 0.4, 0.1)
aseq <- c(seq(0.5, 0.6, 0.05)) 
deltaseq <- seq(0.4, 0.4, 0.1) 
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
results <- matrix(NA, length(aseq)*length(deltaseq), 7)
colnames(results) <- c("a","delta","BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                    a = aseq[i], delta = deltaseq[j], method = "GEB")
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits))
    res <- get.weighted.power(res.post, Q = Q0)
    results[indx, ]<- sprintf(c(aseq[i], deltaseq[j], res$error.tw, 
                                mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                                res$power.cdr, res$power.ccr), fmt = '%#.3f')
  }
}
(res.GEB <- results[order(as.numeric(results[,5])),])

####################################################################################
## Selecting hyperparameter tau for JSD
####################################################################################
# ranges of a and delta have been reduced to save time; please use the following in practice
# as a starting point and then refine it using finer grid points by fixing delta at one or two values.
# aseq <- seq(1, 7, 0.5)
# deltaseq <- seq(0, 1, 0.1)
aseq <- c(seq(6, 7, 0.5))  ## values for epsilon
deltaseq <- seq(0.5, 0.5, 0.1) ## values for tau
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
results <- matrix(NA, length(aseq)*length(deltaseq), 7)
colnames(results) <- c("a","delta","BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "ccr-avg")
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", 
                                    epsilon = aseq[i], tau = deltaseq[j])
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits))
    res <- get.weighted.power(res.post, Q = Q0)
    results[indx, ]<- sprintf(c(aseq[i], deltaseq[j], res$error.tw, 
                                mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                                res$power.cdr, res$power.ccr), fmt = '%#.3f')
  }
}
(res.JSD <- results[order(as.numeric(results[,5])),])

