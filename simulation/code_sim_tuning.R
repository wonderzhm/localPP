# Replicate the simulation results

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
library(BasketTrial) # install it by devtools::install_github("wonderzhm/BasketTrial")
source("../functions/functions.R")
source("../functions/functions_parallel.R")

nperclust <- 500 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 5000 MC replicates
N <- rbind(c(10, 25),
           c(10, 25),
           c(10, 25),
           c(10, 25),
           c(10, 25)
) # interim sample size and total sample size for each indication
B <- nrow(N) # total number of baskets
p0 <- 0.15
p1 <- 0.30
p2 <- 0.45
pnull <- c(p0, p0, p0, p0, p0) # null response rate for each indication
ptarget <- c(p1, p1, p1, p1, p1) # target response rate for each indication

## BOP2 for each indication with error rate sig.level=0.1
sig.level <- 0.1 # type I error
stopbounds <- rbind(1, 
                    1, 
                    1, 
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
BFs <- rep(NA, B)
for(i in 1:B) BFs[i] <- sum(N[-i,ncol(N)])/N[i,ncol(N)]
a.max <- max(BFs) ## maximum borrowing factor
aseq = c(seq(0.25, 0.45, 0.05), seq(0.8, 1, 0.05)) # show only a subset
deltaseq <- seq(0.2, 0.4, 0.1)
list.metrics0 <- list()
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                    a = aseq[i], delta = deltaseq[j], method = "PEB")
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
    res0 <- get.weighted.power(res.post, Q = Q0)
    list.metrics0[[indx]] <- res0
  }
}
save.image("intermediate_results/localPP_EBpaiwise_a_main.RData")

####################################################################################
## Selecting hyperparameter a for local PP with global EB weights
####################################################################################
aseq = c(seq(0.35, 0.55, 0.05), seq(2.9, 3.1, 0.05)) # show only a subset
deltaseq <- seq(0.2, 0.4, 0.1)
list.metrics0 <- list()
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                    a = aseq[i], delta = deltaseq[j], method = "GEB")
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
    res0 <- get.weighted.power(res.post, Q = Q0)
    list.metrics0[[indx]] <- res0
  }
}
save.image("intermediate_results/localPP_EBglobal_a_main.RData")

####################################################################################
## Selecting hyperparameter tau for JSD
####################################################################################
aseq <- seq(3, 7, 0.5) ## values for epsilon
deltaseq = seq(0.1, 1, 0.1) ## values for tau 
list.metrics0 <- list()
a.delta <- matrix(NA, nrow = length(aseq)*length(deltaseq), 2)
for(i in 1:length(aseq)){
  for(j in 1:length(deltaseq)){
    indx <- j+(i-1)*length(deltaseq)
    a.delta[indx,] <- c(aseq[i], deltaseq[j])
    print(c(aseq[i], deltaseq[j]))
    res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                    pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                    beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", 
                                    epsilon = aseq[i], tau = deltaseq[j])
    (Q0 <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
    res0 <- get.weighted.power(res.post, Q = Q0)
    list.metrics0[[indx]] <- res0
  }
}
save.image("intermediate_results/JSD_select_tau_main.RData")
