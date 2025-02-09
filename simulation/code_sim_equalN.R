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
## Simulation results for comparing the  following eight methods:
## Independent, local PP, JSD, EXNEX, BHM, BCHM, local MEM and MEM
####################################################################################

####################### Independent ########################
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "Independent")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.Independent <- end_time - start_time)
## save results
results.Independent <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.Independent <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP pairwise EB weights ########################
## prior values
a = 0.9
delta = 0.4
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "PEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBpairwise <- end_time - start_time)
## save results
results.localPP_EBpairwise <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBpairwise <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP pairwise EB weights matching EXNEX type I error ########################
## prior values
a = 0.35
delta = 0.40
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "PEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBpairwise2 <- end_time - start_time)
## save results
results.localPP_EBpairwise2 <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBpairwise2 <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP pairwise EB weights unadjusted ########################
## prior values
a = 4
delta = 1
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "PEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBpairwise3 <- end_time - start_time)
## save results
results.localPP_EBpairwise3 <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBpairwise3 <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP global EB weights ########################
## prior values
a = 3
delta = 0.4
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "GEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBglobal <- end_time - start_time)
## save results
results.localPP_EBglobal <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBglobal <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP global EB weights matching exnex type I error ########################
## prior values
a = 0.45
delta = 0.4
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "GEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBglobal2 <- end_time - start_time)
## save results
results.localPP_EBglobal2 <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBglobal2 <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local PP global EB weights unadjusted ########################
## prior values
a = 4
delta = 1
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "GEB")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP_EBglobal3 <- end_time - start_time)
## save results
results.localPP_EBglobal3 <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP_EBglobal3 <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### JSD ########################
## prior values
epsilon = 3
tau = 0.5
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", 
                                tau = tau, epsilon = epsilon)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.JSD <- end_time - start_time)
## save results
results.JSD <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.JSD <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### JSD matching exnex type I error ########################
## prior values
epsilon = 6.5
tau = 0.5
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "JSD", 
                                tau = tau, epsilon = epsilon)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.JSD2 <- end_time - start_time)
## save results
results.JSD2 <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.JSD2 <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### local MEM ########################
## prior values
delta = 2 
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "localMEM", delta = delta)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localMEM <- end_time - start_time)
## save results
results.localMEM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localMEM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### BHM uniform ########################
## prior values
u0 = 100
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "BHMunif", u0 = u0)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.BHMunif <- end_time - start_time)
## save results
results.BHMunif <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.BHMunif <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### EXNEX ########################
## prior values
nex.w = 0.5
weight = c((1-nex.w)/2, (1-nex.w)/2, nex.w)
pguess = c(pnull[1], ptarget[1], (pnull[1] + ptarget[1])/2)
tau.HN.scale = rep(1, length(weight)-1)
Nmix <- length(weight)
Nexch <- Nmix - 1
mu.mean <- rep(NA, Nexch)
mu.prec <- rep(NA, Nexch)
for (i in 1:Nexch) {
  prior <- getPriorParameters("exnex", target_rates = pguess[i], tau_scale = tau.HN.scale[i])
  mu.mean[i] <- prior$exnex$mu_mean
  mu.prec[i] <- 1 / prior$exnex$mu_sd^2
}
prior <- getPriorParameters("exnex", target_rates = pguess[Nmix])
nex.mean <- prior$exnex$mu_j
nex.prec <- 1 / prior$exnex$tau_j^2

## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "EXNEX", 
                                weight = weight, nex.mean = nex.mean, nex.prec = nex.prec, 
                                mu.mean = mu.mean, mu.prec = mu.prec, tau.HN.scale = tau.HN.scale)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.EXNEX <- end_time - start_time)
## save results
results.EXNEX <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.EXNEX <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### BCHM ########################
## prior values
alpha = 1e-40 # precision for DP clustering
alpha1 = 50 # gamma(alpha1, beta1) for 1/variance
beta1 = 10 # gamma(alpha1, beta1) for 1/variance
mu0 = logit(pnull[1]) # hyperprior mean
tau2 = 0.01 # hyperprior Precision 
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "BCHM", 
                                alpha = alpha, alpha1 = alpha1, beta1 = beta1, mu0 = mu0, tau2 = tau2)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.BCHM <- end_time - start_time)
## save results
results.BCHM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.BCHM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

####################### MEM ########################
## prior values
shape1 <- pnull[1]
shape2 <- 1- pnull[1]
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "MEM", 
                                shape1 = shape1, shape2 = shape2)
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits, Qclust = rep(1, B)))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.MEM <- end_time - start_time)
## save results
results.MEM <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.MEM <- apply(res.post$postprob>Qmat, c(1,3), mean))
save.image("intermediate_results/AllMethods_main.RData")

