# Data Analysis for BRAF-V600 study
# Generates Tables 4.1, 4.2 and A3. 

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
source(file.path("..", "functions", "functions.R"))
source(file.path("..", "functions", "functions_parallel.R"))

### BRAF V600 trial data and design settings
N <- rbind(19, 10, 26, 8, 14, 7) #  total sample size for each indication
responses <- rbind(8, 0, 1, 1, 6, 2) # number of responses for each indication
nperclust <- 10000 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 100000 MC replicates
B <- nrow(N) # total number of baskets
p0 <- 0.15 # null response rate
p1 <- 0.35 # target response rate
pnull <- rep(p0, B)  
ptarget <- rep(p1, B) 
BFs <- rep(NA, B)
for(i in 1:B) BFs[i] <- sum(N[-i,ncol(N)])/N[i,ncol(N)]
a.max <- max(BFs) ## maximum borrowing factor

## BOP2 for each indication with error rate sig.level=0.05
sig.level <- 0.05 # type I error
stopbounds <- NULL # since no interim
beta.a0 <- pnull # default beta prior 
beta.b0 <- 1-pnull # default beta prior 
ndigits = 3 ## number of digits for Q

## scenarios 
scenarios <- rbind(pnull) # global null

## generate data under global null to calibrate Q
seed <- 2024
data.object<-generate.data(N, scenarios, ntrial = nperclust*nclust, seed = seed)

####################### Independent ########################
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = beta.a0, 
                                beta.b0 = beta.b0, seed = seed, ModelFit = "Independent")
Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits)
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.Independent <- end_time - start_time)
## Q
(Q.Independent <- Q)
## type I error
results.Independent <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.Independent <- apply(res.post$postprob>Qmat, c(1,3), mean))
## posterior probabilities
fit <- Independent(nDat = c(19, 10, 26, 8, 14, 7), 
                   yDat = c(8, 0, 1, 1, 6, 2), 
                   be.a0 = beta.a0, be.b0 = beta.b0)
pp.Independent <- pbeta(0.15, shape1 = fit$a.post, 
                        shape2 = fit$b.post, lower.tail = FALSE)
pp.Independent

####################### local PP ########################
a = 1
delta = 0.4
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = beta.a0, 
                                beta.b0 = beta.b0, seed = seed, ModelFit = "localPP", 
                                a = a, delta=delta, method = "PEB")
Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits)
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.localPP <- end_time - start_time)
## Q
(Q.localPP <- Q)
## type I error
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.localPP <- apply(res.post$postprob>Qmat, c(1,3), mean))
## posterior probabilities
fit <- localPP(nDat = c(19, 10, 26, 8, 14, 7), yDat = c(8, 0, 1, 1, 6, 2),
               be.a0 = beta.a0, be.b0 = beta.b0, a = a, delta = delta, method = "PEB")
pp.localPP <- pbeta(0.15, shape1 = fit$a.post, 
                    shape2 = fit$b.post, lower.tail = FALSE)
pp.localPP

save.image(file.path(".", "intermediate_results", "case_study.RData"))

####################################################################################
## Steps below read the case_study.RData file and 
## create Tables 4.1, 4.2 and A3 in the paper
####################################################################################
source(file.path(".", "code", "get_tables.R"))
