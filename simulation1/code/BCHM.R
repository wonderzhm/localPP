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
