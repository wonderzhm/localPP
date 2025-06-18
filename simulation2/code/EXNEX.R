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
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.EXNEX <- end_time - start_time)
## save results
results.EXNEX <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.EXNEX <- apply(res.post$postprob>Qmat, c(1,3), mean))
