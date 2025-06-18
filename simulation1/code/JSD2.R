####################### JSD2 ########################
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
