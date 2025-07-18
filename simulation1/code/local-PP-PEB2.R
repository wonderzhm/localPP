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
