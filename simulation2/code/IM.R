####################### Independent ########################
## start simulation
seed <- 2024
start_time <- Sys.time()
res.post <- post.infer.parallel(nclust = nclust, nperclust = nperclust, data.object = data.object, 
                                pnull = pnull, stopbounds =  stopbounds, beta.a0 = pnull, 
                                beta.b0 = 1-pnull, seed = seed, ModelFit = "Independent")
(Q <- get.Q.bwer(res.post, alpha = sig.level, digits = ndigits))
res <- get.weighted.power(res.post, Q = Q)
end_time <- Sys.time()
(time.Independent <- end_time - start_time)
## save results
results.Independent <- list(res.post, res)
Qmat <- array(NA, dim = dim(res.post$postprob))
for(i in 1:B) Qmat[,,i] <- res$Q[i]
(oc.Independent <- apply(res.post$postprob>Qmat, c(1,3), mean))
