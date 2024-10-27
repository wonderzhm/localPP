library(xtable)
####################################################################################
## Create tables for appendix
####################################################################################
load("intermediate_results/AllMethods_supp.RData")

##### Table A3 The efficacy cutoff Q_is for all considered methods.
method.names <- c("IM", "local-PP-PEB", "local-PP-GEB", "JSD",  "EXNEX")
nmethods <- length(method.names)
ii <- 2; 
res.post.all <- list(results.Independent[[ii]], results.localPP_EBpairwise[[ii]],
                     results.localPP_EBglobal[[ii]], results.JSD[[ii]],
                     results.EXNEX[[ii]])
Q.by.method <- matrix(NA, nmethods, B)
rownames(Q.by.method) <- method.names
colnames(Q.by.method) <- c("Basket1", "Basket2", "Basket3", "Basket4", "Basket5")
for(i in 1:nmethods){
  res <- res.post.all[[i]]
  Q.by.method[i,] <- sprintf(res$Q, fmt = '%#.3f')
}
print(xtable(Q.by.method))

##### Table A4 Overall performance for different methods
results <- matrix(NA, nmethods, 5)
rownames(results) <- method.names
colnames(results) <- c("BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "CCR-avg")
ii <- 2
res.oc.all <- list(results.Independent[[ii]], results.localPP_EBpairwise[[ii]],
                   results.localPP_EBglobal[[ii]], results.JSD[[ii]],
                   results.EXNEX[[ii]])
for(i in 1:nmethods){
  res <- res.oc.all[[i]]
  results[i, ]<- sprintf(c(res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr), 
                         fmt = '%#.3f')
}
#results
print(xtable(results))

## Table A5  Results by scenarios
method.indx <- 1:5; method.names[method.indx]
method.names_sub <- c("IM", "local-PP-PEB", "local-PP-GEB", "JSD", "EXNEX")
nmethods_sub <- length(method.names_sub)
res.allmethods <- res.oc.all[method.indx]
oc.allmethods <- list(oc.Independent, oc.localPP_EBpairwise, oc.localPP_EBglobal,
                      oc.JSD, oc.EXNEX)
oc.by.sc <- matrix(NA, nmethods_sub, 9)

rownames(oc.by.sc) <- method.names_sub
colnames(oc.by.sc) <- c("Basket1", "Basket2", "Basket3", "Basket4", "Basket5", "FPR", "FDR", "TPR", "CCR")
res.by.sc <- list(scenario1=NULL, scenario2=NULL, scenario3=NULL,
                  scenario4=NULL, scenario5=NULL, scenario6=NULL)
for(k in 1:nrow(scenarios)){
  for(i in 1:nmethods_sub){
    res <- res.allmethods[[i]]
    oc.i <- oc.allmethods[[i]]
    oc.by.sc[i,1:5]<- sprintf(oc.i[k,], fmt = '%#.3f')
    oc.by.sc[i,6:9]<- sprintf(c(res$ind.error.tw[k], res$ind.error.fdr[k], 
                                res$ind.power.cdr[k], res$ind.power.ccr[k]), fmt = '%#.3f')
  }
  res.by.sc[[k]] <- oc.by.sc
  print(paste("Scenario", k, sep=""))
  print(xtable(oc.by.sc))
  cat("\n")
}
#res.by.sc
