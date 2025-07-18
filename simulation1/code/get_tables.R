#### Tables 3.1, 3.2, 3.3 and A2 in the paper from the intermediate results.

library(xtable)
load(file.path(".", "intermediate_results", "AllMethods_main.RData"))

Table3.1 <- scenarios
rownames(Table3.1) <- paste("S", 1:6, sep = "")
colnames(Table3.1) <- paste("Basket", 1:5, sep = " ")
write.csv(Table3.1, file.path(".", "results", "Table 3.1.csv"))

##### Table 3.2 Overall performance for different methods
method.names <- c("IM", "PP-PEB", "local-PP-PEB1", "local-PP-PEB2", 
                  "PP-GEB", "local-PP-GEB1", "local-PP-GEB2",
                  "JSD1", "JSD2", "EXNEX", "BHM", "BCHM", "local-MEM", "MEM")
nmethods <- length(method.names)
results <- matrix(NA, nmethods, 7)
rownames(results) <- method.names
colnames(results) <- c("Q", "BWER_S1", "BWER-avg", "BWER-max", "TPR-avg", "CCR-avg", "Time (hours)")
ii <- 2
res.oc.all <- list(results.Independent[[ii]],  results.localPP_EBpairwise3[[ii]],
                   results.localPP_EBpairwise[[ii]], results.localPP_EBpairwise2[[ii]],
                   results.localPP_EBglobal3[[ii]], results.localPP_EBglobal[[ii]], 
                   results.localPP_EBglobal2[[ii]], results.JSD[[ii]],
                   results.JSD2[[ii]], results.EXNEX[[ii]], results.BHMunif[[ii]], 
                   results.BCHM[[ii]], results.localMEM[[ii]], results.MEM[[ii]])
time.all <- c(time.Independent, time.localPP_EBpairwise3, time.localPP_EBpairwise, 
              time.localPP_EBpairwise2, time.localPP_EBglobal3, time.localPP_EBglobal,
              time.localPP_EBglobal2, time.JSD, time.JSD2, time.EXNEX, time.BHMunif, time.BCHM, 
              time.localMEM, time.MEM)/3600
for(i in 1:nmethods){
  res <- res.oc.all[[i]]
  results[i, ]<- sprintf(c(res$Q[1], res$error.tw, mean(res$bwer, na.rm = TRUE), max(res$bwer, na.rm = TRUE),
                           res$power.cdr, res$power.ccr, time.all[i]), 
                         fmt = '%#.3f')
}
#results
indx.order <- c(1,2,3,5,6,8,10:14, 4,7,9)
print(xtable(results[indx.order,]))
write.csv(results[indx.order,], file.path(".", "results", "Table 3.2.csv"))

## Table 3.3  Results by scenarios
method.indx <- c(1,4,7,9,10); method.names[method.indx]
method.names_sub <- c("IM", "local-PP-PEB", "local-PP-GEB",
                      "JSD", "EXNEX")
nmethods_sub <- length(method.names_sub)
res.allmethods <- res.oc.all[method.indx]
oc.allmethods <- list(oc.Independent, oc.localPP_EBpairwise2, oc.localPP_EBglobal2,
                      oc.JSD2, oc.EXNEX)
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
write.csv(res.by.sc, file.path(".", "results", "Table 3.3.csv"))

## Table A2  Results by scenarios for remaining methods
method.indx <- c(2,3,5,6,8,11,12,13,14); method.names[method.indx]
method.names_sub <- c("PP-PEB", "local-PP-PEB",
                      "PP-GEB", "local-PP-GEB","JSD", 
                      "BHM", "BCHM", "local-MEM", "MEM")
nmethods_sub <- length(method.names_sub)
res.allmethods <- res.oc.all[method.indx]
oc.allmethods <- list(oc.localPP_EBpairwise3, oc.localPP_EBpairwise, 
                      oc.localPP_EBglobal3, oc.localPP_EBglobal, 
                      oc.JSD, oc.BHMunif, oc.BCHM, oc.localMEM, oc.MEM)
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
write.csv(res.by.sc, file.path(".", "results", "Table A2.csv"))
