# Load required packages
library(BasketTrial) # install it by devtools::install_github("wonderzhm/BasketTrial")
library(xtable)

### Similarity Matrix: element (i,j) is the power parameter between baskets i and j
### pairwise empirical Bayes (PEB) raw weight s_ij: i.e. set a at any value >=4
res <- localPP(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20),
               be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5), a = 4, delta = 1, 
               method = "PEB")
print(xtable(round(res$sm,3)))

### global empirical Bayes (GEB) raw weight s_ij: i.e. set a at any value >=4
res <- localPP(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20),
               be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5), a = 4, delta = 1, 
               method = "GEB")
print(xtable(round(res$sm,3)))
# after apply a=1 and delta = 0.3
res <- localPP(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20),
               be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5), a = 1, delta = 0.3, 
               method = "GEB")
print(xtable(round(res$sm,3)))
