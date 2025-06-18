####### Simulation Settings and Generate Data ##############
N <- rbind(c(10, 26),
           c(10, 16),
           c(4, 8), # the first number does not matter since no interim for this basket
           c(10, 17),
           c(10, 22)
) # interim sample size and total sample size for each indication
B <- nrow(N) # total number of baskets
p0 <- 0.15
p1 <- 0.30
p2 <- 0.45
pnull <- c(p0, p0, p0, p0, p0) # null response rate for each indication
ptarget <- c(p1, p1, p1, p1, p1) # target response rate for each indication

## type I error rate sig.level=0.1
sig.level <- 0.1 # type I error
stopbounds <- rbind(1, 
                    1, 
                    -1, # no interim for this basket due to the small maximum sample size
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
