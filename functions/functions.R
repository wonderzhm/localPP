
### Upload Required Packages
library(parallel)
library(doParallel)
library(R2jags)
library(foreach)
library(bhmbasket) ## required for EXNEX
library(basket) ## required for MEM (Kane et al., 2020)
library(BasketTrial) # for local-PP, Independent, and JSD methods
source(file.path("..", "functions", "utils.R"))   # required for BCHM (Chen and Lee, 2020)
source(file.path("..", "functions", "cluster.R")) # required for local MEM (Liu et al. 2022)

#' logit function
#'
#' @param x a value between 0 and 1
#'
#' @return returns log(x/(1-x)).
#' @export
#'
#' @examples
logit<-function(x)
{
  log(x/(1-x))
}

#' MEM (Kane et al., 2020)
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param p0 the response rate under the null hypothesis.
#' @param shape1 the first shape parameter for the prior of each basket.
#' @param shape2 the second shape parameter for the prior of each basket.
#'
#' @return It returns the estimated response rate and the posterior probability
#' of response rate greater than p0 for each basket.  
#' @export
#'
#' @examples
MEM <- function(nDat, yDat, p0 = 0.15, shape1 = 0.15, shape2 = 0.85){
  n <- nDat
  y <- yDat
  B <- length(n)
  res <- mem_exact(responses = y, size = n, name = letters[1:B], p0 = p0, 
                   shape1 = shape1,  shape2 = shape2)
  return(list(phat = as.vector(res$basket$mean_est), 
              post_prob = as.vector(res$basket$post_prob)))
}

#' Local MEM (Liu et al., 2022)
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param be.a0 the first shape parameters for the prior of each basket.
#' @param be.b0 the second shape parameters for the prior of each basket.
#' @param delta the delta parameter used in the prior; see Liu et al. (2022). 
#' @param max_C the maximum # of clusters. 
#'
#' @return It returns the two shape parameters of the posterior beta distribution of response rate for each basket. 
#' @export
#'
#' @examples
localMEM <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, 
                     delta = 2, # Delta Specification for the Prior
                     max_C = length(nDat) # Max # of Clusters
                     ){
  n <- nDat
  y <- yDat
  B <- length(n)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  if(B==1){
    return(list(a.post = be.a0 + y, b.post = be.b0 + n - y))
  }else{
    ## Partition
    part <- get.part(R = B, max_cl = max_C)
    K <- nrow(part)
    # Number of Blocks/Unique Response Rate in Each Partition
    n_bk <- apply(part, 1, function(x){
      length(unique(x))
    })
    prior_part <- n_bk^delta/sum(n_bk^delta)
    post_part <- update.part(x = y, n = n, prior_part = prior_part, 
                             part = part, a0 = 0.5, b0 = 0.5)
    fpart <- unlist( part[which.max(post_part), ])
    fmat <- NULL
    for(r in 1:B){
      fmat <- rbind(fmat, ifelse(fpart[r]==fpart, 1, 0))
    }
    diag(fmat) <- 0
    pOmega <- max(post_part)
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+y[i]+sum(fmat[i,]*y)*pOmega
      ab.post[i,2] <- be.b0[i]+n[i]-y[i]+sum(fmat[i,]*(n-y))*pOmega
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2]))
  }
}

#' BHM uniform (Cunanan et al., 2019)
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param mu0 hyperprior for mu~N(mu0, sig02)
#' @param sig02 hyperprior for mu~N(mu0, sig02)
#' @param u0 hyperprior for sig~ U(0, u0)
#' @param be.a0 Beta(a0, b0) on ORR for each basket for clustering
#' @param be.b0 Beta(a0, b0) on ORR for each basket for clustering
#' @param n.chains number of Markov chains. 
#' @param n.iter number of total iterations per chain (including burn in)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning. 
#' @param n.thin thinning rate. Must be a positive integer. 
#'
#' @return It results the posterior samples of response rate for each basket. 
#' @export
#'
#' @examples
BHMunif <- function(nDat,
                yDat,
                mu0 = 0, # hyperprior for mu~N(mu0, sig02)
                sig02 = 100, #hyperprior for mu~N(mu0, sig02)                 
                u0 = 100, # hyperprior for sig~ U(0, u0) 
                be.a0 = NULL, # Beta(a0, b0) on ORR for each basket for clustering
                be.b0 = NULL, # Beta(a0, b0) on ORR for each basket for clustering
                n.chains = 2,
                n.iter = 12000, 
                n.burnin = 2000,
                n.thin = 5)
{
  n <- nDat
  y <- yDat
  B <- length(n)
  nsave <- (n.iter-n.burnin)/n.thin*n.chains
  
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  
  psamples <- matrix(NA, nsave, B)
  if(B==1){
    psamples[,1] <- rbeta(nsave, be.a0+y, be.b0+n-y)
  }else{
    jagsmodel <- function()
    {
      for (j in 1:B){
        y[j] ~ dbinom(p[j], n[j]) # Likelihood
        logit(p[j]) <- theta[j]
        theta[j] ~ dnorm(mu, tau2) 
      }
      mu ~ dnorm(mu0, 1/sig02)         # Prior on mu
      sig ~ dunif(0, u0)    # Prior on sig
      tau2 <- 1/sig^2
    }
    # Keep track of p
    jags.params<-c("p")
    # Initial Value for Sampling Scheme
    jags.inits <-function(){
      list("mu" = 0, "sig" = 1)
    }
    jags.data <- list("B", "y", "n", "mu0", 
                      "sig02", "u0")
    jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                    n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                    model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
    psamples <- jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE]
  }
  return(psamples)
}

#' BHM gamma prior (Berry et al., 2013)
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param mu0 Hyperprior Mean
#' @param sig02 Hyperprior Variance
#' @param a0 Hyperprior Gamma(a0, b0) for 1/Variance
#' @param b0 Hyperprior Gamma(a0, b0) for 1/Variance
#' @param n.chains number of Markov chains. 
#' @param n.iter number of total iterations per chain (including burn in)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning. 
#' @param n.thin thinning rate. Must be a positive integer. 
#'
#' @return It results the posterior samples of response rate for each basket. 
#' @export
#'
#' @examples
BHM <- function(nDat,
                yDat,
                mu0 = -1.735, # Hyperprior Mean
                sig02 = 10, # Hyperprior Variance                 
                a0 = 2, # Hyperprior Gamma(a0, b0) for 1/Variance
                b0 = 20, # Hyperprior Gamma(a0, b0) for 1/Variance
                n.chains = 2,
                n.iter = 12000, 
                n.burnin = 2000,
                n.thin = 5)
{
  n <- nDat
  y <- yDat
  numGroups <- length(n)
  
  MeanMu <- mu0       
  PercMu <- 1/sig02  
  TauAlpha <- a0
  TauBeta <- b0
  
  jagsmodel <- function()
  {
    mu ~ dnorm(MeanMu, PercMu)         # Prior on Mu
    tau ~ dgamma(TauAlpha, TauBeta)    # Prior on Tau
    
    for (j in 1:numGroups){
      y[j] ~ dbinom(p[j], n[j]) # Likelihood
      theta[j] ~ dnorm(mu, tau) 
      p[j] <- 1/(1 + exp(-theta[j]))
    }
    
  }
  # Keep track of p
  jags.params<-c("p")
  # Initial Value for Sampling Scheme
  mu.init <- (sum(y)/sum(n) + 0.01)
  jags.inits <-function(){
    list("mu" = ifelse(mu.init < 0.99, logit(mu.init), 4.59512), "tau" = 1)
  }
  jags.data <- list("numGroups", "y", "n", "MeanMu", 
                    "PercMu", "TauAlpha", "TauBeta")
  jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                  n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
                  model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
  return(jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE])
}

#' This function generates samples from BCHM (Chen and Lee 2020).
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param mu Hyperprior mean for the cluster.
#' @param sigma02 Hyperprior variance for the cluster.
#' @param sigmaD2 Variance of subgroup response rate.
#' @param alpha Alpha value of the Dirichlet Process determining number of clusters.
#' @param d0 Minimum value for the similarity matrix.
#' @param alpha1 Prior for borrowing strength gamma(alpha1, beta1) in the hierarchical model.
#' @param beta1 Prior for borrowing strength gamma(alpha1, beta1) in the hierarchical model.
#' @param mu0 Hyperprior mean parameter of subgroup means in the hierarchical model
#' @param tau2 Hyperprior precision parameter of subgroup means in the hierarchical model
#' @param n.chains number of Markov chains. 
#' @param n.iter number of total iterations per chain (including burn in)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning. 
#' @param n.thin thinning rate. Must be a positive integer. 
#'
#' @return It results the posterior samples of response rate for each basket. 
#' @export
#'
#' @examples
BCHM <- function(nDat, yDat,
                 mu = 0.2,
                 sigma02 = 10, # Variance of Prior
                 sigmaD2 = 0.001, # Variance of Each Data Point 
                 alpha = 1e-40, # Similarity Constant
                 d0 = 0, # Lower Bound on Similarity Matrix
                 alpha1 = 50, # for Gamma(alpha1, beta1) in BHM
                 beta1 = 10, # for Gamma(alpha1, beta1) in BHM
                 mu0 = -1.735, # for the Hyperprior Mean Mean in BHM
                 tau2 = 0.1, # for the Hyperprior Mean in BHM Variance/Precision
                 n.chains = 2,
                 n.iter = 12000,
                 n.burnin = 2000,
                 n.thin = 5) {
  numArm <- length(nDat)
  nsave <- (n.iter-n.burnin)/n.thin*n.chains
  if (numArm != length(yDat)) {
    stop("Numbers of subgroups in nDat and yDat are not equal.")
  }
  if (numArm > 20) {
    stop("Numbers of subgroups is more than 20.")
  }
  if (numArm==1){
    allPost <- BHM(nDat = nDat, yDat = yDat, mu0 = mu0, sig02 = 1/tau2, a0 = alpha1,  
                   b0 = beta1, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                   n.thin = n.thin)
  } else{
    weight <- nDat
    alphaP <- alpha
    alpha <- alpha1
    beta <- beta1
    posi <- yDat
    priorMean <- mu
    priorVar <- sigma02
    res <- posi/weight
    x <- as.matrix(res)
    estVarGroup <- sigmaD2
    result <- gibbsSampler(x, alphaP, priorMean, priorVar, estVarGroup, 
                           weight, 1000, 2000)
    tables <- result$tables
    sm <- matrix(0, numArm, numArm)
    tSize <- dim(tables)[1]
    for (i in 1:tSize) {
      rr <- tables[i, ]
      for (j in 1:numArm) {
        for (k in 1:numArm) {
          if (rr[j] == rr[k]) {
            sm[j, k] <- sm[j, k] + 1
          }
        }
      }
    }
    sm <- sm/tSize
    smMinR <- d0
    sm[sm < smMinR] <- smMinR
    SMatrix <- sm
    jagsmodel <- function()
    {
      for (i in 1:numGroups)
      {
        y[i] ~ dbin(p[i],n[i]) 
        logit(p[i]) <- eta[i] 
        rr[i]<-tau1 * m[i]
        eta[i] ~ dnorm(mu,rr[i])
      }
      # Priors
      mu ~ dnorm(mu0, tau2)
      tau1 ~ dgamma(alpha, beta)
    }
    
    if(all(sm > 0.99)){
      allPost <- BHM(nDat = weight, yDat = posi, mu0 = mu0, sig02 = 1/tau2, a0 = alpha1,  
                     b0 = beta1, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
                     n.thin = n.thin)
    } else {
      for (i in 1:numArm) {
        smP <- sm[i, ]
        ind <- which(smP>0.01)
        n.ind <- length(ind)
        mydata <- list(y = posi[ind], n = weight[ind], m = smP[ind], 
                       numGroups = n.ind, mu0 = mu0, 
                       tau2 = tau2, alpha = alpha, beta = beta)
        # Parameters of Interest
        jags.params<-c("p")
        # Initial Value 
        mu.init <- sum(posi[ind])/sum(weight[ind]) + 0.01
        jags.inits<-function(){
          list("mu" = ifelse(mu.init<0.99, logit(mu.init), 4.59512), "tau1" = 1)
        }
        jagsfit <- jags(data = mydata, inits = jags.inits, parameters.to.save = jags.params, n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = jagsmodel, quiet = TRUE, progress.bar = "none")
        samples <- jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE]
        sampledP <- samples[, which(ind==i)]
        
        if (i == 1) {
          allPost <- sampledP
        }
        else {
          allPost <- cbind(allPost, sampledP)
        }
      }
    }
  }
  allPost
}

#' This function generates samples from EXNEX (Neuenschwander et al., 2016).
#'
#' @param nDat a vector with each element being the sample size for each basket. 
#' @param yDat a vector with each element being the number of responses for each basket. 
#' @param weight the weights for EX and NEX components
#' @param nex.mean mean of the normal prior distribution of the NEX part.
#' @param nex.prec precision of the normal prior distribution of the NEX part.
#' @param mu.mean mean of the normal prior distribution of the EX part.
#' @param mu.prec precision of the normal prior distribution of the EX part.
#' @param tau.HN.scale Prior Tau Precision Power Parameter
#' @param n.chains number of Markov chains. 
#' @param n.iter number of total iterations per chain (including burn in)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning. 
#' @param n.thin thinning rate. Must be a positive integer. 
#'
#' @return It results the posterior samples of response rate for each basket. 
#' @export
#'
#' @examples
EXNEX <- function(nDat, yDat, weight = c(0.25, 0.25, 0.5), 
                  nex.mean = -1.73, # NeX.Mean
                  nex.prec = 0.09, # Nex.Prec
                  mu.mean = c(-2.20, -2.20), # Mu.Mean in EX
                  mu.prec = c(0.099, 0.099), # Mu.Prec in EX
                  tau.HN.scale = c(1,1), # from the Paper - Prior Tau Precision Power Parameter
                  n.chains = 2,
                  n.iter = 12000,
                  n.burnin = 2000,
                  n.thin = 5) {
  Nexch = length(mu.mean)
  Nmix = Nexch + 1
  
  EXNEX_Model <- function(){
    
    # prior distributions for EX-parameters
    for (jj in 1:Nexch){
      mu[jj] ~ dnorm(mu.mean[jj], mu.prec[jj])
      prior.tau.prec[jj] <- pow(tau.HN.scale[jj],-2)
      tau[jj] ~ dnorm(0, prior.tau.prec[jj]); T(0.001, )
      prec.tau[jj] <- pow(tau[jj],-2)
    }
    
    # log-odds parameters under EX
    for (jj in 1:Nexch){
      for (j in 1:Nstrata){
        re[jj,j] ~ dnorm(0, prec.tau[jj])
        LogOdds[jj,j] <- mu[jj]+re[jj,j]
      }
    }
    
    # log-odds parameters under NEX
    for (j in 1:Nstrata){
      LogOdds[Nmix,j] ~ dnorm(nex.mean, nex.prec)
    }
    
    # latent mixture indicators:
    # exch.index: categorial 1,...,Nmix=Nexch+1
    # exch: Nstrata x Nmix matrix of 0/1 elements
    for (j in 1:Nstrata){
      exch.index[j] ~ dcat(pMix[1:Nmix])
      ## Below is different from Appendix of Neuenschwander, since they are not used.  
      #for (jj in 1:Nmix) {
      #  exch[j, jj] <- equals(exch_index[j], jj)
      #}
    }
    
    # pick theta
    for (j in 1:Nstrata){
      theta[j] <- LogOdds[exch.index[j],j]
    }
    
    # likelihood part
    for (i in 1:Nstrata){                
      logit(p[i]) <- theta[i]
      # p_success[i] <- step(p[i] - p_cut) ## Deviation from Appendix of Neuenschwander: since we don't need this. 
      y[i] ~ dbin(p[i],n[i])
    }
  }
  
  # Initial Values
  jags.params <- c("p")
  
  jags.inits <- function(){
    list(mu = mu.mean)
  }
  
  Nstrata <- length(nDat)
  n <- nDat
  y <- yDat
  pMix <- weight
  
  jags.data <- list("y","n","Nstrata","Nexch","Nmix","nex.mean","nex.prec", "pMix", "mu.mean", "mu.prec", "tau.HN.scale")
  
  jagsfit <- jags(data = jags.data, inits = jags.inits, parameters.to.save = jags.params,
                  n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, model.file = EXNEX_Model, quiet = TRUE, progress.bar = "none")
  
  return (jagsfit$BUGSoutput$sims.matrix[,-1, drop = FALSE])
}




