#' Posterior Inference for Simulated Basket Trial Data
#' 
#' Generates posterior probabilities P(p_j > pnull) after all interim analysis and 
#' calculates rates for early stopping, number of patients and estimated ORR. It is a parallel computing version of 
#' \code{BasketTrial::post.infer}. 
#' @param nclust number of clusters for parallel computing. 
#' @param nperclust number of simulated trials per cluster. 
#' @param data.object object returned from \code{BasketTrial::generate.data}.
#' @param pnull B by 1 vector of null response rates, where B is the number of baskets.
#' @param stopbounds B by (stage-1) matrix: stopping boundaries for each basket at each interim.
#' @param beta.a0 a vector of length B for beta prior parameter a0 in each basket.
#' @param beta.b0 a vector of length B for beta prior parameter b0 in each basket.
#' @param seed random seed for reproducibility.
#' @param ModelFit	the method function, e.g., localPP, JSD, and other user defined methods.
#' @param ... additional arguments passed to the method function defined by \code{ModelFit}.
#'
#' @return  It returns a list including early stopping, number of patients and estimated ORR.
#' @export
#'
#' @examples
post.infer.parallel <- function(nclust, nperclust, data.object, pnull, stopbounds = NULL, 
                                beta.a0 = pnull, beta.b0 = 1-pnull, seed = 38098, ModelFit, ...) {
  if(nclust > detectCores()) stop("nclust needs to be smaller than detectCores()")
  
  ## Create Clusters
  cl <- makeCluster(nclust, type = "PSOCK")
  registerDoParallel(cl)
  
  ## Export Functions to the Cluster
  clusterCall(cl, function() {source(file.path("..", "functions", "functions.R"))})
  clusterExport(cl, list("nperclust", "data.object", "pnull", "stopbounds", "beta.a0", "beta.b0", "seed"))
  
  ## Parallel Computing
  results <- parLapply(cl, (1:nclust), fun = post.infer, object = data.object, 
                       pnull = pnull, stopbounds = stopbounds, nperclust = nperclust, 
                       beta.a0 = beta.a0, beta.b0 = beta.b0, seed = seed,
                       ModelFit = ModelFit, ...)
  
  ## Close Cluster
  stopCluster(cl)
  
  ## Combine the Results
  resultsAll <- results[[1]]
  for(i in 2:nclust){
    resultsAll$earlystop <- abind::abind(resultsAll$earlystop, results[[i]]$earlystop, along=2)
    resultsAll$postprob <- abind::abind(resultsAll$postprob, results[[i]]$postprob, along=2)
    resultsAll$npts <- abind::abind(resultsAll$npts, results[[i]]$npts, along=2)
    resultsAll$est <- abind::abind(resultsAll$est, results[[i]]$est, along=2)
  }
  return(resultsAll)
}

