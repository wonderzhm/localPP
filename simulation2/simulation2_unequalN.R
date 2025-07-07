# Replicate the simulation results under the unequal basket size setting
# It re-generates Tables 3.4, 3.5, and 3.6 in the paper

# Load required packages
library(parallel)
library(doParallel)
library(R2jags)
library(foreach)
library(bhmbasket) 
library(basket) 
library(partitions)
library(plyr)
library(cluster)
library(coda)
library(rjags)
library(BasketTrial)
source(file.path("..", "functions", "functions.R"))
source(file.path("..", "functions", "functions_parallel.R"))

####### Simulation Settings and Generate Data ##############
nperclust <- 500 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 5000 MC replicates
source(file.path(".", "code", "sim_settings.R"))

####################################################################################
## Simulation results for comparing the  following methods:
## Independent, local-PP-PEB, local-PP-GEB, JSD, EXNEX 
## The intermediate RData has been saved into "intermediate_results" which
## can be used directly to generate the Tables in the paper. 
####################################################################################

####################### Independent Model ########################
source(file.path(".", "code", "IM.R"))

####################### local-PP-PEB ########################
source(file.path(".", "code", "local-PP-PEB.R"))

####################### local-PP-GEB ########################
source(file.path(".", "code", "local-PP-GEB.R"))

####################### JSD ########################
source(file.path(".", "code", "JSD.R"))

####################### EXNEX ########################
source(file.path(".", "code", "EXNEX.R"))

######## save above methods into RData file.
save.image(file.path(".", "intermediate_results", "AllMethods_supp.RData"))

####################################################################################
## Steps below read the AllMethods_supp.RData file and 
## create Tables 3.4, 3.5, and 3.6 in the paper
####################################################################################
source(file.path(".", "code", "get_tables.R"))

