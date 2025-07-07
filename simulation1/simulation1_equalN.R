# Replicate the simulation results under the equal basket size setting
# It re-generates Tables 3.1, 3.2, 3.3 and A2 in the paper

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
source("../functions/functions.R")
source("../functions/functions_parallel.R")

####### Simulation Settings and Generate Data ##############
nperclust <- 500 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 5000 MC replicates
source("./code/sim_settings.R")

####################################################################################
## Simulation results for comparing the  following eight methods:
## Independent, local PP, JSD, EXNEX, BHM, BCHM, local MEM and MEM
## This will take about 20 hours depending on the computer power.
## The intermediate RData has been saved into "intermediate_results" which
## can be used directly to generate all the Tables in the paper. 
####################################################################################

####################### Independent Model ########################
source("./code/IM.R")

####################### local-PP-PEB1 ########################
source("./code/local-PP-PEB1.R")

####################### local-PP-PEB2 ########################
source("./code/local-PP-PEB2.R")

####################### PP-PEB ########################
source("./code/PP-PEB.R")

####################### local-PP-GEB1 ########################
source("./code/local-PP-GEB1.R")

####################### local-PP-GEB2 ########################
source("./code/local-PP-GEB2.R")

####################### PP-GEB ########################
source("./code/PP-GEB.R")

####################### JSD1 ########################
source("./code/JSD1.R")

####################### JSD2 ########################
source("./code/JSD2.R")

####################### local MEM ########################
source("./code/local-MEM.R")

####################### BHM ########################
source("./code/BHM.R")

####################### EXNEX ########################
source("./code/EXNEX.R")

####################### BCHM ########################
source("./code/BCHM.R")

####################### MEM ########################
source("./code/MEM.R")

######## save above methods into RData file. 
save.image("./intermediate_results/AllMethods_main.RData")


####################################################################################
## Steps below read the AllMethods_main.RData file and 
## create Tables 3.2, 3.3, and A2 in the paper
####################################################################################
source("./code/get_tables.R")