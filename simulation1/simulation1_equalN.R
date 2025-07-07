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
source(file.path("..", "functions", "functions.R"))
source(file.path("..", "functions", "functions_parallel.R"))

####### Simulation Settings and Generate Data ##############
nperclust <- 500 # number of simulated trials per cluster 
detectCores()
nclust <- 10 # Total 5000 MC replicates
source(file.path(".", "code", "sim_settings.R"))

####################################################################################
## Simulation results for comparing the  following eight methods:
## Independent, local PP, JSD, EXNEX, BHM, BCHM, local MEM and MEM
## This will take about 20 hours depending on the computer power.
## The intermediate RData has been saved into "intermediate_results" which
## can be used directly to generate all the Tables in the paper. 
####################################################################################

####################### Independent Model ########################
source(file.path(".", "code", "IM.R"))

####################### local-PP-PEB1 ########################
source(file.path(".", "code", "local-PP-PEB1.R"))

####################### local-PP-PEB2 ########################
source(file.path(".", "code", "local-PP-PEB2.R"))

####################### PP-PEB ########################
source(file.path(".", "code", "PP-PEB.R"))

####################### local-PP-GEB1 ########################
source(file.path(".", "code", "local-PP-GEB1.R"))

####################### local-PP-GEB2 ########################
source(file.path(".", "code", "local-PP-GEB2.R"))

####################### PP-GEB ########################
source(file.path(".", "code", "PP-GEB.R"))

####################### JSD1 ########################
source(file.path(".", "code", "JSD1.R"))

####################### JSD2 ########################
source(file.path(".", "code", "JSD2.R"))

####################### local MEM ########################
source(file.path(".", "code", "local-MEM.R"))

####################### BHM ########################
source(file.path(".", "code", "BHM.R"))

####################### EXNEX ########################
source(file.path(".", "code", "EXNEX.R"))

####################### BCHM ########################
source(file.path(".", "code", "BCHM.R"))

####################### MEM ########################
source(file.path(".", "code", "MEM.R"))

######## save above methods into RData file. 
save.image(file.path(".", "intermediate_results", "AllMethods_main.RData"))


####################################################################################
## Steps below read the AllMethods_main.RData file and 
## create Tables 3.2, 3.3, and A2 in the paper
####################################################################################
source(file.path(".", "code", "get_tables.R"))

