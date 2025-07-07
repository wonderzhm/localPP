################################################################################################
# This main R file re-produces all tables of the article by executing one single file. 
# Option 1: re-run all the models and re-produce the tables. 
# Option 2: Use the authors' previously saved intermediate workspace results to generate the tables. 
################################################################################################

######################################## Option 1 ##############################################
# If one wants to reduce the run time, please reduce the value of nperclust in source R files. 
# If the machine used for computation is not able to compute on 10 cores simultaneously, 
# the number of cores used for simulation should be reduced by reducing the value of nclust.
####### Generate Tables 3.1, 3.2, 3.3 and A2 ##############
## Results are saved to ./simulation1/results
setwd(file.path(".", "simulation1")) # set the current working directory to folder simulation1
source(file.path(".", "simulation1_equalN.R"))
setwd("..") # reset the current working directory to the parent folder that contains main.R. 

####### Generate Tables 3.4, 3.5, and 3.6 ##############
## Results are saved to ./simulation2/results
setwd(file.path(".", "simulation2"))  # Go to simulation2 folder
source(file.path(".", "simulation2_unequalN.R"))
setwd("..")  # Back to parent folder

####### Generate Tables 4.1, 4.2, and A3 ##############
## Results are saved to ./case_study/results
setwd(file.path(".", "case_study"))  # Go to case_study folder
source(file.path(".", "code_case_study.R"))
setwd("..")  # Back to parent folder

######################################## Option 2 ##############################################
####### Generate Tables 3.1, 3.2, 3.3 and A2 ##############
## Results are saved to ./simulation1/results
setwd(file.path(".", "simulation1"))  # Go to simulation1 folder
source(file.path(".", "code", "get_tables.R"))
setwd("..")  # Back to parent folder

####### Generate Tables 3.4, 3.5, and 3.6 ##############
## Results are saved to ./simulation2/results
setwd(file.path(".", "simulation2"))  # Go to simulation2 folder
source(file.path(".", "code", "get_tables.R"))
setwd("..")  # Back to parent folder

####### Generate Tables 4.1, 4.2, and A3 ##############
## Results are saved to ./case_study/results
setwd(file.path(".", "case_study"))  # Go to case_study folder
source(file.path(".", "code", "get_tables.R"))
setwd("..")  # Back to parent folder
