Supplementary information / reproducible research files for
Title: "A Bayesian Basket Trial Design Using Local Power Prior"

Authors: Haiming Zhou, Rex Shen, Sutan Wu and Philip He
Code was written by Haiming Zhou and Rex Shen. 
In case of questions or comments please contact haiming2019@gmail.com

Please install the R package BasketTrial by devtools::install_github("wonderzhm/BasketTrial").

The code was written/evaluated in R with the following software versions:
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] BasketTrial_0.1.0 cluster_2.1.6     plyr_1.8.9        partitions_1.10-7 basket_0.10.11    bhmbasket_0.9.5  
 [7] R2jags_0.8-5      rjags_4-16        coda_0.19-4.1     doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2    

loaded via a namespace (and not attached):
 [1] gtable_0.3.5       ggplot2_3.5.1      ggrepel_0.9.5      lattice_0.22-6     mathjaxr_1.6-0     Rdpack_2.6.1      
 [7] vctrs_0.6.5        tools_4.4.1        generics_0.1.3     tibble_3.2.1       fansi_1.0.6        pkgconfig_2.0.3   
[13] RColorBrewer_1.1-3 lifecycle_1.0.4    compiler_4.4.1     farver_2.1.2       stringr_1.5.1      munsell_0.5.1     
[19] ggforce_0.4.2      graphlayouts_1.2.0 codetools_0.2-20   gmp_0.7-5          pillar_1.9.0       crayon_1.5.3      
[25] tidyr_1.3.1        MASS_7.3-60.2      cachem_1.1.0       viridis_0.6.5      boot_1.3-30        abind_1.4-5       
[31] tidyselect_1.2.1   stringi_1.8.4      dplyr_1.1.4        purrr_1.0.2        polyclip_1.10-7    fastmap_1.2.0     
[37] grid_4.4.1         colorspace_2.1-1   cli_3.6.3          magrittr_2.0.3     ggraph_2.2.1       tidygraph_1.3.1   
[43] utf8_1.2.4         GenSA_1.1.14.1     withr_3.0.1        scales_1.3.0       igraph_2.0.3       gridExtra_2.3     
[49] memoise_2.0.1      rbibutils_2.2.16   viridisLite_0.4.2  rlang_1.1.4        itertools_0.1-3    Rcpp_1.0.13       
[55] glue_1.7.0         polynom_1.4-1      tweenr_2.0.3       rstudioapi_0.16.0  R6_2.5.1           R2WinBUGS_2.1-22.1


This folder contains the following data and files that can be used to reproduce all analysis results of the manuscript.
It contains four subfolders containing the following files:

./case_study/:
    code_case_study_alpha=0.05.R
    An R script that contains the code of the analysis reported in the paper (section 4) when type I error is controlled at 0.05. The resulted workspace is saved in case_study_alpha=0.05.RData.

./functions/:
    functions.R
    An R script that contains functions to implement MEM, local-MEM, BHM, BCHM and EXNEX on a single dataset. The methods for Independent, local-PP and JSD methods are implemented in the R package BasketTrial. 

    cluster.R
    An R script that contains functions required for local-MEM method, copied from Liu et al (2022).

    utils.R
    An R script that contains functions required for BCHM, copied from Chen and Lee (2020).

    functions_parallel.R
    An R script that contains functions for parallel computing in simulations.

./misc/:
    localPP_SimilarityMatrix.R
    An R script that was used to generate the similarity matrix of local-PP in the manuscript (Section 2.2).
    
./simulation/
    code_sim_tuning.R
    An R script that performs model parameter tuning for local-PP-PEB, local-PP-GEB and JSD methods (Section 3.2).
    
    code_sim_main.R
    An R script that performs the simulations reported in the manuscript (Section 3). Simulation was performed in parallel on 10 cores on Windows 10. Results of the simulation are saved into the folder ./intermediate_results/ since computation of the whole simulation can take around 2 days since some MCMC methods are time-consuming.
    For faster computation the number of replications within each simulation scenario could be reduced from nperclust=500. If the machine used for computation is not able to compute on 10 cores simultaneously the number of cores used for simulation should be reduced by changing the nclust argument.
    
    code_sim_supp.R
    An R script that performs the simulations with unequal basket sample sizes (Section 3.4). 
    For faster computation the number of replications within each simulation scenario could be reduced from nperclust=500. If the machine used for computation is not able to compute on 10 cores simultaneously the number of cores used for simulation should be reduced by changing the nclust argument.
    
    ./intermediate_results/
    A folder containing the results of code_sim_tuning.R, code_sim_main.R and code_sim_supp.R
    
    results_sim_main.R and results_sim_supp.R
    R scripts that take the RData files from the intermediate_results and creates the
    results tables presented in the manuscript.
    
