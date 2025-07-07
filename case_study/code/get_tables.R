# Data Analysis for BRAF-V600 study
# Generates Tables 4.1, 4.2 and A3 from the intermediate results. 

# Load required packages
library(xtable)
load("./intermediate_results/case_study.RData")

#################################
### Table 4.1 BRAF V600 trial data
################################
Table4.1 <- cbind(N, responses, round(responses/N, 3))
rownames(Table4.1) <- c("NSCLC", "CRC vemu", "CRC vemu+cetu", 
                        "Bile duct", "ECD or LCH", "ATC")
colnames(Table4.1) <- c("Sample Size", "Number of Responses", "Response rate")
write.csv(Table4.1, "./results/Table 4.1.csv")

#################################
### Table 4.2 Results on BRAF V600 trial data
################################
results <- matrix(NA, nrow = B, ncol = 6)
rownames(results) <- c("NSCLC", "CRC vemu", "CRC vemu+cetu", "Bile duct", "ECD or LCH", "ATC")
colnames(results) <- c("Q_ind", "Q_localPP", "TypeI_ind", "TypeI_localPP", "pp_ind", "pp_localPP")
results[,1] <- sprintf(fmt = '%#.3f', Q.Independent)
results[,2] <- sprintf(fmt = '%#.3f', Q.localPP)
results[,3] <- sprintf(fmt = '%#.3f', oc.Independent[1,])
results[,4] <- sprintf(fmt = '%#.3f', oc.localPP[1,])
results[,5] <- sprintf(fmt = '%#.3f', pp.Independent)
results[,6] <- sprintf(fmt = '%#.3f', pp.localPP)
results
print(xtable(results))
write.csv(results, "./results/Table 4.2.csv")

#################################
### Table A3: Estimated similarity matrix on BRAF V600 trial data
################################
sss <- round(fit$sm, 2)
rownames(sss) <- c("NSCLC", "CRC vemu", "CRC vemu+cetu", "Bile duct", "ECD or LCH", "ATC")
print(xtable(sss)) # Appendix Table A3
fit$sm%*%N-N ## borrowing amount in term of number of subjects
write.csv(sss, "./results/Table A3.csv")
