########################################################################################
# PGC Methylation Meta-Analysis
########################################################################################

rm(list=ls())

setwd("/Users/ar3054/Dropbox/PTSD methylation workgroup Results Files/Gender-Stratified/Non-Smoking/Male/")
homeDir<-"/Users/ar3054/Documents/R/PGC_EWAS/Current/PGC_EWAS/Civilian/Male"

########################################################################################
# Step 1: Load Summary Statistics 
########################################################################################

# Step 1A: Load DNHS Results
load("DNHS100_limma_results_malesOnly_08.09.16.Rdata")
DNHS.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
DNHS.results<-results
rm(results, rank)
DNHS.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1B: Load GTP Results
load("GTP_limma_results_males_02.26.17.Rdata")
GTP.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
GTP.results<-results
rm(results)
rank<-c(1:nrow(GTP.results))
GTP.results<-cbind(GTP.results, rank)
rm(rank)
GTP.ebayes<-fit2.ebayes
rm(fit2.ebayes)

# Step 1C: Load WTC Data
load("/Users/ar3054/Dropbox/PTSD methylation workgroup Results Files/Non_Smoking/WTC180_meta_limma_results_08.09.16.Rdata")
WTC.coef<-fit.coef
rm(fit.coef)
results<-results[order(results[, "P.Value"]), ]
rank<-c(1:nrow(results))
results<-cbind(results, rank)
WTC.results<-results
rm(results, rank)
WTC.ebayes<-fit2.ebayes
rm(fit2.ebayes)

setwd(homeDir)
save.image("PGC_EWAS_combinedData_civilian_males.Rdata")
rm(list=ls())

########################################################################################
# # Step 2: QQ Plots of Individual Data Sources
########################################################################################

load("PGC_EWAS_combinedData_civilian_males.Rdata")

studies<-c("DNHS", "GTP", "WTC")
tab<-matrix(nrow=length(studies), ncol=9)
colnames(tab)<-c("Study", "Sites", "Min N", "Max N", "Lambda", "FDR", "p<5x10^-5", "p<5x10^-6", "p<5x10^-7")

for(ii in 1:length(studies)){
  tab[ii, "Study"]<-studies[ii]
  res<-get(paste(studies[ii], ".results", sep=""))
  res<-res[!is.na(res[, "logFC"]),] # this is for VA Boston. Ther are no p-values for some sites. I'm guess there weren't any cases
  tab[ii, "p<5x10^-5"]<-sum(res[, "P.Value"]<(5*10^(-5))) 
  tab[ii, "p<5x10^-6"]<-sum(res[, "P.Value"]<(5*10^(-6))) 
  tab[ii, "p<5x10^-7"]<-sum(res[, "P.Value"]<(5*10^(-7))) 
  tab[ii, "FDR"]<-sum(res[, "adj.P.Val"]<=0.05)
  tab[ii, "Sites"]<-nrow(res) 
  coef<-get(paste(studies[ii], ".coef", sep=""))
  tab[ii, "Min N"]<-range(coef[, "N.subjects"])[1] 
  tab[ii, "Max N"]<-range(coef[, "N.subjects"])[2] 
  p<-res[, "P.Value"]
  chisq<-qchisq(1-p,1)
  tab[ii, "Lambda"]<-round(median(chisq)/qchisq(0.5,1),3)
  rm(p, chisq, coef, res)
}
tab
write.csv(tab, "PGC_EWAS_individualCohort_results_civilian_males.csv")
