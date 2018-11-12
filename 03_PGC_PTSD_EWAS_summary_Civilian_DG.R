########################################################################################
# PGC EWAS Meta-Analysis: Summarizing Results
########################################################################################

rm(list=ls())
library(rmeta)
library(ChAMP)
library(forestplot)
library(qqman)
library(ggplot2)

########################################################################################
# Step 1: Loading Data
########################################################################################

# Combining results files
load("PGC_EWAS_nonSmoke_inverseNorm_intersectSites_civilian_DG_10.05.17.Rdata")
load("PGC_EWAS_nonSmoke_inverseNorm_remainingSites_civilian_DG_10.05.17.Rdata")

Studies<-rep("DNHS, GTP", nrow(res))
res<-cbind(res, Studies)
rm(Studies)
CpG<-rownames(res)
res<-cbind(CpG, res)
rm(CpG)

# Removing CpG sites not in at least two studies
dim(results)
results<-results[results[, "Studies"]!="DNHS",]
results<-results[results[, "Studies"]!="GTP",]

res<-res[, c("CpG", "Studies", "z.combined", "p")]
rownames(results)<-results[, "CpG"]
results<-results[, colnames(res)]
results<-rbind(res, results)
rm(res)

results<-results[order(as.numeric(results[, "p"])), ]
FDR<-p.adjust(as.numeric(results[, "p"]), method="fdr", n=nrow(results))
results<-cbind(results, FDR)

pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1) 
lambda # 1.052623

save(results, file="PGC_EWAS_civilian_DG_inverseNorm_allResults.Rdata")

# Summary info
tab<-matrix(nrow=6, ncol=2)
colnames(tab)<-c("Parameter", "Value")
tab[, "Parameter"]<-c("CpG sites", "N sites p<5x10^-5", "N sites p<5x10^-6",
                      "N sites p<5x10^-7", "N sites FDR < 0.05", "lambda")
rownames(tab)<-tab[,"Parameter"]
pvalue <- as.numeric(results[, "p"])
chisq <- qchisq(1-pvalue,1)
lambda = median(chisq)/qchisq(0.5,1)
tab["lambda", "Value"]<-round(lambda,6)

tab["N sites FDR < 0.05", "Value"]<-sum(FDR<=0.05)
tab["N sites p<5x10^-5", "Value"]<-sum(pvalue<=5*10^-5)
tab["N sites p<5x10^-6", "Value"]<-sum(pvalue<=5*10^-6) # 7
tab["N sites p<5x10^-7", "Value"]<-sum(pvalue<=5*10^-7) # 3
tab["CpG sites", "Value"]<-nrow(results)
rownames(tab)<-c(1:nrow(tab))

write.csv(tab, file="PGC_EWAS_civilian_DG_inverseNorm_allResults_summary.csv")
rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 2: QQ plot
########################################################################################

ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=2, cex.lab=4.5, cex.axis=4, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}
pvalue <- as.numeric(results[, "p"])

png("PGC_EWAS_inverseNorm_QQplot_civilian_DG.png", width=1600, height=800, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()

rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 3: Manhattan plot
########################################################################################

rm(list=ls())

load("PGC_EWAS_civilian_DG_inverseNorm_allResults.Rdata")
data(probe.features)
results<-results[which(rownames(results)%in%rownames(probe.features)),]
probe.features$CHR<-as.character(probe.features$CHR)
probe.features<-probe.features[rownames(results),]
probe.features$CHR[probe.features$CHR=="X"]<-23
probe.features$CHR[probe.features$CHR=="Y"]<-24
probe.features$CHR<-as.numeric(probe.features$CHR)
probe.features<-as.matrix(probe.features[, c("CHR", "MAPINFO")])
all(rownames(results)==rownames(probe.features))
results<-cbind(results, probe.features)

head(results)
results<-results[, c("p", "FDR", "CHR", "MAPINFO")]
colnames(results)<-c("P", "FDR", "CHR", "BP")
head(results)
str(results)
class(results)<-"numeric"
head(results)
str(results)

cutpoint<--log10(1.210462*10^-06)
df<-data.frame(results)
labs<-c(1:22, "X", "Y")
df<-df[-which(is.na(df[, "CHR"])),]
opar<-par()

png("PGC_EWAS_manhattan_civilan_DG_wINTRuST.png",height=600, width=1200, units="px")
par(mar=c(5.1,6.1, 4.1, 2.1))
manhattan(df, suggestiveline=FALSE, genomewideline=cutpoint, ylim=c(0,8.5),
          chrlabs=labs, cex=2, cex.axis=2, cex.lab=3, col=c("blue4", "red4"))
dev.off()

rm(list=ls())

########################################################################################
# Step 4: Top Results Table
########################################################################################
rm(list=ls())
load("PGC_EWAS_civilian_DG_inverseNorm_allResults.Rdata")

rownames(results)<-results[, "CpG"]
cpgs<-rownames(results[as.numeric(results[, "p"])<=5*10^-5, ]) # subsetting only most significant results
topResults<-read.csv("/Users/aratanat/Documents/R/PGC_EWAS/Current/PGC_EWAS/nonSmoke/wINTRuST/PGC_EWAS_topresults.csv",
                     stringsAsFactors = F, row.names=1)
topResults<-rownames(topResults)
topResults<-topResults[!topResults%in%cpgs]
results<-results[c(cpgs, topResults), ] # subsetting only most significant results
results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

load("PGC_EWAS_DataPrep_nonSmoke_civilian_DG.Rdata")
rm(list=ls()[grep("ebayes", ls())])

colnames(DNHS.coef)
DNHS.coef<-DNHS.coef[rownames(DNHS.coef)%in%sites,c("PTSDpm", "N.subjects")]
DNHS.results<-DNHS.results[rownames(DNHS.coef),]
colnames(DNHS.coef)<-c("PTSD", "N.subjects") # need consistent column headings
all(rownames(DNHS.coef)==rownames(DNHS.results))
DNHS<-data.frame(DNHS.coef, DNHS.results)
str(DNHS) # should all be numbers not factors
DNHS$s.e.<-DNHS$PTSD/DNHS$t # calculate the standard error
DNHS$weight<-1/(DNHS$s.e.^2) # calculate weight
rm(DNHS.results, DNHS.coef, DNHS.oneSided)

colnames(GTP.coef)
GTP.coef<-GTP.coef[rownames(GTP.coef)%in%sites,c("PTSDcurr", "N.subjects")]
GTP.results<-GTP.results[rownames(GTP.coef),]
colnames(GTP.coef)<-c("PTSD", "N.subjects")
all(rownames(GTP.coef)==rownames(GTP.results))
GTP<-data.frame(GTP.coef, GTP.results)
str(GTP)
GTP$s.e.<-GTP$PTSD/GTP$t # calculate the standard error
GTP$weight<-1/(GTP$s.e.^2) # calculate weight
rm(GTP.results, GTP.coef, GTP.oneSided)

studies<-c("DNHS", "GTP")

for(ii in 1:length(studies)){
  mat<-matrix(nrow=nrow(results), ncol=4)
  rownames(mat)<-rownames(results)
  colnames(mat)<-paste(studies[ii], c(".beta", ".se", ".p", ".N"), sep="")
  results<-cbind(results, mat)
  rm(mat)
}

results$variance<-results$beta<-NA

for(ii in 1:nrow(results)){
  cpg<-rownames(results)[ii]
  weights<-NULL
  betas<-NULL
  studies<-unlist(strsplit(results[ii, "Studies"], ", "))
  for(jj in 1:length(studies)){
    temp<-get(paste(studies[jj]))
    weights<-append(weights, temp[cpg, "weight"])
    betas<-append(betas, temp[cpg, "PTSD"])
    results[cpg, paste(studies[jj], ".beta", sep="")]<-temp[cpg, "PTSD"]
    results[cpg, paste(studies[jj], ".se", sep="")]<-temp[cpg, "s.e."]
    results[cpg, paste(studies[jj], ".p", sep="")]<-temp[cpg, "P.Value"]
    results[cpg, paste(studies[jj], ".N", sep="")]<-temp[cpg, "N.subjects"]
  }
  results[cpg, "beta"]<-sum(betas*weights)/sum(weights)
  results[cpg, "variance"]<-1/sum(weights)
}

data(probe.features)
probe.features<-probe.features[rownames(results), c("CHR", "MAPINFO", "gene", "feature")]
probe.features<-cbind(results, probe.features)
probe.features<-probe.features[, c("Studies", "CpG", "CHR", "MAPINFO", "gene", "feature",
                                   "beta", "variance", "p", "FDR",
                                   paste(studies, ".beta", sep=""),
                                   paste(studies, ".se", sep=""),
                                   paste(studies, ".p", sep=""),
                                   paste(studies, ".N", sep=""))]
                                   
colnames(probe.features)<-c("Studies", "CpG", "CHR", "Position", "Gene", "Feature",
                            "beta", "variance", "p-value", "FDR",
                            paste(studies, ".beta", sep=""),
                            paste(studies, ".se", sep=""),
                            paste(studies, ".p", sep=""),
                            paste(studies, ".N", sep=""))

rownames(probe.features)<-c(1:nrow(probe.features))

write.csv(probe.features, "PGC_EWAS_civilian_DG_TopResults.csv")
rm(list=ls())

dg<-read.csv("/Users/aratanat/Documents/R/PGC_EWAS/Current/PGC_EWAS/Civilian/nonSmoke/DNHS_GTP/PGC_EWAS_civilian_DG_TopResults.csv", row.names=1)
rownames(dg)<-dg$CpG
load("/Users/aratanat/Dropbox/PTSD methylation workgroup Results Files/Non_Smoking/WTC180_meta_limma_results_08.09.16.Rdata")
Rank<-c(1:nrow(results))
results<-cbind(results, Rank)
res<-results[which(as.numeric(results[, "P.Value"])<=5*10^-5), ]
write.csv(res, "/Users/aratanat/Documents/R/WTC_Meta/WTC224_meta/Final_08.16/WTC180_topResults.csv")
rm(res)

results<-results[rownames(dg),]
results<-results[rownames(dg),c("logFC", "t", "P.Value", "adj.P.Val", "Rank")]
head(results)
colnames(results)<-c("WTC_beta", "WTC_t", "WTC_p", "WTC_FDR", "WTC_Rank")
dg<-cbind(dg, results)
View(dg)
write.csv(dg, "PGC_EWAS_civilian_DG_TopResults_withWTC.csv")
