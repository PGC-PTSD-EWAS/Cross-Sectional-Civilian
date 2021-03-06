########################################################################################
# PGC EWAS Civilian Female Meta-Analysis: Summarizing Results
########################################################################################

rm(list=ls())
library(rmeta)
library(ChAMP)
library(forestplot)
library(qqman)

setwd("/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/Civilian/Female/")

########################################################################################
# Step 1: Loading Data
########################################################################################

# Combining results files
load("PGC_EWAS_nonSmoke_inverseNorm_intersectSites_civilian_female.Rdata")
load("PGC_EWAS_nonSmoke_inverseNorm_remainingSites_civilian_female.Rdata")

all<-pval.combined.marot
rm(pval.combined.marot)
colnames(all)<-c("CpG", "p")
Studies<-rep(c("DNHS, GTP"), nrow(all))
all<-cbind(all, Studies)
all<-all[, colnames(results)]
all(colnames(all)==colnames(results))
results<-rbind(all, results)

# Removing CpG sites not in at least two studies
dim(results)
results<-results[results[, "Studies"]!="DNHS",]
results<-results[results[, "Studies"]!="GTP",]
dim(results)

results<-results[order(as.numeric(results[, "p"])), ]
FDR<-p.adjust(as.numeric(results[, "p"]), method="fdr", n=nrow(results))
results<-cbind(results, FDR)
save(results, file="PGC_EWAS_inverseNorm_allResults_civilian_female.Rdata")

# Summary info
tab<-matrix(nrow=9, ncol=2)
colnames(tab)<-c("Parameter", "Value")
tab[, "Parameter"]<-c("CpG sites", "N sites in 3 studies", "N sites in 2 studies",
                      "N sites in 1 study", "N sites p<5x10^-5", "N sites p<5x10^-6",
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

nStudies<-unlist(lapply(strsplit(results[, "Studies"], ", "), function(x) length(x)))
tab["N sites in 3 studies", "Value"]<-sum(nStudies==3)
tab["N sites in 2 studies", "Value"]<-sum(nStudies==2)
tab["N sites in 1 study", "Value"]<-sum(nStudies==1)
tab["CpG sites", "Value"]<-nrow(results)
rownames(tab)<-c(1:nrow(tab))

write.csv(tab, file="PGC_EWAS_inverseNorm_allResults_civilian_female_summary.csv")
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

png("PGC_EWAS_inverseNorm_QQplot_civilian_female.png", width=1280, height=960, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()

png("PGC_EWAS_inverseNorm_QQplot_civilian_forPublication.png", width=1600, height=800, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()


rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 3: Manhattan plot
########################################################################################

rm(list=ls())

load("PGC_EWAS_inverseNorm_allResults_civilian_female.Rdata")
rownames(results)<-results[, "CpG"]
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

cutpoint<--log10(0.05/nrow(results))
df<-data.frame(results)
labs<-c(1:22, "X", "Y")

opar<-par()

png("PGC_EWAS_manhattan_civilian_female.png",height=600, width=1200, units="px")
par(mar=c(5.1,6.1, 4.1, 2.1))
manhattan(df, suggestiveline=FALSE, genomewideline=cutpoint, ylim=c(0,8),
          chrlabs=labs, cex=1, cex.axis=1.5, cex.lab=2, col=c("blue4", "red4"))
dev.off()

rm(list=ls())

########################################################################################
# Step 4: Top Results Table
########################################################################################

load("PGC_EWAS_inverseNorm_allResults_civilian_female.Rdata")
results<-results[as.numeric(results[, "p"])<=5*10^-5, ] # subsetting only most significant results
rownames(results)<-results[, "CpG"]
results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

load("PGC_EWAS_DataPrep_nonSmoke_civilian_female.Rdata")
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

results$GTP.beta<-results$DNHS.beta<-NA
results$GTP.se<-results$DNHS.se<-NA
results$GTP.N<-results$DNHS.N<-NA
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
                                   "DNHS.beta", "GTP.beta", 
                                   "DNHS.se", "GTP.se", 
                                   "DNHS.N", "GTP.N")]
colnames(probe.features)<-c("Studies", "CpG", "CHR", "Position", "Gene", "Feature",
                            "beta", "variance", "p-value", "FDR",
                            "DNHS.beta", "GTP.beta",
                            "DNHS.se", "GTP.se", 
                            "DNHS.N", "GTP.N")

rownames(probe.features)<-c(1:nrow(probe.features))
write.csv(probe.features, "PGC_EWAS_TopResults_civilian_female.csv")
rm(list=ls())

results<-read.csv("PGC_EWAS_TopResults_civilian_female.csv", row.names=1, stringsAsFactors=F)
studies<-c("DNHS", "GTP")

pdf("PGC_EWAS_Civilian_female_forestPlots.pdf")
for(ii in 1:nrow(results)){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", studies, NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  forestplot(summary, cochrane,
             #boxsize=0.5, lwd.ci=1.5,
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,rep(FALSE,4), TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="black",line="black", summary="royalblue"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=1)),
                               ticks = gpar(fontfamily = "", cex=1)))
}
dev.off()

# Significant Sites
sigs<-length(results[as.numeric(results[, "FDR"]<=0.05), "FDR"])

for(ii in 1:sigs){
  betas<-as.numeric(results[ii,paste(studies, ".beta", sep="")])
  stderr<-as.numeric(results[ii,paste(studies, ".se", sep="")])
  lower<-betas-1.96*stderr
  upper<-betas+1.96*stderr
  betaC<-results[ii, "beta"]
  stderrC<-results[ii, "variance"]
  lowerC<-betaC-(1.96*sqrt(stderrC))
  upperC<-betaC+(1.96*sqrt(stderrC))
  beta.table<-c("", "Beta", round(betas,3), NA, round(betaC,3))
  subjs<-results[ii,paste(studies, ".N", sep="")]
  Ns<-as.character(unlist(c("", "N", subjs, NA, sum(subjs))))
  betas<-c(NA, NA, round(betas,3), NA, round(betaC,3))
  lower<-c(NA, NA, lower, NA, lowerC)
  upper<-c(NA, NA, upper, NA, upperC)
  stud<-c(results[ii, "CpG"], "Study", studies, NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  png(paste("PGC_EWAS_Civilian_female_", results[ii, "CpG"], ".png", sep=""), width=960, height=640,units="px")
  forestplot(summary, cochrane,
             boxsize=c(0.8, 0.8, 0.4, 0.6, 0.5), lwd.ci=1.5,
             new_page = TRUE,
             is.summary=c(TRUE,TRUE,rep(FALSE,3),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="black",line="black", summary="royalblue"),
             xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=3)),
                               ticks = gpar(fontfamily = "", cex=2)))
  dev.off()
}