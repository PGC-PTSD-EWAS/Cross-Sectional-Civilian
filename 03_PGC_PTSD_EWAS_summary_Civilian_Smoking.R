########################################################################################
# PGC EWAS Civilian Meta-Analysis: Summarizing Smoking Results
########################################################################################

rm(list=ls())
library(rmeta)
library(ChAMP)
library(forestplot)
library(qqman)

setwd("/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/Civilian/Smoke/")

########################################################################################
# Step 1: Loading Data
########################################################################################

# Combining results files
load("PGC_EWAS_smoking_inverseNorm_intersectSites_civilian.Rdata")
load("PGC_EWAS_smoking_inverseNorm_remainingSites_civilian.Rdata")

all<-pval.combined.marot
rm(pval.combined.marot)
colnames(all)<-c("CpG", "p")
Studies<-rep(c("DNHS, GTP, WTC"), nrow(all))
all<-cbind(all, Studies)
all<-all[, colnames(results)]
all(colnames(all)==colnames(results))
results<-rbind(all, results)

# Removing CpG sites not in at least two studies
dim(results)
results<-results[results[, "Studies"]!="DNHS",]
results<-results[results[, "Studies"]!="GTP",]
results<-results[results[, "Studies"]!="WTC",]
dim(results)

results<-results[order(as.numeric(results[, "p"])), ]
FDR<-p.adjust(as.numeric(results[, "p"]), method="fdr", n=nrow(results))
results<-cbind(results, FDR)
save(results, file="PGC_EWAS_inverseNorm_allResults_civilian_smoking.Rdata")

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

write.csv(tab, file="PGC_EWAS_inverseNorm_allResults_civilian_smoking_summary.csv")
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

png("PGC_EWAS_inverseNorm_QQplot_civilian_smoking.png", width=1280, height=960, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()

png("PGC_EWAS_inverseNorm_QQplot_civilian_smoking_forPublication.png", width=1600, height=800, units="px", bg="white")
par(mar=c(10,10,4,2), mgp=c(6,2,0))
ggd.qqplot(pvalue, main = paste(""))
dev.off()


rm(list=ls()[-match("results", ls())])

########################################################################################
# Step 3: Manhattan plot
########################################################################################

rm(list=ls())

load("PGC_EWAS_inverseNorm_allResults_civilian_smoking.Rdata")
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

cutpoint<--log10((3*0.05)/nrow(results))
df<-data.frame(results)
labs<-c(1:22, "X", "Y")

opar<-par()

png("PGC_EWAS_manhattan_civilian_smoking.png",height=600, width=1200, units="px")
par(mar=c(5.1,6.1, 4.1, 2.1))
manhattan(df, suggestiveline=FALSE, genomewideline=cutpoint, ylim=c(0,8),
          chrlabs=labs, cex=1, cex.axis=1.5, cex.lab=2, col=c("blue4", "red4"))
dev.off()

png("PGC_EWAS_manhattan_civilian_smoking_forPublication.png",height=600, width=1200, units="px")
par(mar=c(5.1,6.1, 4.1, 2.1))
manhattan(df, suggestiveline=FALSE, genomewideline=cutpoint, ylim=c(0,8.5),
          chrlabs=labs, cex=2, cex.axis=2, cex.lab=3, col=c("blue4", "red4"))
dev.off()

postscript("PGC_EWAS_manhattan_civilian_smoking.eps", horizontal=FALSE, width = 6, height = 3, 
           onefile = FALSE, paper = "special", colormodel = "cmyk")
#par(mar=c(5.1,6.1, 4.1, 2.1))
par(opar)
manhattan(df, suggestiveline=FALSE, genomewideline=FALSE, ylim=c(0,8.5),
          chrlabs=labs, cex=0.5, cex.axis=0.5, cex.lab=1)#, col=c("blue4", "red4"))
abline(h=cutpoint, lty=2)
dev.off()

png("PGC_EWAS_manhattan_civilian_smoking_BW.png", width = 600, height = 300, units="px")
#tiff("PGC_EWAS_manhattan_wINTRuST_BW.tiff", width = 4, height = 3, units="in", res=300, 
#     compression = "lzw")
#par(mar=c(5.1,6.1, 4.1, 2.1))
#par(opar)
manhattan(df, suggestiveline=FALSE, genomewideline=FALSE, ylim=c(0,8.5),
          chrlabs=labs, cex=0.5, cex.axis=0.75, cex.lab=1)#, col=c("blue4", "red4"))
abline(h=cutpoint, lty=2)
dev.off()


rm(list=ls())

########################################################################################
# Step 4: Top Results Table
########################################################################################

load("PGC_EWAS_inverseNorm_allResults_civilian_smoking.Rdata")
rownames(results)<-results[, "CpG"]
cpgs<-rownames(results[as.numeric(results[, "p"])<=5*10^-5, ]) # subsetting only most significant results
topResults<-read.csv("/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/nonSmoke/wINTRuST/PGC_EWAS_topresults.csv",
                     stringsAsFactors = F, row.names=1)
topResults<-rownames(topResults)
topResults<-topResults[!topResults%in%cpgs]
results<-results[c(cpgs, topResults), ] # subsetting only most significant results
results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

load("PGC_EWAS_DataPrep_smoking_civilian.Rdata")
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

colnames(WTC.coef)
WTC.coef<-WTC.coef[rownames(WTC.coef)%in%sites,c("PTSD", "N.subjects")]
WTC.results<-WTC.results[rownames(WTC.coef),]
colnames(WTC.coef)<-c("PTSD", "N.subjects")
all(rownames(WTC.coef)==rownames(WTC.results))
WTC<-data.frame(WTC.coef, WTC.results)
str(WTC)
WTC$s.e.<-WTC$PTSD/WTC$t # calculate the standard error
WTC$weight<-1/(WTC$s.e.^2) # calculate weight
rm(WTC.results, WTC.coef, WTC.oneSided)

results$WTC.beta<-results$GTP.beta<-results$DNHS.beta<-NA
results$WTC.se<-results$GTP.se<-results$DNHS.se<-NA
results$WTC.N<-results$GTP.N<-results$DNHS.N<-NA
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
                                   "DNHS.beta", "GTP.beta", "WTC.beta",
                                   "DNHS.se", "GTP.se", "WTC.se",
                                   "DNHS.N", "GTP.N", "WTC.N")]
colnames(probe.features)<-c("Studies", "CpG", "CHR", "Position", "Gene", "Feature",
                            "beta", "variance", "p-value", "FDR",
                            "DNHS.beta", "GTP.beta", "WTC.beta",
                            "DNHS.se", "GTP.se", "WTC.se",
                            "DNHS.N", "GTP.N", "WTC.N")

rownames(probe.features)<-c(1:nrow(probe.features))
write.csv(probe.features, "PGC_EWAS_TopResults_civilian_smoking.csv")
rm(list=ls())

results<-read.csv("PGC_EWAS_TopResults_civilian_smoking.csv", row.names=1, stringsAsFactors=F)
studies<-c("DNHS", "GTP", "WTC")

pdf("PGC_EWAS_Civilian_smoking_forestPlots.pdf")
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
rm(list=ls())

results<-read.csv("PGC_EWAS_TopResults_civilian_smoking.csv", row.names=1, stringsAsFactors=F)

rownames(results)<-results[, "CpG"]
cpgs<-rownames(results[results[, "FDR"]<=0.05,])
topResults<-read.csv("/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/nonSmoke/wINTRuST/PGC_EWAS_topresults.csv",
                     stringsAsFactors = F, row.names=1)
topResults<-rownames(topResults)
topResults<-topResults[!topResults%in%cpgs]

results<-results[c(cpgs, topResults), ] # subsetting only most significant results
results<-data.frame(results, stringsAsFactors=F)
results$p<-as.numeric(results$p)
results$FDR<-as.numeric(results$FDR)
sites<-rownames(results)

studies<-unlist(strsplit(results[1, "Studies"], ","))
studies<-gsub(" ", "", studies)
names<-matrix(nrow=length(studies), ncol=2)
colnames(names)<-c("Old", "Correct")
names[, "Old"]<-studies
rownames(names)<-names[, "Old"]

names[match("VA", rownames(names)), "Correct"]<-"VA-NCPTSD"
names[match("MIR", rownames(names)), "Correct"]<-"VA-M-AA"
names[match("DUKE", rownames(names)), "Correct"]<-"VA-M-EU"
names[match("TRUST", rownames(names)), "Correct"]<-"INTRuST"
names[match("AS", rownames(names)), "Correct"]<-"Army STARRS"
names[which(is.na(names[, "Correct"])), "Correct"]<-names[which(is.na(names[, "Correct"])), "Old"]

names<-names[order(names[, "Correct"]), ]
studies<-names[, "Old"]

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
  stud<-c(results[ii, "CpG"], "Study", names[, "Correct"], NA, "Summary")
  summary<-cbind(stud, beta.table, Ns)
  summary[which(is.na(summary))]<-""
  summary<-summary[, 1]
  
  # Call forestplot
  cochrane<-data.frame(betas, lower, upper)
  summary<-summary[-1]
  cochrane<-cochrane[-1,]
  postscript(paste("PGC_EWAS_civilian_smoking_", results[ii, "CpG"], ".eps", sep=""), horizontal=FALSE, width = 3, height = 3, 
             onefile = FALSE, paper = "special", colormodel = "cmyk")
  
  forestplot(summary, cochrane,
             # boxsize=c(0.8, 0.8, 0.4, 0.6, 0.5), lwd.ci=1.5,
             new_page = TRUE,
             is.summary=c(TRUE, rep(FALSE,3),TRUE),
             xlog=FALSE, #lineheight=unit(1,"cm"),
             col=fpColors(box="grey48",line="grey48", summary="black"),
             #xticks=(c(-0.4, -0.2, 0.1)),
             txt_gp = fpTxtGp( label = list(gpar(fontfamily="", cex=0.5)),
                               ticks = gpar(fontfamily = "", cex=0.5)))
  dev.off()
  
}
