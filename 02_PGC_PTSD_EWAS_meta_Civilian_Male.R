########################################################################################
# PGC Meta-Analysis: CpG Sites in All Studies
########################################################################################

rm(list=ls())

########################################################################################
# Step 1: Load and subset data
########################################################################################

load("PGC_EWAS_combinedData_civilian_males.Rdata")

# The minimum adjusted p-values are:
min(DNHS.results[,"adj.P.Val"]) #  
min(GTP.results[,"adj.P.Val"]) # 
min(WTC.results[,"adj.P.Val"]) #  

# Step 1A: Convert data.frames to matrices. This speeds up the later loops
str(DNHS.results) 
DNHS.results.m<-as.matrix(DNHS.results)
all(DNHS.results==DNHS.results.m)
DNHS.results<-DNHS.results.m
rm(DNHS.results.m)

str(GTP.results)
GTP.results.m<-as.matrix(GTP.results)
all(GTP.results==GTP.results.m)
GTP.results<-GTP.results.m
rm(GTP.results.m)

str(WTC.results)
WTC.results.m<-as.matrix(WTC.results)
all(WTC.results==WTC.results.m)
WTC.results<-WTC.results.m
rm(WTC.results.m)

# Step 1B: Check that all the rownames 
sum(is.na(match(rownames(DNHS.coef), rownames(DNHS.results))))
sum(is.na(match(rownames(DNHS.coef), rownames(DNHS.ebayes))))
head(DNHS.results)
DNHS.results<-DNHS.results[rownames(DNHS.ebayes),]
head(DNHS.results)
head(DNHS.coef)
DNHS.coef<-DNHS.coef[rownames(DNHS.ebayes),]
head(DNHS.coef)
all(rownames(DNHS.results)==rownames(DNHS.coef))
all(rownames(DNHS.results)==rownames(DNHS.ebayes))

sum(is.na(match(rownames(GTP.coef), rownames(GTP.results))))
sum(is.na(match(rownames(GTP.coef), rownames(GTP.ebayes))))
head(GTP.results)
GTP.results<-GTP.results[rownames(GTP.ebayes),]
head(GTP.results)
head(GTP.coef)
GTP.coef<-GTP.coef[rownames(GTP.ebayes),]
head(GTP.coef)
all(rownames(GTP.results)==rownames(GTP.coef))
all(rownames(GTP.results)==rownames(GTP.ebayes))

sum(is.na(match(rownames(WTC.coef), rownames(WTC.results))))
sum(is.na(match(rownames(WTC.coef), rownames(WTC.ebayes))))
head(WTC.results)
WTC.results<-WTC.results[rownames(WTC.ebayes),]
head(WTC.results)
head(WTC.coef)
WTC.coef<-WTC.coef[rownames(WTC.ebayes),]
head(WTC.coef)
all(rownames(WTC.results)==rownames(WTC.coef))
all(rownames(WTC.results)==rownames(WTC.ebayes))

# Step 1C: Determine the number of probes available in all studies
dnhs.sites<-rownames(DNHS.coef)
gtp.sites<-rownames(GTP.coef)
wtc.sites<-rownames(WTC.coef)

sites<-append(dnhs.sites, gtp.sites)
sites<-append(sites, wtc.sites)
length(sites) 

all<-intersect(dnhs.sites, gtp.sites)
all<-intersect(all, wtc.sites)

sum(is.na(match(all, rownames(DNHS.coef)))) # these should all be 0 
sum(is.na(match(all, rownames(GTP.coef))))
sum(is.na(match(all, rownames(WTC.coef))))
length(all) 

rm(dnhs.sites, gtp.sites, wtc.sites)

########################################################################################
# Step 2: Data Prep
########################################################################################

# Step 2A: Calculate one-sided p-values for each CpG site's t-statistic
DNHS.oneSided<-pt(DNHS.ebayes$t, df = (DNHS.ebayes$df.prior+DNHS.ebayes$df.residual))
all(rownames(DNHS.oneSided)==rownames(DNHS.results))
all(rownames(DNHS.oneSided)==rownames(DNHS.coef))
all(rownames(DNHS.oneSided)==rownames(DNHS.ebayes))

GTP.oneSided<-pt(GTP.ebayes$t, df = (GTP.ebayes$df.prior+GTP.ebayes$df.residual))
all(rownames(GTP.oneSided)==rownames(GTP.results))
all(rownames(GTP.oneSided)==rownames(GTP.coef))
all(rownames(GTP.oneSided)==rownames(GTP.ebayes))

WTC.oneSided<-pt(WTC.ebayes$t, df = (WTC.ebayes$df.prior+WTC.ebayes$df.residual))
all(rownames(WTC.oneSided)==rownames(WTC.results))
all(rownames(WTC.oneSided)==rownames(WTC.coef))
all(rownames(WTC.oneSided)==rownames(WTC.ebayes))

save.image("PGC_EWAS_DataPrep_nonSmoke_civilian_males.Rdata")

# Step 2B: Subset all and non-all sites

studies<-c("DNHS", "GTP", "WTC")

for(ii in 1:length(studies)){
  coef<-get(paste(studies[ii], ".coef", sep=""))
  coef<-coef[all, ]
  assign(paste(studies[ii], ".coef", sep=""), coef)
  rm(coef)
  results<-get(paste(studies[ii], ".results", sep=""))
  results<-results[all, ]
  assign(paste(studies[ii], ".results", sep=""), results)
  rm(results)
  ebayes<-get(paste(studies[ii], ".ebayes", sep=""))
  ebayes<-ebayes[all, ]
  assign(paste(studies[ii], ".ebayes", sep=""), ebayes)
  rm(ebayes)
  oneSided<-get(paste(studies[ii], ".oneSided", sep=""))
  oneSided<-oneSided[all, ]
  assign(paste(studies[ii], ".oneSided", sep=""), oneSided)
  rm(oneSided)
}

save.image("PGC_EWAS_DataPrep_nonSmoke_intersectSites_civilian_males.Rdata")

rm(list=ls())

load("PGC_EWAS_DataPrep_nonSmoke_civilian_males.Rdata")

studies<-c("DNHS", "GTP", "WTC")

for(ii in 1:length(studies)){
  coef<-get(paste(studies[ii], ".coef", sep=""))
  coef<-coef[!rownames(coef)%in%all,]
  assign(paste(studies[ii], ".coef", sep=""), coef)
  rm(coef)
  results<-get(paste(studies[ii], ".results", sep=""))
  results<-results[!rownames(results)%in%all,]
  assign(paste(studies[ii], ".results", sep=""), results)
  rm(results)
  ebayes<-get(paste(studies[ii], ".ebayes", sep=""))
  ebayes<-ebayes[!rownames(ebayes)%in%all,]
  assign(paste(studies[ii], ".ebayes", sep=""), ebayes)
  rm(ebayes)
  oneSided<-get(paste(studies[ii], ".oneSided", sep=""))
  oneSided<-oneSided[!rownames(oneSided)%in%all, ]
  assign(paste(studies[ii], ".oneSided", sep=""), oneSided)
  rm(oneSided)
}

save.image("PGC_EWAS_DataPrep_nonSmoke_remainingSites_civilian_males.Rdata")

rm(list=ls())

########################################################################################
# Step 3: Meta-Analysis of Sites in All Studies
########################################################################################

load("PGC_EWAS_DataPrep_nonSmoke_intersectSites_civilian_males.Rdata")

# Step 3A: Check that all rownames line up
all(rownames(DNHS.coef)==rownames(DNHS.results))
all(rownames(DNHS.coef)==rownames(DNHS.ebayes))
all(rownames(DNHS.coef)==names(DNHS.oneSided))

all(rownames(DNHS.coef)==rownames(GTP.coef))
all(rownames(DNHS.coef)==rownames(GTP.results))
all(rownames(DNHS.coef)==rownames(GTP.ebayes))
all(rownames(DNHS.coef)==names(GTP.oneSided))

all(rownames(DNHS.coef)==rownames(WTC.coef))
all(rownames(DNHS.coef)==rownames(WTC.results))
all(rownames(DNHS.coef)==rownames(WTC.ebayes))
all(rownames(DNHS.coef)==names(WTC.oneSided))

pval.combined.marot<-NULL

for(ii in 1:nrow(DNHS.coef)){
  probe.subjs<-c(DNHS.coef[ii, "N.subjects"], GTP.coef[ii, "N.subjects"], 
                 WTC.coef[ii, "N.subjects"]) 
  weight<-sqrt(probe.subjs/sum(probe.subjs))
  DNHS.Zi<-qnorm(1-as.numeric(DNHS.oneSided[ii]));DNHS.Zi
  DNHS.Z<-DNHS.Zi*weight[1]; DNHS.Z # multiply by weight 1
  GTP.Zi<-qnorm(1-as.numeric(GTP.oneSided[ii]))
  GTP.Z<-GTP.Zi*weight[2] # multiply by weight 2
  WTC.Zi<-qnorm(1-as.numeric(WTC.oneSided[ii]))
  WTC.Z<-WTC.Zi*weight[3] # multiply by weight 3
  Sg<-DNHS.Z+GTP.Z+WTC.Z # add studies together
  p<-2*(1-pnorm(abs(Sg))) # calculated two-tailed p-value
  pval.combined.marot<-append(pval.combined.marot, p)
  if(ii%%10000==0){
    print(ii)
  }
}

pval.combined.marot<-cbind(rownames(DNHS.results), pval.combined.marot)

save(pval.combined.marot, file="PGC_EWAS_nonSmoke_inverseNorm_intersectSites_civilian_males.Rdata")

rm(list=ls())

########################################################################################
# Step 4: Meta-analysis in remaining sties
########################################################################################

load("PGC_EWAS_DataPrep_nonSmoke_remainingSites_civilian_males.Rdata")
sites<-unique(sites)
sites<-sites[!sites%in%all]

# Step 2B: Add to list of sites the Study IDs
studies<-matrix(FALSE, nrow=length(sites), ncol=3)
colnames(studies)<-c("DNHS", "GTP", "WTC")
rownames(studies)<-sites

# DNHS
studies[rownames(studies)%in%rownames(DNHS.results), "DNHS"]<-TRUE
miss<-names(which(studies[, "DNHS"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(DNHS.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# GTP
studies[rownames(studies)%in%rownames(GTP.results), "GTP"]<-TRUE
miss<-names(which(studies[, "GTP"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(GTP.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

# WTC
studies[rownames(studies)%in%rownames(WTC.results), "WTC"]<-TRUE
miss<-names(which(studies[, "WTC"]==FALSE)) 
unmat<-rownames(studies)[is.na(match(rownames(studies), rownames(WTC.results)))]
sum(is.na(match(miss, unmat))) # should be 0
rm(miss, unmat)

head(studies)
all(apply(studies, 1, sum)>0) # should be true

(dim(studies)[1]*dim(studies)[2])-sum(studies) # this is the number of missing values overall

sum(apply(studies, 1, function(x) 3-sum(x))>0) # this is the number of CpG sites with less than 3 studies
sum(apply(studies, 1, function(x) 2-sum(x))>0) # this is the number of CpG sites with less than 2 studies
sum(apply(studies, 1, function(x) 1-sum(x))>0) # this is the number of CpG sites with less than 1 studies, should be 0

# Step 2C: Loop to calculate p-values
results<-matrix(nrow=nrow(studies), ncol=3)
colnames(results)<-c("CpG", "Studies", "p")

start<-proc.time()[3]
for(ii in 1:nrow(studies)){
  cpg<-rownames(studies)[ii]  # CpG site of analysis
  inc<-colnames(studies)[studies[ii,]] # studies included in analysis
  results[ii, "CpG"]<-cpg
  results[ii, "Studies"]<-paste(inc, collapse=", ")
  
  # Table of parameters for the inverse normal method
  tab<-matrix(nrow=length(inc), ncol=4)
  rownames(tab)<-inc
  colnames(tab)<-c("N", "Weight", "Zi", "Z")
  
  # Getting the number of subjects and calculating weights
  for(jj in 1:length(inc)){
    temp<-get(paste(inc[jj], ".coef", sep=""))
    tab[inc[jj], "N"]<-temp[cpg, "N.subjects"]
    rm(temp)
  }
  tab[, "Weight"]<-sqrt(tab[, "N"]/sum(tab[, "N"]))
  
  # Getting the Z and weighted Z values for each included study
  
  for(ll in 1:length(inc)){
    temp<-get(paste(inc[ll],".oneSided", sep=""))
    tab[inc[ll], "Zi"]<-qnorm(1-as.numeric(temp[cpg]))
    rm(temp)
  }
  tab[, "Z"]<-tab[, "Zi"]*tab[, "Weight"]
  
  # Calculating effect and two-sided p-value
  Sg<-sum(tab[, "Z"]) # add studies together
  results[ii, "p"]<-2*(1-pnorm(abs(Sg))) # calculated two-tailed p-value
  rm(tab)
  
  if(ii%%1000==0){
    print(paste(ii, proc.time()[3]-start, sep=": "))
  }
}

head(results) # check that things look good

save(results, file="PGC_EWAS_nonSmoke_inverseNorm_remainingSites_civilian_males.Rdata")
rm(list=ls())



