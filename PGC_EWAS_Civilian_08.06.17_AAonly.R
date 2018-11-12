w########################################################################################################
# PGC EWAS: Civilan- and AA-only Results
########################################################################################################

rm(list=ls())
setwd("/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/08.06.17/Civilian/AA-only/")
dataDir<-"/Users/ly2207/Documents/Andrew/R/PGC/Current/PGC_EWAS/08.06.17/Civilian/AA-only/"
studies<-c("DNHS_AA-only", "GTP_AA-only")
mods<-c("main_results", "main_results_smokersOnly", "main_results_nonSmokersOnly",
  "main_results_malesOnly", "main_results_femalesOnly")

for(jj in 1:length(mods)){
  m<-mods[jj]
  for(ll in 1:length(studies)){
    # Load Data
    load(paste(dataDir, studies[ll], "/", studies[ll], "_", m, ".Rdata", sep=""))
    oneSided<-pt(fit2.ebayes$t, df = (fit2.ebayes$df.prior+fit2.ebayes$df.residual))
    assign(paste(studies[ll], ".fit.coef", sep=""), fit.coef)
    assign(paste(studies[ll], ".results", sep=""), results)
    assign(paste(studies[ll], ".fit2.ebayes", sep=""), fit2.ebayes)
    assign(paste(studies[ll], ".oneSided", sep=""), oneSided)
    rm(fit.coef, results, fit2.ebayes)
  }
  
  tab<-matrix(nrow=length(studies), ncol=9)
  colnames(tab)<-c("Study", "Sites", "Min N", "Max N", "Lambda", "FDR", "p<5x10^-5", "p<5x10^-6", "p<5x10^-7")
  
  # Summary Table
  for(ii in 1:length(studies)){
    tab[ii, "Study"]<-studies[ii]
    res<-get(paste(studies[ii], ".results", sep=""))
    res<-res[!is.na(res[, "logFC"]),] # this is for VA Boston. Ther are no p-values for some sites. I'm guess there weren't any cases
    tab[ii, "p<5x10^-5"]<-sum(res[, "P.Value"]<(5*10^(-5))) 
    tab[ii, "p<5x10^-6"]<-sum(res[, "P.Value"]<(5*10^(-6))) 
    tab[ii, "p<5x10^-7"]<-sum(res[, "P.Value"]<(5*10^(-7))) 
    tab[ii, "FDR"]<-sum(res[, "adj.P.Val"]<=0.05)
    tab[ii, "Sites"]<-nrow(res) 
    coef<-get(paste(studies[ii], ".fit.coef", sep=""))
    tab[ii, "Min N"]<-range(coef[, "N.subjects"])[1] 
    tab[ii, "Max N"]<-range(coef[, "N.subjects"])[2] 
    p<-res[, "P.Value"]
    chisq<-qchisq(1-p,1)
    tab[ii, "Lambda"]<-round(median(chisq)/qchisq(0.5,1),3)
    rm(p, chisq, coef, res)
  }
  write.csv(tab, paste(paste(studies, collapse="."), m, "summary.csv", sep="_"))
  
  ### Meta-Analysis ###
  
  ## Intersecting Sites ##
  for(kk in 1:length(studies)){
    if(kk==1){
      temp<-get(paste(studies[kk], ".fit.coef", sep=""))
      sites<-rownames(temp)
      all<-rownames(temp)
      rm(temp)
    }else{
      temp<-get(paste(studies[kk], ".fit.coef", sep=""))
      sites<-intersect(sites, rownames(temp))
      all<-append(all, rownames(temp))
      rm(temp)
    }
  }
  all<-unique(all)
  
  # t-statistics from 1 sided p-values
  tstats<-matrix(nrow=length(sites), ncol=length(studies))
  colnames(tstats)<-studies
  rownames(tstats)<-sites
  
  for(nn in 1:ncol(tstats)){
    temp<-get(paste(studies[nn], ".oneSided", sep=""))
    if(!all(sites%in%rownames(temp))){
      print("t-stats don't match")
      break
    }
    temp<-temp[sites, ]
    tstats[, studies[nn]]<-qnorm(1-as.numeric(temp)) # t-statistic from 1-sided p-value
    rm(temp)
  }
  
  # Weights
  weights<-matrix(nrow=length(sites), ncol=length(studies))
  colnames(weights)<-studies
  rownames(weights)<-sites
  
  for(mm in 1:length(studies)){
    temp<-get(paste(studies[mm], ".fit.coef", sep=""))
    if(!all(sites%in%rownames(temp))){
      print("Weights don't match")
      break
    }
    temp<-temp[sites, ]
    weights[, studies[mm]]<-as.numeric(temp[, "N.subjects"])
    rm(temp)
  }
  
  for(ww in 1:nrow(weights)){
    weights[ww, ]<-sqrt(weights[ww, ]/sum(weights[ww, ]))
  }
  
  # Calculated Z-score and p-value
  w.tstats<-weights*tstats
  z.combined<-apply(w.tstats, 1, sum)
  p<-2*(1-pnorm(abs(z.combined)))
  res<-cbind(z.combined, p)
  combinedStudies<-rep(paste(studies, collapse="_"), times=nrow(res))
  res<-cbind(res, combinedStudies)
  cpgs<-rownames(res)
  res<-cbind(res, cpgs)
  res<-res[, c("cpgs", "combinedStudies", "p")]
  colnames(res)<-c("CpG", "Studies", "p")
  rm(combinedStudies, cpgs)
  
  ## Remaining Sites ##
  sites<-all[!all%in%sites]
  sTemp<-matrix(FALSE, nrow=length(sites), ncol=length(studies))
  colnames(sTemp)<-studies
  rownames(sTemp)<-sites
  
  for(pp in 1:length(studies)){
    results<-get(paste(studies[pp], ".results", sep=""))
    sTemp[rownames(sTemp)%in%rownames(results), studies[pp]]<-TRUE
    miss<-names(which(sTemp[, studies[pp]]==FALSE)) 
    unmat<-rownames(sTemp)[is.na(match(rownames(sTemp), rownames(results)))]
    
    if(sum(is.na(match(miss, unmat)))!=0){
      print(paste("Error in Remaining Sites: ", studies[pp]), sep="")
    }
    rm(miss, unmat, results)
  }
  
  # Step 2C: Loop to calculate p-values
  results<-matrix(nrow=nrow(sTemp), ncol=3)
  colnames(results)<-c("CpG", "Studies", "p")
  z.combined<-NULL
  
  for(qq in 1:nrow(sTemp)){
    cpg<-rownames(sTemp)[qq]  # CpG site of analysis
    inc<-colnames(sTemp)[sTemp[qq,]] # sTemp included in analysis
    results[qq, "CpG"]<-cpg
    results[qq, "Studies"]<-paste(inc, collapse=", ")
    
    # Table of parameters for the inverse normal method
    tab<-matrix(nrow=length(inc), ncol=4)
    rownames(tab)<-inc
    colnames(tab)<-c("N", "Weight", "Zi", "Z")
    
    # Getting the number of subjects and calculating weights
    for(ww in 1:length(inc)){
      temp<-get(paste(inc[ww], ".fit.coef", sep=""))
      tab[inc[ww], "N"]<-temp[cpg, "N.subjects"]
      rm(temp)
    }
    tab[, "Weight"]<-sqrt(tab[, "N"]/sum(tab[, "N"]))
    
    # Getting the Z and weighted Z values for each included study
    
    for(vv in 1:length(inc)){
      temp<-get(paste(inc[vv],".oneSided", sep=""))
      tab[inc[vv], "Zi"]<-qnorm(1-as.numeric(temp[cpg,]))
      rm(temp)
    }
    tab[, "Z"]<-tab[, "Zi"]*tab[, "Weight"]
    
    # Calculating effect and two-sided p-value
    Sg<-sum(tab[, "Z"]) # add sTemp together
    z.combined<-append(z.combined, Sg)
    results[qq, "p"]<-2*(1-pnorm(abs(Sg))) # calculated two-tailed p-value
    rm(tab)
  }
  
  rownames(results)<-results[, "CpG"]
  if(!all(colnames(results)==colnames(res))){
    print("Error Combining results")
    break
  }
  
  res<-rbind(res, results)
  save(res, file=paste(paste(studies, collapse="_"), "_", m, ".Rdata", sep=""))
  print(jj)
}

