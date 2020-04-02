#library
library(dplyr)
# Create pthe first function 

CrTab_DamageHg <- function(tpedtab, posfile){
  
  # It match the position in hg37 with the tped
  concat1 <- tpedtab[na.omit(match(posfile[,3], tpedtab[,4])),c(4,5)]
  # It matches back to the previous file for the Hg
  concat2  <- posfile[match(concat1$V4, posfile$V3),]
  # It generates the table for the prediction
  concat <- cbind(concat1, concat2)
  #Give the colname
  colnames(concat) <- c("posPlink", "sample","Hg","BranchNumb","posPoznik","anc","der","accordance")
  concat$HgYes <- NA
  concat$HgNo <- NA
  concat$Damage <- NA
  
  rownames(concat) <- c(1:nrow(concat))
  #### seewp up the x coz they are not consistent as said in Poznik et al. 2016 ######
  
  concat <- concat[concat$accordance!="x",]
  
  # return(concat)
  
  #
  for(r in 1:nrow(concat)){
    
    # See if we have the hg
    concat$HgYes[r] <- ifelse(as.character(concat$sample[r])==as.character(concat$der[r]), "der", "")
    concat$HgNo[r] <- ifelse(as.character(concat$sample[r])==as.character(concat$anc[r]), "anc", "")
    
    # Look at the damage C-->T
    concat$Damage[r] <- ifelse(as.character(concat$anc[r])=="C" & as.character(concat$der[r])=="T", "dam", "")
    # Look at the damage G -->A
    concat$Damage[r] <- ifelse(as.character(concat$anc[r])=="G" & as.character(concat$der[r])=="A", "dam", "") }
  
  # Two data set divided by hg YES and hg NO
  concat_4plotHgYes.tmp1 <-   concat %>% group_by(Hg, HgYes, Damage) %>% 
    summarise(n = n())
  concat_4plotHgYes <- concat_4plotHgYes.tmp1[concat_4plotHgYes.tmp1$HgYes=="der",]
  
  concat_4plotHgNo.temp1 <-   concat %>% group_by(Hg, HgNo, Damage) %>% 
    summarise(n = n())
  concat_4plotHgNo <- concat_4plotHgNo.temp1[concat_4plotHgNo.temp1$HgNo=="anc",]
  
  # Prepare the 
  concat_4plotHgYes$Damage[concat_4plotHgYes$Damage=="dam"] <- c("der_dam")
  concat_4plotHgYes$Damage[concat_4plotHgYes$Damage==""] <- c("der")
  concat_4plotHgNo$Damage[concat_4plotHgNo$Damage=="dam"] <- c("anc_dam")
  concat_4plotHgNo$Damage[concat_4plotHgNo$Damage==""] <- c("anc")
  
  colnames(concat_4plotHgYes) <- c("Hg", "HgPred", "Damage", "SNP")
  colnames(concat_4plotHgNo) <-  c("Hg", "HgPred", "Damage", "SNP")
  concat_4plotHgFinal <- rbind(as.data.frame(concat_4plotHgYes),as.data.frame(concat_4plotHgNo))
  
  return(concat_4plotHgFinal[complete.cases(concat_4plotHgFinal),])
  
}


tpedapply <- function(tpedtab, p, namesamples){
  
  # tpedtab= tped.panama.ancient
  
  strcln = 6:ncol(tped.panama.ancient)
  
  clm2take = strcln[strcln %% 2 == 0]
  
  myl <- list()
  
  for(i in 1:length(clm2take)){
      
      tab = cbind(tpedtab[1:4],tpedtab[clm2take[i]-1:clm2take[i]])
      
      myl[[i]] <- CrTab_DamageHg(tab, pos)
    
  }
  
  names(myl) <- namesamples
  return(myl)
}









