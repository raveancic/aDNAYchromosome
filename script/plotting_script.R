#### Prediction of the main Ychr haplogroup in aDNA samples######
library(ggplot2)
library(dplyr)
library(cowplot)
library(jcolors)
library(scales)

source("script/function.R")

#The input file for the sample you want to test has to be tfile plink format tped ecc ecc.

files_tped <- list.files("data/Panama_chrY.aDNA/",pattern = "*aDNA*.tped")
tped.panama.ancient <- read.table(paste0("data/Panama_chrY.aDNA/", files_tped), header = F, na.strings = "N")

tfam.ancient <- read.table(paste0("data/Panama_chrY.aDNA/", gsub("tped", "tfam", files_tped)), header = F, na.strings = "N")

# Take the file lists with snps
# This is the list taken by Poznick et al., 2016
pos <- read.table("data/all.haplogroups.snps.txt")

# This is the list taken by Grugni&Raveane et al., 2019
pos_spec <- read.table("data/pos.hgQ.haplogroup", header = T, sep = "\t")

# This is the list taken by Pinotti et al., 2018 look at the Supp Tables
pos_spec_Pinotti_M242 <- read.csv("data/1-s2.0-S0960982218314957-mmc3.csv")

# Run the function on all the tped samples

tabs <- tpedapply(tped.panama.ancient, pos, as.character(tfam.ancient$V1))

plots <- list()

for(i in 1:length(tabs)){
  
  tabggplot <- tabs[[i]]
  
  col_bp <- jcolors(palette = "rainbow")[3:6] 
  names(col_bp) <- sort(unique(tabggplot$Damage))
  
  
  plots[[i]] <- ggplot(data=as.data.frame(tabggplot), aes(x=Hg, y=SNP,fill=Damage)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = col_bp) +
    labs(title= names(tabs[i]),
         x= "Y haplogroups (as in Poznick et al., 2016)",
         y= "number of SNPs") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank()
    )
  
  
}


p2save <- plot_grid(plotlist = plots)  

ggsave("figures/samplesplot.png", p2save, width = 21, height = 7)
ggsave("figures/samplesplot.pdf", p2save, width = 21, height = 7)




##### FOR SPECIFIC HAPLOGROUP##################
  
  
  ########The input file for the sample you want to test has to be in plink format tped###########################Not yet implemented in a loop so the analyses have to be sample per sample######################
  
  
  concat1 <- tped.panama.ancient[na.omit(match(pos_spec$GrCh37.hg19, tped.panama.ancient$V4)),c(4,5)]
  concat2  <- pos_spec[match(concat1$V4, pos_spec$GrCh37.hg19),]
  concat <- cbind(concat1, concat2)
  colnames(concat) <- c("posPlink", "Pa16","Ch37","Ch38","Hg","anc","der")
  
  concat$HgYes <- NA
  concat$HgNo <- NA
  concat$Damage <- NA
  
  
  rownames(concat) <- c(1:nrow(concat))
  #### seewp up the x coz they are not consistent as said in Poznik et al. 2016 ######
  
  concat[concat$accordance=="x",] <-   concat[concat$accordance=="x", c(1,2,3,4,5,7,6,8,9)]  
  
  
  for(r in 1:nrow(concat)){
    concat$HgYes[r] <- ifelse(as.character(concat$Pa16[r])==as.character(concat$der[r]), "der", "")
    concat$HgNo[r] <- ifelse(as.character(concat$Pa16[r])==as.character(concat$anc[r]), "anc", "")
    concat$Damage[r] <- ifelse(as.character(concat$anc[r])=="C" & as.character(concat$der[r])=="T", "dam", "")
    concat$Damage[r] <- ifelse(as.character(concat$anc[r])=="G" & as.character(concat$der[r])=="A", "dam", "")
    
  }
  #### Compose two dataset divided by HG yes and Hg no #######
  
  concat_4plotHgYes.tmp1 <-   concat %>% group_by(Hg, HgYes, Damage) %>% 
    summarise(n = n())
  concat_4plotHgYes <- concat_4plotHgYes.tmp1[concat_4plotHgYes.tmp1$HgYes=="der",]
  
  
  concat_4plotHgNo.temp1 <-   concat %>% group_by(Hg, HgNo, Damage) %>% 
    summarise(n = n())
  concat_4plotHgNo <- concat_4plotHgNo.temp1[concat_4plotHgNo.temp1$HgNo=="anc",]
  
  ### Make the barplot#####
  concat_4plotHgYes$Damage[concat_4plotHgYes$Damage=="dam"] <- c("der_dam")
  concat_4plotHgYes$Damage[concat_4plotHgYes$Damage==""] <- c("der")
  concat_4plotHgNo$Damage[concat_4plotHgNo$Damage=="dam"] <- c("anc_dam")
  concat_4plotHgNo$Damage[concat_4plotHgNo$Damage==""] <- c("anc")
  
  colnames(concat_4plotHgYes) <- c("Hg", "HgPred", "Damage", "snps")
  colnames(concat_4plotHgNo) <-  c("Hg", "HgPred", "Damage", "snps")
  concat_4plotHgFinal <- rbind(as.data.frame(concat_4plotHgYes),as.data.frame(concat_4plotHgNo))
  
  p<-ggplot(data=as.data.frame(concat_4plotHgFinal), aes(x=Hg, y=snps,fill=Damage)) +
    geom_bar(stat="identity")
  p <- p + scale_fill_brewer(palette="Paired", direction = -1) + labs(title= gsub(".Y.final.bam.30.vcf.gz.plink.tped", "",files_tped[i])+ scale_fill_brewer(palette="Paired", direction = -1) ) + theme(axis.text.x = element_text(size = rel(0.8),angle = 90, hjust = 1))
  
  assign(paste0(gsub(".Y.final.bam.30.vcf.gz.plink.tped", "",files_tped[i]),"_spec"), p)
  
  
}



pdf("predictor_panama_antico_generalHg.pdf", width = 13, height = 11)

plot_grid(Pan1102_gen,Pa09_gen,Pa10_gen,Pa16A_gen,Pa24A_gen,Pa30A_gen,labels= c('A','B','C','D','E','F' ), ncol = 2)

dev.off()

pdf("predictor_panama_antico_specificHg.pdf", width = 13, height = 11)

plot_grid(Pa09_spec,Pa10_spec,Pa16A_spec,Pa24A_spec,Pa30A_spec,labels= c('A','B','C','D','E','F'), ncol = 2)

dev.off()


###### Plot them #######

