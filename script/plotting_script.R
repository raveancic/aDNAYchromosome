#### Prediction of the main Ychr haplogroup in aDNA samples
library(ggplot2)
library(dplyr)
library(cowplot)
library(jcolors)
library(scales)


source("script/function.R")

#The input file for the sample you want to test has to be tfile plink format tped ecc ecc.

# This input file is in the folder data/Panama_chrY.aDNA/Panama_chrY.aDNA.tar.gz
system("tar -xvf data/Panama_chrY.aDNA/Panama_chrY.aDNA.tar.gz")

files_tped <- list.files("data/Panama_chrY.aDNA/",pattern = "*aDNA.tped")
tped.panama.ancient <- read.table(paste0("data/Panama_chrY.aDNA/", files_tped), header = F, na.strings = "N")

# tfam.ancient <- read.table(paste0("data/Panama_chrY.aDNA/", gsub("tped", "tfam", files_tped)), header = F, na.strings = "N")
tfam.ancient <- read.table("data/Panama_chrY.aDNA/new_namesaDNA.txt", header = F, na.strings = "N")

tfam.ancient$sex <- rep("M", nrow(tfam.ancient))
tfam.ancient$sex[tfam.ancient$V4 == "PAPV172"] <- "F"

# Give names to the tped columns

colnames_tped <- paste(as.character(unlist(lapply(tfam.ancient$V4, function(x) rep(x,2)))), 
      rep(c("A","B"), times=length(tfam.ancient$V1)), sep = "_")

colnames(tped.panama.ancient) <- c("CHR", "rs", "CM", "pos", colnames_tped)
 

# select only the M
tfam.ancient_2test <- tfam.ancient[tfam.ancient$sex=="M",]


# Take the file lists with snps
# This is the list taken by Poznick et al., 2016
pos <- read.table("data/all.haplogroups.snps.txt")

# This is the list taken by Grugni&Raveane et al., 2019
pos_spec <- read.table("data/pos.hgQ.haplogroup", header = T, sep = "\t")
# Little modification as table before to make it attractive for the function
#1 Hg #2Branch #3Pos #4Anc #5Der #6aff
pos_spec_mod <- data.frame("Hg"= pos_spec$NOMI.GRUPPI, 
                           "Branch" = NA,
                           "Pos"= pos_spec$GrCh37.hg19,
                           "Anc"= pos_spec$ANC,
                           "Der"= pos_spec$DER,
                           "acc"= ".")

# # This is the list taken by Pinotti et al., 2018
# pos_spec_Pinotti_M242 <- read.table("data/1-s2.0-S0960982218314957-mmc3.csv", header = T, sep = ",")
# # Little modification as table before to make it attractive for the function
# #1 Hg #2Branch #3Pos #4Anc #5Der #6aff
# pos_spec_Pin_mod <- data.frame("Hg"= pos_spec_Pinotti_M242$SNP_eq, 
#                            "Branch" = NA,
#                            "Pos"= pos_spec_Pinotti_M242$position_GRCh37,
#                            "Anc"= pos_spec_Pinotti_M242$reference,
#                            "Der"= pos_spec_Pinotti_M242$alternate,
#                            "acc"= gsub(".x.", "x" ,gsub("", ".", gsub("\\*", "x", pos_spec_Pinotti_M242$X.note.2.))))

# Run the function on all the tped samples

# Define the colors 

col_bp <- jcolors(palette = "rainbow")[3:6] 
names(col_bp) <- c("anc_dam", "anc", "der_dam", "der")

lab <- c("anc.", "anc. damage", "der.","der. damage")
names(lab) <- c("anc", "anc_dam", "der", "der_dam")

# Run the function

tab <- list()

tabs <- tpedapply(tped.panama.ancient, CrTab_DamageHg, pos, as.character(tfam.ancient_2test$V4))


plots <- list()

for(i in 1:length(tabs)){
  
  tabggplot <- tabs[[i]]
  
  
  
  plots[[i]] <- ggplot(data=as.data.frame(tabggplot), aes(x=Hg, y=SNP,fill=Damage)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = col_bp
                      , labels=lab
                      ) +
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


n = length(plots)/2

p2save_gen <- plot_grid(plotlist = plots, nrow = n, labels = "AUTO")  

ggsave("figures/samplesplot_genv2.png", p2save_gen, width = 15, height = 13)
ggsave("figures/samplesplot_genv2.pdf", p2save_gen, width = 15, height = 13)
# ggsave("~/Desktop/samplesplot_gen.pdf", p2save, width = 21, height = 10)

# For specific Hg Tomorrow with our position.
  
# Run the function on all the tped samples

tabs_spec <- list()

tabs_spec <- tpedapply(tped.panama.ancient, CrTab_DamageHg, pos_spec_mod, as.character(tfam.ancient$V4))

tab_Q <- list()

tab_Q <- tabs_spec[unlist(lapply(tabs_spec, function(x) any(x[,2]=="der")))]

plots <- list()

for(i in 1:length(tab_Q)){
  
  tabggplot <- tab_Q[[i]]
  # 
  # col_bp <- jcolors(palette = "rainbow")[3:6] 
  # names(col_bp) <- sort(unique(tabggplot$Damage))
  # 
  
  plots[[i]] <- ggplot(data=as.data.frame(tabggplot), aes(x=Hg, y=SNP,fill=Damage)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(values = col_bp, 
                      labels= lab
                      ) +
    labs(title= names(tab_Q[i]),
         x= "Spec Y haplogroups (as in Grugni et al., 2019)",
         y= "number of SNPs") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank()
    )
  
  
}

p2save_spec <- plot_grid(plotlist = plots, nrow = 2, labels= "AUTO")  

ggsave("figures/samplesplot_Specificv2.png", p2save_spec, width = 15, height = 8)
ggsave("figures/samplesplot_Specificv2.pdf", p2save_spec, width = 15, height = 8)
ggsave("~/Desktop/samplesplot_Grugni2019.pdf", p2save, width = 21, height = 10)
