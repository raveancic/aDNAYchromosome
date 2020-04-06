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

tabs <- tpedapply(tped.panama.ancient, CrTab_DamageHg, pos, as.character(tfam.ancient$V4))


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

<<<<<<< HEAD
=======

>>>>>>> 2d07c55f3c615dcc93b6fd6ee012b9ac4495b7b2
p2save_gen <- plot_grid(plotlist = plots)  

ggsave("figures/samplesplot_gen.png", p2save, width = 21, height = 10)
ggsave("figures/samplesplot_gen.pdf", p2save, width = 21, height = 10)
<<<<<<< HEAD
=======
# ggsave("~/Desktop/samplesplot_gen.pdf", p2save, width = 21, height = 10)
>>>>>>> 2d07c55f3c615dcc93b6fd6ee012b9ac4495b7b2

# For specific Hg Tomorrow with our position.
  
# Run the function on all the tped samples

tabs <- tpedapply(tped.panama.ancient, CrTab_DamageHg, pos_spec_mod, as.character(tfam.ancient$V4))

tab_Q <- tabs[unlist(lapply(tabs, function(x) any(x[,2]=="der")))]

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
    labs(title= names(tabs[i]),
         x= "Spec Y haplogroups (as in Grugni et al., 2019)",
         y= "number of SNPs") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank()
    )
  
  
}

p2save <- plot_grid(plotlist = plots)  

ggsave("figures/samplesplot_Specific.png", p2save, width = 21, height = 10)
ggsave("figures/samplesplot_Specific.pdf", p2save, width = 21, height = 10)
ggsave("~/Desktop/samplesplot_Grugni2019.pdf", p2save, width = 21, height = 10)
