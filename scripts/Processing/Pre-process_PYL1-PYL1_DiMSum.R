# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#############################################
## Pre-processing of PYL1-PYL1 DiMSum data ##
#############################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input DiMSum tables ##
############################

## 0 µM ABA

### block 1
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/block1_0uM_ABA_fitness_replicates.RData")
block1.0uM.ABA <- rbind(all_variants, synonymous)
block1.0uM.ABA <- block1.0uM.ABA[-which(block1.0uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 2
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/block2_0uM_ABA_fitness_replicates.RData")
block2.0uM.ABA <- rbind(all_variants, synonymous)
block2.0uM.ABA <- block2.0uM.ABA[-which(block2.0uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 3
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/block3_0uM_ABA_fitness_replicates.RData")
block3.0uM.ABA <- rbind(all_variants, synonymous)
block3.0uM.ABA <- block3.0uM.ABA[-which(block3.0uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 4
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/block4_0uM_ABA_fitness_replicates.RData")
block4.0uM.ABA <- rbind(all_variants, synonymous)
block4.0uM.ABA <- block4.0uM.ABA[-which(block4.0uM.ABA$WT == T)[-1],] ## remove duplicated WT

## 250 µM ABA

### block 1
load("../../data/DiMSum/PYL1-PYL1/250uM_ABA/block1_250uM_ABA_fitness_replicates.RData")
block1.250uM.ABA <- rbind(all_variants, synonymous)
block1.250uM.ABA <- block1.250uM.ABA[-which(block1.250uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 2
load("../../data/DiMSum/PYL1-PYL1/250uM_ABA/block2_250uM_ABA_fitness_replicates.RData")
block2.250uM.ABA <- rbind(all_variants, synonymous)
block2.250uM.ABA <- block2.250uM.ABA[-which(block2.250uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 3
load("../../data/DiMSum/PYL1-PYL1/250uM_ABA/block3_250uM_ABA_fitness_replicates.RData")
block3.250uM.ABA <- rbind(all_variants, synonymous)
block3.250uM.ABA <- block3.250uM.ABA[-which(block3.250uM.ABA$WT == T)[-1],] ## remove duplicated WT

### block 4
load("../../data/DiMSum/PYL1-PYL1/250uM_ABA/block4_250uM_ABA_fitness_replicates.RData")
block4.250uM.ABA <- rbind(all_variants, synonymous)
block4.250uM.ABA <- block4.250uM.ABA[-which(block4.250uM.ABA$WT == T)[-1],] ## remove duplicated WT

## Clean up environment
rm(packages, install_if_missing, doubles, singles, all_variants, synonymous, wildtype)


## 2. Normalise growth rates between blocks, using stops & synonymous variants ##
#################################################################################

## 0 µM ABA

### Error-weighted mean of growth rate estimates
block1.0uM.ABA$gr_over_sigmasquared <- block1.0uM.ABA$growthrate/(block1.0uM.ABA$growthrate_sigma)**2
block1.0uM.ABA$one_over_sigmasquared <- 1/(block1.0uM.ABA$growthrate_sigma)**2
block2.0uM.ABA$gr_over_sigmasquared <- block2.0uM.ABA$growthrate/(block2.0uM.ABA$growthrate_sigma)**2
block2.0uM.ABA$one_over_sigmasquared <- 1/(block2.0uM.ABA$growthrate_sigma)**2
block3.0uM.ABA$gr_over_sigmasquared <- block3.0uM.ABA$growthrate/(block3.0uM.ABA$growthrate_sigma)**2
block3.0uM.ABA$one_over_sigmasquared <- 1/(block3.0uM.ABA$growthrate_sigma)**2
block4.0uM.ABA$gr_over_sigmasquared <- block4.0uM.ABA$growthrate/(block4.0uM.ABA$growthrate_sigma)**2
block4.0uM.ABA$one_over_sigmasquared <- 1/(block4.0uM.ABA$growthrate_sigma)**2

### Error-weighted mean of stop mutations
stops.block1.0uM.ABA <- sum(block1.0uM.ABA[which(block1.0uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block1.0uM.ABA[which(block1.0uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block2.0uM.ABA <- sum(block2.0uM.ABA[which(block2.0uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block2.0uM.ABA[which(block2.0uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block3.0uM.ABA <- sum(block3.0uM.ABA[which(block3.0uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block3.0uM.ABA[which(block3.0uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block4.0uM.ABA <- sum(block4.0uM.ABA[which(block4.0uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block4.0uM.ABA[which(block4.0uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)

### Error-weighted mean of WT/synonymous mutations
syn.block1.0uM.ABA <- sum(block1.0uM.ABA[which(block1.0uM.ABA$Nham_aa == 0 & (block1.0uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block1.0uM.ABA[which(block1.0uM.ABA$Nham_aa == 0 & (block1.0uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block2.0uM.ABA <- sum(block2.0uM.ABA[which(block2.0uM.ABA$Nham_aa == 0 & (block2.0uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block2.0uM.ABA[which(block2.0uM.ABA$Nham_aa == 0 & (block2.0uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block3.0uM.ABA <- sum(block3.0uM.ABA[which(block3.0uM.ABA$Nham_aa == 0 & (block3.0uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block3.0uM.ABA[which(block3.0uM.ABA$Nham_aa == 0 & (block3.0uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block4.0uM.ABA <- sum(block4.0uM.ABA[which(block4.0uM.ABA$Nham_aa == 0 & (block4.0uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block4.0uM.ABA[which(block4.0uM.ABA$Nham_aa == 0 & (block4.0uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)

### Calculate coefficients for linear transformation
scaling.0uM.ABA <- cbind(data.frame(rbind(stops.block1.0uM.ABA, syn.block1.0uM.ABA)),
                         data.frame(rbind(stops.block2.0uM.ABA, syn.block2.0uM.ABA)),
                         data.frame(rbind(stops.block3.0uM.ABA, syn.block3.0uM.ABA)),
                         data.frame(rbind(stops.block4.0uM.ABA, syn.block4.0uM.ABA)))
rownames(scaling.0uM.ABA) <- c("stops", "syn")
colnames(scaling.0uM.ABA) <- c("block1", "block2", "block3", "block4")

### Linear transformation
coefficients.0uM.ABA <- matrix(NA, ncol = 4, nrow = 2)
colnames(coefficients.0uM.ABA) <-  colnames(scaling.0uM.ABA)
for (i in 1:4){

  tmp.lm <- lm(formula = c(0, 100) ~ scaling.0uM.ABA[,i])
  tmp.lm <- summary(tmp.lm)
  coefficients.0uM.ABA[1,i] <- tmp.lm$coefficients[[2]]
  coefficients.0uM.ABA[2,i] <- tmp.lm$coefficients[[1]]  
  
}
block1.0uM.ABA$gr_normalised <- block1.0uM.ABA$growthrate*coefficients.0uM.ABA[1,1] + coefficients.0uM.ABA[2,1]
block2.0uM.ABA$gr_normalised <- block2.0uM.ABA$growthrate*coefficients.0uM.ABA[1,2] + coefficients.0uM.ABA[2,2]
block3.0uM.ABA$gr_normalised <- block3.0uM.ABA$growthrate*coefficients.0uM.ABA[1,3] + coefficients.0uM.ABA[2,3]
block4.0uM.ABA$gr_normalised <- block4.0uM.ABA$growthrate*coefficients.0uM.ABA[1,4] + coefficients.0uM.ABA[2,4]
block1.0uM.ABA$gr_sigma_normalised <- block1.0uM.ABA$growthrate_sigma*coefficients.0uM.ABA[1,1]
block2.0uM.ABA$gr_sigma_normalised <- block2.0uM.ABA$growthrate_sigma*coefficients.0uM.ABA[1,2]
block3.0uM.ABA$gr_sigma_normalised <- block3.0uM.ABA$growthrate_sigma*coefficients.0uM.ABA[1,3]
block4.0uM.ABA$gr_sigma_normalised <- block4.0uM.ABA$growthrate_sigma*coefficients.0uM.ABA[1,4]

## 250 µM ABA

### Error-weighted mean of growth rate estimates
block1.250uM.ABA$gr_over_sigmasquared <- block1.250uM.ABA$growthrate/(block1.250uM.ABA$growthrate_sigma)**2
block1.250uM.ABA$one_over_sigmasquared <- 1/(block1.250uM.ABA$growthrate_sigma)**2
block2.250uM.ABA$gr_over_sigmasquared <- block2.250uM.ABA$growthrate/(block2.250uM.ABA$growthrate_sigma)**2
block2.250uM.ABA$one_over_sigmasquared <- 1/(block2.250uM.ABA$growthrate_sigma)**2
block3.250uM.ABA$gr_over_sigmasquared <- block3.250uM.ABA$growthrate/(block3.250uM.ABA$growthrate_sigma)**2
block3.250uM.ABA$one_over_sigmasquared <- 1/(block3.250uM.ABA$growthrate_sigma)**2
block4.250uM.ABA$gr_over_sigmasquared <- block4.250uM.ABA$growthrate/(block4.250uM.ABA$growthrate_sigma)**2
block4.250uM.ABA$one_over_sigmasquared <- 1/(block4.250uM.ABA$growthrate_sigma)**2

### Error-weighted mean of stop mutations
stops.block1.250uM.ABA <- sum(block1.250uM.ABA[which(block1.250uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block1.250uM.ABA[which(block1.250uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block2.250uM.ABA <- sum(block2.250uM.ABA[which(block2.250uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block2.250uM.ABA[which(block2.250uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block3.250uM.ABA <- sum(block3.250uM.ABA[which(block3.250uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block3.250uM.ABA[which(block3.250uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
stops.block4.250uM.ABA <- sum(block4.250uM.ABA[which(block4.250uM.ABA[,"STOP"] == T),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block4.250uM.ABA[which(block4.250uM.ABA[,"STOP"] == T),,drop = F]$one_over_sigmasquared, na.rm = TRUE)

### Error-weighted mean of WT/synonymous mutations
syn.block1.250uM.ABA <- sum(block1.250uM.ABA[which(block1.250uM.ABA$Nham_aa == 0 & (block1.250uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block1.250uM.ABA[which(block1.250uM.ABA$Nham_aa == 0 & (block1.250uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block2.250uM.ABA <- sum(block2.250uM.ABA[which(block2.250uM.ABA$Nham_aa == 0 & (block2.250uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block2.250uM.ABA[which(block2.250uM.ABA$Nham_aa == 0 & (block2.250uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block3.250uM.ABA <- sum(block3.250uM.ABA[which(block3.250uM.ABA$Nham_aa == 0 & (block3.250uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block3.250uM.ABA[which(block3.250uM.ABA$Nham_aa == 0 & (block3.250uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)
syn.block4.250uM.ABA <- sum(block4.250uM.ABA[which(block4.250uM.ABA$Nham_aa == 0 & (block4.250uM.ABA$Nham_nt != 1)),,drop = F]$gr_over_sigmasquared, na.rm = TRUE) / sum(block4.250uM.ABA[which(block4.250uM.ABA$Nham_aa == 0 & (block4.250uM.ABA$Nham_nt != 1)),,drop = F]$one_over_sigmasquared, na.rm = TRUE)

### Calculate coefficients for linear transformation
scaling.250uM.ABA <- cbind(data.frame(rbind(stops.block1.250uM.ABA, syn.block1.250uM.ABA)),
                           data.frame(rbind(stops.block2.250uM.ABA, syn.block2.250uM.ABA)),
                           data.frame(rbind(stops.block3.250uM.ABA, syn.block3.250uM.ABA)),
                           data.frame(rbind(stops.block4.250uM.ABA, syn.block4.250uM.ABA)))
rownames(scaling.250uM.ABA) <- c("stops", "syn")
colnames(scaling.250uM.ABA) <- c("block1", "block2", "block3", "block4")

### Scale all of the concentrations
coefficients.250uM.ABA <- matrix(NA, ncol = 4, nrow = 2)
colnames(coefficients.250uM.ABA) <-  colnames(scaling.250uM.ABA)
for (i in 1:4){
  
  tmp.lm <- lm(formula = c(0, 100) ~ scaling.250uM.ABA[,i])
  tmp.lm <- summary(tmp.lm)
  coefficients.250uM.ABA[1,i] <- tmp.lm$coefficients[[2]]
  coefficients.250uM.ABA[2,i] <- tmp.lm$coefficients[[1]]  
  
}
block1.250uM.ABA$gr_normalised <- block1.250uM.ABA$growthrate*coefficients.250uM.ABA[1,1] + coefficients.250uM.ABA[2,1]
block2.250uM.ABA$gr_normalised <- block2.250uM.ABA$growthrate*coefficients.250uM.ABA[1,2] + coefficients.250uM.ABA[2,2]
block3.250uM.ABA$gr_normalised <- block3.250uM.ABA$growthrate*coefficients.250uM.ABA[1,3] + coefficients.250uM.ABA[2,3]
block4.250uM.ABA$gr_normalised <- block4.250uM.ABA$growthrate*coefficients.250uM.ABA[1,4] + coefficients.250uM.ABA[2,4]
block1.250uM.ABA$gr_sigma_normalised <- block1.250uM.ABA$growthrate_sigma*coefficients.250uM.ABA[1,1]
block2.250uM.ABA$gr_sigma_normalised <- block2.250uM.ABA$growthrate_sigma*coefficients.250uM.ABA[1,2]
block3.250uM.ABA$gr_sigma_normalised <- block3.250uM.ABA$growthrate_sigma*coefficients.250uM.ABA[1,3]
block4.250uM.ABA$gr_sigma_normalised <- block4.250uM.ABA$growthrate_sigma*coefficients.250uM.ABA[1,4]

## Clean up environment
rm(i, tmp.lm)


## 3. Combine PYL1 blocks ##
############################

## 0 µM ABA

### Extend the WT aa sequence of each PYL1 block by the other blocks
block1.0uM.ABA$aa_seq <- paste0(block1.0uM.ABA$aa_seq,
                                "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block2.0uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                block2.0uM.ABA$aa_seq,
                                "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block3.0uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                block3.0uM.ABA$aa_seq,
                                "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block4.0uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                block4.0uM.ABA$aa_seq)

### Merge PYL1 blocks
PYL1.PYL1.0uM.ABA <- rbind(block1.0uM.ABA, block2.0uM.ABA, block3.0uM.ABA, block4.0uM.ABA)

### Annotate the individual mutation types
PYL1.PYL1.0uM.ABA <- cbind("Pos" = rep(NA, nrow(PYL1.PYL1.0uM.ABA)), 
                      "WT_AA" = rep(NA, nrow(PYL1.PYL1.0uM.ABA)), 
                      "Mut" = rep(NA, nrow(PYL1.PYL1.0uM.ABA)),PYL1.PYL1.0uM.ABA)
for(i in 1:nrow(PYL1.PYL1.0uM.ABA)){
  
  #### Do not act on WT seqs
  if(PYL1.PYL1.0uM.ABA[i,"Nham_aa"] == 0){
    
    next
    
  }
  
  #### Find the position and mutation
  WT.seq <- PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA[,"WT"] == T)[1],"aa_seq"]
  WT.seq <- as.character(str_split_fixed(WT.seq, "", 177))
  mut.seq <- PYL1.PYL1.0uM.ABA[i,"aa_seq"]
  mut.seq <- as.character(str_split_fixed(mut.seq, "", 177))
  PYL1.PYL1.0uM.ABA[i,"Pos"] <- which(mut.seq != WT.seq)
  PYL1.PYL1.0uM.ABA[i,"WT_AA"] <- WT.seq[PYL1.PYL1.0uM.ABA[i,"Pos"]]
  PYL1.PYL1.0uM.ABA[i,"Mut"] <- mut.seq[PYL1.PYL1.0uM.ABA[i,"Pos"]]
  
} 

### Positional adjustment with respect to the full-length PYL1 protein
PYL1.PYL1.0uM.ABA$Pos <- PYL1.PYL1.0uM.ABA$Pos + 32

### WT specification
PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0),"Pos"] <- "WT"
PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0),"WT_AA"] <- "WT"
PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0),"Mut"] <- "WT"

## 250 µM ABA

### Extend the WT aa sequence of each PYL1 block by the other blocks
block1.250uM.ABA$aa_seq <- paste0(block1.250uM.ABA$aa_seq,
                                  "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                  "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                  "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block2.250uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                  block2.250uM.ABA$aa_seq,
                                  "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                  "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block3.250uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                  "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                  block3.250uM.ABA$aa_seq,
                                  "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN")

block4.250uM.ABA$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                  "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                  "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                  block4.250uM.ABA$aa_seq)

### Merge PYL1 blocks
PYL1.PYL1.250uM.ABA <- rbind(block1.250uM.ABA, block2.250uM.ABA, block3.250uM.ABA, block4.250uM.ABA)

### Annotate the individual mutation types
PYL1.PYL1.250uM.ABA <- cbind("Pos" = rep(NA, nrow(PYL1.PYL1.250uM.ABA)), 
                             "WT_AA" = rep(NA, nrow(PYL1.PYL1.250uM.ABA)), 
                             "Mut" = rep(NA, nrow(PYL1.PYL1.250uM.ABA)),PYL1.PYL1.250uM.ABA)
for(i in 1:nrow(PYL1.PYL1.250uM.ABA)){
  
  #### do not act on WT seqs
  if(PYL1.PYL1.250uM.ABA[i,"Nham_aa"] == 0){
    
    next
    
  }
  
  #### find the position and mutation
  WT.seq <- PYL1.PYL1.250uM.ABA[which(PYL1.PYL1.250uM.ABA[,"WT"] == T)[1],"aa_seq"]
  WT.seq <- as.character(str_split_fixed(WT.seq, "", 177))
  mut.seq <- PYL1.PYL1.250uM.ABA[i,"aa_seq"]
  mut.seq <- as.character(str_split_fixed(mut.seq, "", 177))
  PYL1.PYL1.250uM.ABA[i,"Pos"] <- which(mut.seq != WT.seq)
  PYL1.PYL1.250uM.ABA[i,"WT_AA"] <- WT.seq[PYL1.PYL1.250uM.ABA[i,"Pos"]]
  PYL1.PYL1.250uM.ABA[i,"Mut"] <- mut.seq[PYL1.PYL1.250uM.ABA[i,"Pos"]]
  
} 

### Positional adjustment with respect to the full-length PYL1 protein
PYL1.PYL1.250uM.ABA$Pos <- PYL1.PYL1.250uM.ABA$Pos + 32

### WT specification
PYL1.PYL1.250uM.ABA[which(PYL1.PYL1.250uM.ABA[,"Nham_aa"] == 0),"Pos"] <- "WT"
PYL1.PYL1.250uM.ABA[which(PYL1.PYL1.250uM.ABA[,"Nham_aa"] == 0),"WT_AA"] <- "WT"
PYL1.PYL1.250uM.ABA[which(PYL1.PYL1.250uM.ABA[,"Nham_aa"] == 0),"Mut"] <- "WT"


## 4. Scale growth rates and growth rate errors by wildtype ##
##############################################################

## Scale abundance fitness to WT
PYL1.PYL1.0uM.ABA$gr_normalised_WTscaled <- PYL1.PYL1.0uM.ABA$gr_normalised
PYL1.PYL1.0uM.ABA$gr_sigma_normalised_WTscaled <- PYL1.PYL1.0uM.ABA$gr_sigma_normalised
PYL1.PYL1.250uM.ABA$gr_normalised_WTscaled <- PYL1.PYL1.250uM.ABA$gr_normalised
PYL1.PYL1.250uM.ABA$gr_sigma_normalised_WTscaled <- PYL1.PYL1.250uM.ABA$gr_sigma_normalised

## Save PYL1.PYL1 dataframes as .Rdata files
save(PYL1.PYL1.0uM.ABA, file = "../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
save(PYL1.PYL1.250uM.ABA, file = "../../data/DiMSum/PYL1-PYL1/250uM_ABA/PYL1-PYL1_250uM_ABA_preprocessed.RData")


## 5. Version ##
################

# sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Madrid
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] stringr_1.5.2
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1  magrittr_2.0.4  cli_3.6.5       tools_4.5.1     glue_1.8.0      vctrs_0.6.5     stringi_1.8.7   lifecycle_1.0.4 rlang_1.1.6   