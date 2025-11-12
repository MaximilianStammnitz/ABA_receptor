# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 12.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#####################################################################
## Figure 5C - aPCA residuals to Hill parameters on PYL1 structure ##
#####################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "bio3d")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input aPCA data and bPCA Hill parameters ##
#################################################

load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Mut"] == "*"),]

## input Hill parameter distributions from dose-response curve fits, filtering like in Figure 2D/E
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")
parameters.Hill <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.9),]
parameters.Hill <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]


## 2. Combine data ##
#####################

parameters.Hill <- cbind(parameters.Hill, "aPCA" = rep(NA, nrow(parameters.Hill)))
for (i in 1:nrow(parameters.Hill)){
  
  print(i)
  tmp.res <- rownames(parameters.Hill)[i]
  if(tmp.res == "WT"){
    
    parameters.Hill[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT == T)[1],"gr_normalised_WTscaled"]
    
  }else{
    
    tmp.pos <- as.numeric(substr(tmp.res, 2, nchar(tmp.res)-1))
    tmp.WT <- substr(tmp.res, 1, 1)
    tmp.MUT <- substr(tmp.res, nchar(tmp.res), nchar(tmp.res))
    parameters.Hill[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[which(as.numeric(PYL1.PYL1.0uM.ABA$Pos) == tmp.pos & PYL1.PYL1.0uM.ABA$Mut == tmp.MUT),"gr_normalised_WTscaled"]
    
  }
  
}


## 3. Calculate residuals & significance ##
###########################################

## LOESS of all four parameters vs. pseudo-abundance
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

## summarise by position
loess.res.pos <- vector(mode = "list", length = 4)
names(loess.res.pos) <- colnames(loess.res)
loess.res.pos <- lapply(loess.res.pos, 
                        function(x){x <- matrix(NA, nrow = 176, ncol = 20); 
                        rownames(x) <- unique(substr(rownames(loess.res), 1, nchar(rownames(loess.res))-1))[-1]; 
                        colnames(x) <- c("G", "A", "V", "L", "M", "I", "F", 
                                         "Y", "W", "K", "R", "H", "D", "E", 
                                         "S", "T", "C", "N", "Q", "P"); return(x)})
for(i in 2:nrow(loess.res)){
  print(i)
  tmp.pos <- substr(rownames(loess.res)[i], 1, nchar(rownames(loess.res)[i])-1)
  tmp.mut <- substr(rownames(loess.res)[i], nchar(rownames(loess.res)[i]), nchar(rownames(loess.res)[i]))
  loess.res.pos$B0[tmp.pos,tmp.mut] <- loess.res[i,"B0"]
  loess.res.pos$Binf[tmp.pos,tmp.mut] <- loess.res[i,"Binf"]
  loess.res.pos$EC50[tmp.pos,tmp.mut] <- loess.res[i,"EC50"]
  loess.res.pos$Hill[tmp.pos,tmp.mut] <- loess.res[i,"Hill"]
}

## take the mean
loess.res.pos.mean <- cbind(apply(loess.res.pos$B0, 1, mean, na.rm = T),
                            apply(loess.res.pos$Binf, 1, mean, na.rm = T),
                            apply(loess.res.pos$EC50, 1, mean, na.rm = T),
                            apply(loess.res.pos$Hill, 1, mean, na.rm = T))
colnames(loess.res.pos.mean) <- c("B0", "Binf", "log(EC50)", "log(Hill)")


## 4. Put this on the 3D structure ##
#####################################

PYL1.pdb <- read.pdb(file = "../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/pdb3kdj.ent")
PYL1.B0.pdb <- PYL1.Binf.pdb <- PYL1.EC50.pdb <- PYL1.Hill.pdb <- PYL1.pdb
for (i in 1:nrow(loess.res.pos.mean)){
  
  PYL1.B0.pdb$atom[which(PYL1.B0.pdb$atom[,"resno"] == as.numeric(substr(rownames(loess.res.pos.mean), 2, nchar(rownames(loess.res.pos.mean))))[i]),"b"] <- loess.res.pos.mean[i,"B0"] + 30
  PYL1.Binf.pdb$atom[which(PYL1.Binf.pdb$atom[,"resno"] == as.numeric(substr(rownames(loess.res.pos.mean), 2, nchar(rownames(loess.res.pos.mean))))[i]),"b"] <- loess.res.pos.mean[i,"Binf"] + 30
  PYL1.EC50.pdb$atom[which(PYL1.EC50.pdb$atom[,"resno"] == as.numeric(substr(rownames(loess.res.pos.mean), 2, nchar(rownames(loess.res.pos.mean))))[i]),"b"] <- loess.res.pos.mean[i,"log(EC50)"] + 5.5
  PYL1.Hill.pdb$atom[which(PYL1.Hill.pdb$atom[,"resno"] == as.numeric(substr(rownames(loess.res.pos.mean), 2, nchar(rownames(loess.res.pos.mean))))[i]),"b"] <- loess.res.pos.mean[i,"log(Hill)"] + 0.7
  
}

### overwrite the B-factor values of the corresponding ABI1 (chain B)
PYL1.B0.pdb$atom[which(PYL1.B0.pdb$atom$chain == "A" & c(PYL1.B0.pdb$atom[,"resno"] == 31 | PYL1.B0.pdb$atom[,"resno"] == 32 | PYL1.B0.pdb$atom[,"resno"] == 106)),"b"] <- 0
PYL1.B0.pdb$atom[which(PYL1.B0.pdb$atom$chain == "B"),"b"] <- 0
PYL1.Binf.pdb$atom[which(PYL1.Binf.pdb$atom$chain == "A" & c(PYL1.Binf.pdb$atom[,"resno"] == 31 | PYL1.Binf.pdb$atom[,"resno"] == 32 | PYL1.Binf.pdb$atom[,"resno"] == 106)),"b"] <- 0
PYL1.Binf.pdb$atom[which(PYL1.Binf.pdb$atom$chain == "B"),"b"] <- 0
PYL1.EC50.pdb$atom[which(PYL1.EC50.pdb$atom$chain == "A" & c(PYL1.EC50.pdb$atom[,"resno"] == 31 | PYL1.EC50.pdb$atom[,"resno"] == 32 | PYL1.EC50.pdb$atom[,"resno"] == 106)),"b"] <- 0
PYL1.EC50.pdb$atom[which(PYL1.EC50.pdb$atom$chain == "B"),"b"] <- 0
PYL1.Hill.pdb$atom[which(PYL1.Hill.pdb$atom$chain == "A" & c(PYL1.Hill.pdb$atom[,"resno"] == 31 | PYL1.Hill.pdb$atom[,"resno"] == 32 | PYL1.Hill.pdb$atom[,"resno"] == 106)),"b"] <- 0
PYL1.Hill.pdb$atom[which(PYL1.Hill.pdb$atom$chain == "B"),"b"] <- 0

### output
write.pdb(PYL1.B0.pdb, file = "../../data/PDB/modified/Figure5C_plots/Figure5C_B0.pdb")
write.pdb(PYL1.Binf.pdb, file = "../../data/PDB/modified/Figure5C_plots/Figure5C_Binf.pdb")
write.pdb(PYL1.EC50.pdb, file = "../../data/PDB/modified/Figure5C_plots/Figure5C_EC50.pdb")
write.pdb(PYL1.Hill.pdb, file = "../../data/PDB/modified/Figure5C_plots/Figure5C_Hill.pdb")


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
# [1] bio3d_2.4-5   scales_1.4.0  stringr_1.5.2
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.4     cli_3.6.5          parallel_4.5.1     tools_4.5.1        RColorBrewer_1.1-3 glue_1.8.0         farver_2.1.2       Rcpp_1.1.0        
# [11] stringi_1.8.7      grid_4.5.1         lifecycle_1.0.4    rlang_1.1.6      
