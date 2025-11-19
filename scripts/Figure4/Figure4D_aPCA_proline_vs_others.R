# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#########################################
## Figure 4D - aPCA proline vs. others ##
#########################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "beeswarm")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input aPCA data ##
########################

load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Mut"] == "*"),]


## 2. Stratify data ##
######################

## divide into proline and other
PYL1.list <- vector(mode = "list", length = 2)
names(PYL1.list) <- c("Proline", "Other")
for(i in 1:nrow(PYL1.PYL1.0uM.ABA)){
  
  if(PYL1.PYL1.0uM.ABA[i,"Pos"] == "WT"){
    
    next
    
  }else{
    
    tmp.mut <- PYL1.PYL1.0uM.ABA[i,"gr_normalised_WTscaled"]
    names(tmp.mut) <- paste0(PYL1.PYL1.0uM.ABA[i,"WT_AA"],
                             PYL1.PYL1.0uM.ABA[i,"Pos"],
                             PYL1.PYL1.0uM.ABA[i,"Mut"])
    tmp.type <- PYL1.PYL1.0uM.ABA[i,"Mut"]
    
    if(tmp.type == "P"){
      PYL1.list[[1]] <- c(PYL1.list[[1]], tmp.mut)
    }else{
      PYL1.list[[2]] <- c(PYL1.list[[2]], tmp.mut)
    }

  }

}


## 4. Plot ##
#############

wilcox.test(x = PYL1.list$Proline, y = PYL1.list$Other) ## P = 3.444e-12

pdf('../../results/Figure4/Figure4D_aPCA_proline_vs_others.pdf', height = 6, width = 12)
par(mar = c(6,13,1,0))

names(PYL1.list) <- c("Proline\n(N = 171)", "Other\n(N = 3,187)")
beeswarm(PYL1.list, 
         pch = 16, 
         method = "compactswarm",
         cex = 0.5,
         bty = "n",
         main = "", 
         cex.main = 3,
         col = c("black"),
         ylim = c(-20, 130),
         axes = F,
         cex.lab = 3,
         cex.axis = 2.8, 
         yaxt = "n")
axis(1, at = c(1, 2), labels = names(PYL1.list), 
     cex.axis = 2.5, tick = F, padj = 0.6)
axis(2, 
     at = seq(f = 0, t = 100, length.out = 6), 
     labels = seq(f = 0, t = 100, length.out = 6), las = 2,
     line = 1, cex.axis = 2)
mtext("Relative PYL1 Abundance", side = 2, line = 8, cex = 2.5)

dev.off()


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
# [1] beeswarm_0.4.0 scales_1.4.0   stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.4     cli_3.6.5          tools_4.5.1        RColorBrewer_1.1-3 glue_1.8.0         farver_2.1.2       stringi_1.8.7      lifecycle_1.0.4   
# [11] rlang_1.1.6