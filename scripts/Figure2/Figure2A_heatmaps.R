# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 22.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

################################################
## Figure 2A - binding heatmaps for all conc. ##
################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "pheatmap", "grid", "gridExtra", "cowplot")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processing DiMSum data ##
###################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## Build matrices
PYL1.ABI1.B <- vector(mode = "list", length = 12)
names(PYL1.ABI1.B) <- names(PYL1.ABI1)
PYL1.ABI1.B <- lapply(PYL1.ABI1.B, function(x){x <- matrix(NA, nrow = 177, ncol = 21); 
rownames(x) <- 1:177; colnames(x) <- c("G", "A", "V", "L", "M", "I", "F", 
                                       "Y", "W", "K", "R", "H", "D", "E", 
                                       "S", "T", "C", "N", "Q", "P","*"); 
return(x)})

## fill-in fitness matrix
for(j in 1:length(PYL1.ABI1)){
  
  for(i in 1:nrow(PYL1.ABI1[[j]])){
    
    if(PYL1.ABI1[[j]][i,"Pos"] == "WT"){
      
      next
      
    }else{
      
      PYL1.ABI1.B[[j]][as.character(c(as.numeric(PYL1.ABI1[[j]][i,"Pos"]) - 32)),PYL1.ABI1[[j]][i,"Mut"]] <- as.numeric(PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"])
      
    }
    
  }
  
}

## rename position IDs by the actual positions in PYL1
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))
PYL1.ABI1.B <- lapply(PYL1.ABI1.B, function(x){rownames(x) <- WT.seq; class(x) <- "numeric"; return(x)})

## calculate number of total measurements
sum(sapply(PYL1.ABI1.B, function(x){sum(!is.na(x))}))+12

## calculate range of mutations across concentrations
range(sapply(PYL1.ABI1.B, function(x){sum(!is.na(x))}))

## fill in WTs
for(j in 1:length(PYL1.ABI1.B)){
  
  for(i in 1:nrow(PYL1.ABI1.B[[j]])){
    
    PYL1.ABI1.B[[j]][i,WT.seq[i]] <- as.numeric(PYL1.ABI1[[j]][which(PYL1.ABI1[[j]]$WT == T)[1],"gr_normalised_WTscaled"])
    
  }
  
}

## wildtype "dot" annotations
WT.mat <- t(PYL1.ABI1.B[[1]])
for(i in 1:ncol(WT.mat)){
  WT.mat[,i] <- rep("",nrow(WT.mat))
  WT.mat[WT.seq[i],i] <- "-"
}


## 2. Preparing heatmap annotations ##
######################################

## Annotation layers of PYL1:
## 1. conservation, based on alignment of all Arabidopsis PYR/PYL proteins
## 2. core vs. surface (based on relative SASA)
## 3. alpha-helices & 4. beta-sheets
## 5. turns and loops - highlight gate & latch
## 6. (+)-ABA contact (Miyazono et al., Nature 2009)
## 7. ABI1 contact (Miyazono et al., Nature 2009)

## add secondary structure annotations
alpha.beta <- alpha <- beta <- loop <- ABA_cont <- ABI1_cont <- rep(0, 177)
ABI1_cont[c(c(87,88,90,111:117,142:144,178,180,181,185,186,188,189,192,193,196,199)-32)] <- 1
ABA_cont[c(c(86,88,110,116,121,135,137,143,144,147,149,171,189,193,197)-32)] <- 1
loop[c(c(111:118,140:147)-32)] <- 1
alpha.beta[c(c(36:47,69:76,82:84,184:208)-32)] <- 1
alpha.beta[c(c(55:64,91:94,105:110,119:127,132:137,148:157,165:176)-32)] <- 2

## add annotation layer: relative SASA
SASA <- read.table("../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/SASA_PyMol.txt", header = F)[,1:2]
SASA[,2] <- as.numeric(gsub("%", "", SASA[,2]))/100
SASA <- SASA[-c(1:2),]
SASA <- SASA[grep("/A/", SASA[,1]),]

## summarise annotation layers
meta <- matrix(NA, ncol = 5, nrow = 177)
colnames(meta) <- c("ABA contact", "ABI1 contact", "SASA", "Alpha.Beta", "Loop")
rownames(meta) <- 33:209
meta[,"ABI1 contact"] <- ABI1_cont
meta[,"ABA contact"] <- ABA_cont
meta[,"Loop"] <- loop
meta[,"Alpha.Beta"] <- alpha.beta
meta[,"SASA"] <- SASA[,2]
class(meta) <- "numeric"
meta <- as.data.frame(meta)


## 3. Plot ##
#############

p.out <- vector(mode = "list", length = 12)
for (i in 1:12){
  
  if(i != 12){

    p.out[[i]] <- pheatmap(t(PYL1.ABI1.B[[i]]),
                           color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
                           breaks = seq(f = 0, to = 120, length.out = 1000),
                           show_rownames = T,
                           show_colnames = F,
                           cluster_rows = F,
                           cluster_cols = F,
                           legend = T,
                           legend_breaks = seq(0, 100, by = 20), 
                           legend_labels = seq(0, 100, by = 20),
                           na_col = "grey",
                           display_numbers = WT.mat,
                           border_color = NA, fontsize_row = 11, fontsize_col = 17,
                           angle_col = 0, treeheight_row = 100, treeheight_col = 100,
                           fontsize = 20, height = 3.75, width = 35, cellwidth = 13.5,
                           main = "")

  }else if(i == 12){
    
    p.out[[i]] <- pheatmap(t(PYL1.ABI1.B[[i]]),
                           color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
                           breaks = seq(f = 0, to = 120, length.out = 1000),
                           show_rownames = T,
                           show_colnames = T,
                           cluster_rows = F,
                           cluster_cols = F,
                           legend = T,
                           legend_breaks = seq(0, 100, by = 20), 
                           legend_labels = seq(0, 100, by = 20),
                           na_col = "grey",
                           display_numbers = WT.mat,
                           annotation_legend = T,
                           border_color = NA, fontsize_row = 11, fontsize_col = 17,
                           angle_col = 0, treeheight_row = 100, treeheight_col = 100,
                           fontsize = 20, height = 3.75, width = 35, cellwidth = 13.5,
                           main = "")
    
  }
  
}

pdf("../../results/Figure2/Figure2A_combined_heatmaps.pdf", width = 35, height = 38.5)
grid.arrange(p.out[[1]]$gtable, p.out[[2]]$gtable, p.out[[3]]$gtable, 
             p.out[[4]]$gtable, p.out[[5]]$gtable, p.out[[6]]$gtable, 
             p.out[[7]]$gtable, p.out[[8]]$gtable, p.out[[9]]$gtable, 
             p.out[[10]]$gtable, p.out[[11]]$gtable, p.out[[12]]$gtable, 
             ncol = 1, respect = F)
dev.off()


## 4. Version ##
################

# sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
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
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] cowplot_1.2.0   gridExtra_2.3   pheatmap_1.0.13 scales_1.4.0    stringr_1.5.1  
# 
# loaded via a namespace (and not attached):
# [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.3     gtable_0.3.6      
# [7] glue_1.8.0         tibble_3.3.0       pkgconfig_2.0.3    generics_0.1.4     dplyr_1.1.4        ggplot2_4.0.0     
# [13] lifecycle_1.0.4    cli_3.6.5          S7_0.2.0           vctrs_0.6.5        compiler_4.5.1     rstudioapi_0.17.1 
# [19] tools_4.5.1        pillar_1.11.0      rlang_1.1.6        stringi_1.8.7    
