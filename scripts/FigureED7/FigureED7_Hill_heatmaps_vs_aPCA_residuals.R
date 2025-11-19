# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#####################################################################################
## Extended Data Figure 7 - Hill parameter heatmaps before vs. after normalisation ##
#####################################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "pheatmap", "grid", "gridExtra")

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

## input Hill parameter distributions from dose-response curve fits, filtering like in Figure 2
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
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

### prepare matrices
# WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))
# hill.mat <- vector(mode = "list", length = 4)
# names(hill.mat) <- c("B0", "Binf", "EC50", "Hill")
# for(i in 1:4){
#   
#   hill.mat[[i]] <- matrix(NA, nrow = 177, ncol = 21)
#   rownames(hill.mat[[i]]) <- paste0(WT.seq, 33:209)
#   colnames(hill.mat[[i]]) <- c("G", "A", "V", "L", "M", "I", "F", 
#                                "Y", "W", "K", "R", "H", "D", "E", 
#                                "S", "T", "C", "N", "Q", "P","*")
#   
# }
# 
# ### fill-in
# for(i in 2:nrow(parameters.Hill)){
#   
#   tmp.WT <- substr(rownames(parameters.Hill)[i], 1, nchar(rownames(parameters.Hill)[i]) - 1)
#   tmp.MUT <- substr(rownames(parameters.Hill)[i], nchar(rownames(parameters.Hill)[i]), nchar(rownames(parameters.Hill)[i]))
#   hill.mat[["B0"]][tmp.WT,tmp.MUT] <- parameters.Hill[i,"B[0]"]
#   hill.mat[["Binf"]][tmp.WT,tmp.MUT] <- parameters.Hill[i,"B[inf]"]
#   hill.mat[["EC50"]][tmp.WT,tmp.MUT] <- log10(parameters.Hill[i,"EC50"])
#   hill.mat[["Hill"]][tmp.WT,tmp.MUT] <- log10(parameters.Hill[i,"Hill"])
#   
# }
# 
# ### add in WT
# for(i in 1:length(WT.seq)){
#   
#   hill.mat[["B0"]][i,WT.seq[i]] <- parameters.Hill["WT","B[0]"]
#   hill.mat[["Binf"]][i,WT.seq[i]] <- parameters.Hill["WT","B[inf]"]
#   hill.mat[["EC50"]][i,WT.seq[i]] <- log10(parameters.Hill["WT","EC50"])
#   hill.mat[["Hill"]][i,WT.seq[i]] <- log10(parameters.Hill["WT","Hill"])
#   
# }


## 3. Calculate aPCA residuals and prepare heatmaps ##
######################################################

## LOESS of all four parameters vs. pseudo-abundance
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

### prepare matrices
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))
loess.res.mat <- vector(mode = "list", length = 4)
names(loess.res.mat) <- colnames(loess.res)
for(i in 1:4){
  
  loess.res.mat[[i]] <- matrix(NA, nrow = 177, ncol = 21)
  rownames(loess.res.mat[[i]]) <- paste0(WT.seq, 33:209)
  colnames(loess.res.mat[[i]]) <- c("G", "A", "V", "L", "M", "I", "F", 
                                    "Y", "W", "K", "R", "H", "D", "E", 
                                    "S", "T", "C", "N", "Q", "P","*")
  
}

### fill-in
for(i in 2:nrow(loess.res)){
  
  tmp.WT <- substr(rownames(loess.res)[i], 1, nchar(rownames(loess.res)[i]) - 1)
  tmp.MUT <- substr(rownames(loess.res)[i], nchar(rownames(loess.res)[i]), nchar(rownames(loess.res)[i]))
  loess.res.mat[["B0"]][tmp.WT,tmp.MUT] <- loess.res[i,"B0"]
  loess.res.mat[["Binf"]][tmp.WT,tmp.MUT] <- loess.res[i,"Binf"]
  loess.res.mat[["EC50"]][tmp.WT,tmp.MUT] <- loess.res[i,"EC50"]
  loess.res.mat[["Hill"]][tmp.WT,tmp.MUT] <- loess.res[i,"Hill"]
  
}

### add in WT
for(i in 1:length(WT.seq)){

  loess.res.mat[["B0"]][i,WT.seq[i]] <- loess.res["WT","B0"]
  loess.res.mat[["Binf"]][i,WT.seq[i]] <- loess.res["WT","Binf"]
  loess.res.mat[["EC50"]][i,WT.seq[i]] <- loess.res["WT","EC50"]
  loess.res.mat[["Hill"]][i,WT.seq[i]] <- loess.res["WT","Hill"]
  
}

### WT dashes
WT.mat <- t(loess.res.mat[["B0"]])
for(i in 1:ncol(WT.mat)){
  
  WT.mat[,i] <- rep("",nrow(WT.mat))
  WT.mat[WT.seq[i],i] <- "-"
  
}
colnames(WT.mat) <- paste0(WT.seq, 33:209)

### rename position IDs by the actual positions in PYL1
class(loess.res.mat[["B0"]]) <- "numeric"
class(loess.res.mat[["Binf"]]) <- "numeric"
class(loess.res.mat[["EC50"]]) <- "numeric"
class(loess.res.mat[["Hill"]]) <- "numeric"


## 4. Preparing heatmap annotations ##
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

## summarise annotation layers
meta <- matrix(NA, ncol = 4, nrow = 177)
colnames(meta) <- c("ABA contact", "ABI1 contact      ", "Alpha/Beta", "Loop")
rownames(meta) <- rownames(loess.res.mat$B0)
meta[,"ABI1 contact      "] <- ABI1_cont
meta[,"ABA contact"] <- ABA_cont
meta[,"Loop"] <- loop
meta[,"Alpha/Beta"] <- alpha.beta
class(meta) <- "numeric"
meta <- as.data.frame(meta)


## 5. Plot ##
#############

p.out <- vector(mode = "list", length = 4)

p.out[[1]] <- pheatmap(t(loess.res.mat$B0),
                       color = colorRampPalette(c("red", "gray90", "blue"))(n = 1000),
                       breaks = seq(f = -30, to = 30, length.out = 1000),
                       show_rownames = T,
                       show_colnames = T,
                       cluster_rows = F,
                       cluster_cols = F,
                       legend = F,
                       legend_breaks = seq(-30, 30, by = 15), 
                       legend_labels = seq(-30, 30, by = 15),
                       na_col = "white",
                       display_numbers = WT.mat,
                       annotation_legend = F,
                       annotation_col = meta,
                       annotation_names_col = F,
                       annotation_colors = list(`ABI1 contact      ` = colorRampPalette(c("white", "black"))(n = 2),
                                                `ABA contact` = colorRampPalette(c("white", "black"))(n = 2),
                                                `Alpha/Beta` = colorRampPalette(c("white", "grey30", "grey60"))(n = 3),
                                                `Loop` = colorRampPalette(c("white", "black"))(n = 2)),
                       border_color = "black", fontsize_row = 18, fontsize_col = 16,
                       angle_col = 90, treeheight_row = 100, treeheight_col = 100,
                       fontsize = 20, height = 8, width = 40, cellwidth = 15,
                       main = "")

p.out[[2]] <- pheatmap(t(loess.res.mat$Binf),
                       color = colorRampPalette(c("red", "gray90", "blue"))(n = 1000),
                       breaks = seq(f = -30, to = 30, length.out = 1000),
                       show_rownames = T,
                       show_colnames = T,
                       cluster_rows = F,
                       cluster_cols = F,
                       legend = F,
                       legend_breaks = seq(-30, 30, by = 15), 
                       legend_labels = seq(-30, 30, by = 15),
                       na_col = "white",
                       display_numbers = WT.mat,
                       annotation_legend = F,
                       annotation_col = meta,
                       annotation_names_col = F,
                       annotation_colors = list(`ABI1 contact      ` = colorRampPalette(c("white", "black"))(n = 2),
                                                `ABA contact` = colorRampPalette(c("white", "black"))(n = 2),
                                                `Alpha/Beta` = colorRampPalette(c("white", "grey30", "grey60"))(n = 3),
                                                `Loop` = colorRampPalette(c("white", "black"))(n = 2)),
                       border_color = "black", fontsize_row = 18, fontsize_col = 16,
                       angle_col = 90, treeheight_row = 100, treeheight_col = 100,
                       fontsize = 20, height = 8, width = 40, cellwidth = 15,
                       main = "")

p.out[[3]] <- pheatmap(t(loess.res.mat$EC50),
                       color = colorRampPalette(c("red", "gray90", "blue"))(n = 1000),
                       breaks = seq(f = -5.5, to = 5.5, length.out = 1000),
                       show_rownames = T,
                       show_colnames = T,
                       cluster_rows = F,
                       cluster_cols = F,
                       legend = F,
                       legend_breaks = log10(c(10^-4, 10^-2, 10^0, 10^2, 10^4)), 
                       legend_labels = c(expression(10^-4), 
                                         expression(10^-2), 
                                         expression(10^-0),
                                         expression(10^2), 
                                         expression(10^4)),
                       na_col = "white",
                       display_numbers = WT.mat,
                       annotation_legend = F,
                       annotation_col = meta,
                       annotation_names_col = F,
                       annotation_colors = list(`ABI1 contact      ` = colorRampPalette(c("white", "black"))(n = 2),
                                                `ABA contact` = colorRampPalette(c("white", "black"))(n = 2),
                                                `Alpha/Beta` = colorRampPalette(c("white", "grey30", "grey60"))(n = 3),
                                                `Loop` = colorRampPalette(c("white", "black"))(n = 2)),
                       border_color = "black", fontsize_row = 18, fontsize_col = 16,
                       angle_col = 90, treeheight_row = 100, treeheight_col = 100,
                       fontsize = 20, height = 8, width = 40, cellwidth = 15,
                       main = "")

p.out[[4]] <- pheatmap(t(loess.res.mat$Hill),
                       color = colorRampPalette(c("red", "gray90", "blue"))(n = 1000),
                       breaks = seq(f = -0.7, to = 0.7, length.out = 1000),
                       show_rownames = T,
                       show_colnames = T,
                       cluster_rows = F,
                       cluster_cols = F,
                       legend = F,
                       legend_breaks = log10(c(10^-0.6, 10^-0.3, 10^0, 10^0.3, 10^0.6)), 
                       legend_labels = c(expression(10^-0.6), 
                                         expression(10^-0.3), 
                                         expression(10^0), 
                                         expression(10^0.3), 
                                         expression(10^0.6)),
                       na_col = "white",
                       display_numbers = WT.mat,
                       annotation_legend = F,
                       annotation_col = meta,
                       annotation_names_col = F,
                       annotation_colors = list(`ABI1 contact      ` = colorRampPalette(c("white", "black"))(n = 2),
                                                `ABA contact` = colorRampPalette(c("white", "black"))(n = 2),
                                                `Alpha/Beta` = colorRampPalette(c("white", "grey30", "grey60"))(n = 3),
                                                `Loop` = colorRampPalette(c("white", "black"))(n = 2)),
                       border_color = "black", fontsize_row = 18, fontsize_col = 16,
                       angle_col = 90, treeheight_row = 100, treeheight_col = 100,
                       fontsize = 20, height = 8, width = 40, cellwidth = 15,
                       main = "")

pdf("../../results/FigureED7/FigureED7_aPCA_residual_heatmaps.pdf", width = 40, height = 25)
grid.arrange(p.out[[1]]$gtable, 
             p.out[[2]]$gtable, 
             p.out[[3]]$gtable, 
             p.out[[4]]$gtable,
             ncol = 1, respect = T)
dev.off()


## 6. Version ##
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
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] gridExtra_2.3   pheatmap_1.0.13 scales_1.4.0    stringr_1.5.2  
# 
# loaded via a namespace (and not attached):
# [1] RColorBrewer_1.1-3 R6_2.6.1           farver_2.1.2       magrittr_2.0.4     gtable_0.3.6       glue_1.8.0         lifecycle_1.0.4    cli_3.6.5          vctrs_0.6.5        compiler_4.5.1    
# [11] tools_4.5.1        rlang_1.1.6        stringi_1.8.7   