# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

######################################
## Figure 4C - aPCA fitness heatmap ##
######################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "viridis", "pheatmap", "bio3d")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input pre-processed aPCA data ##
######################################

## Input aPCA
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
#load("../../data/DiMSum/PYL1-PYL1/250uM_ABA/PYL1-PYL1_250uM_ABA_preprocessed.RData")

## remove synonymous variants
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]
# PYL1.PYL1.250uM.ABA <- PYL1.PYL1.250uM.ABA[-which(PYL1.PYL1.250uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.250uM.ABA[,"WT"]) == T),]
# PYL1.PYL1.250uM.ABA <- PYL1.PYL1.250uM.ABA[-which(PYL1.PYL1.250uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.250uM.ABA[,"WT"]) == F)[-1],]

### for Supplementary Table 2B: summarise in one matrix
# PYL1.PYL1.summary <- matrix(NA, nrow = 177*21, ncol = 2)
# colnames(PYL1.PYL1.summary) <- c("0", "250")
# ids <- c()
# for(i in 1:177){
#   tmp.pos <- rep(str_split_fixed(PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT)[1],"aa_seq"], "", 177)[i], 21)
#   ids <- c(ids, paste0(tmp.pos, i + 32))
# }
# rownames(PYL1.PYL1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
#                                                  "I", "F", "Y", "W", "K",
#                                                  "R", "H", "D", "E", "S",
#                                                  "T", "C", "N", "Q", "P",
#                                                  "*"),177))
# WT <- str_split_fixed(PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT)[1],"aa_seq"], "", 177)[1,]

### fill in 0 µM (+)-ABA conc. results
# tmp.vars <- str_split_fixed(PYL1.PYL1.0uM.ABA$aa_seq,"",177)
# for(i in 1:nrow(tmp.vars)){
# 
#   ### only take the "best" WT
#   if(all(tmp.vars[i,] == WT)){
# 
#     PYL1.PYL1.summary[paste0(WT, 33:209, WT),"0"] <- PYL1.PYL1.0uM.ABA[i,"gr_normalised_WTscaled"]
#     next
# 
#   }
# 
#   pos.mut <- which(tmp.vars[i,] != WT)
# 
#   ### skip if higher-order mutant
#   if(length(pos.mut) > 1){
# 
#     next
# 
#   }
# 
#   ### otherwise obtain the corresponding fitness
#   tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
#   PYL1.PYL1.summary[tmp.mut,"0"] <- PYL1.PYL1.0uM.ABA[i,"gr_normalised_WTscaled"]
# 
# }

### fill in 250 µM (+)-ABA conc. results
# tmp.vars <- str_split_fixed(PYL1.PYL1.250uM.ABA$aa_seq,"",177)
# for(i in 1:nrow(tmp.vars)){
# 
#   ### only take the "best" WT
#   if(all(tmp.vars[i,] == WT)){
# 
#     PYL1.PYL1.summary[paste0(WT, 33:209, WT),"250"] <- PYL1.PYL1.250uM.ABA[i,"gr_normalised_WTscaled"]
#     next
# 
#   }
# 
#   pos.mut <- which(tmp.vars[i,] != WT)
# 
#   ### skip if higher-order mutant
#   if(length(pos.mut) > 1){
# 
#     next
# 
#   }
# 
#   ### otherwise obtain the corresponding fitness
#   tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
#   PYL1.PYL1.summary[tmp.mut,"250"] <- PYL1.PYL1.250uM.ABA[i,"gr_normalised_WTscaled"]
# 
# }

# ### remove WT repetition 
# WTs.pos <- match(paste0(WT, 33:209, WT), rownames(PYL1.PYL1.summary))
# PYL1.PYL1.summary.nonWT <- PYL1.PYL1.summary[-WTs.pos[-1],]
# PYL1.PYL1.summary.nonWT <- PYL1.PYL1.summary.nonWT[c(16,1:15,17:nrow(PYL1.PYL1.summary.nonWT)),]
# rownames(PYL1.PYL1.summary.nonWT)[1] <- "WT"
# write.table(PYL1.PYL1.summary.nonWT, file = "../../data/Supplementary Tables/TableS2B.txt", sep = "\t", quote = F, col.names = T, row.names = T)


## 2. Matrix format ##
######################

## WT
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))

### prepare matrix
PYL1.PYL1.mat <- matrix(NA, nrow = 177, ncol = 21)
rownames(PYL1.PYL1.mat) <- paste0(WT.seq, 33:209)
colnames(PYL1.PYL1.mat) <- c("G", "A", "V", "L", "M", "I", "F", 
                             "Y", "W", "K", "R", "H", "D", "E", 
                             "S", "T", "C", "N", "Q", "P","*")

### fill-in aPCA matrix
for(i in 1:nrow(PYL1.PYL1.0uM.ABA)){
  
  if(PYL1.PYL1.0uM.ABA[i,"WT_AA"] == "WT"){
    next
  }
  
  tmp.WT <- paste0(PYL1.PYL1.0uM.ABA[i,"WT_AA"], PYL1.PYL1.0uM.ABA[i,"Pos"])
  tmp.MUT <- PYL1.PYL1.0uM.ABA[i,"Mut"]
  PYL1.PYL1.mat[tmp.WT,tmp.MUT] <- PYL1.PYL1.0uM.ABA[i,"gr_normalised_WTscaled"]
  
}

### add in WT
for(i in 1:length(WT.seq)){
  PYL1.PYL1.mat[i,WT.seq[i]] <- as.numeric(PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT == T),"gr_normalised_WTscaled"])
}

### WT dashes
WT.mat <- t(PYL1.PYL1.mat)
for(i in 1:ncol(WT.mat)){
  WT.mat[,i] <- rep("",nrow(WT.mat))
  WT.mat[WT.seq[i],i] <- "-"
}
colnames(WT.mat) <- paste0(WT.seq, 33:209)

### rename position IDs by the actual positions in PYL1
class(PYL1.PYL1.mat) <- "numeric"


## 3. Heatmap ##
################

### add secondary structure annotations
alpha.beta <- rep(0, 177)
alpha.beta[c(c(36:47,69:76,82:84,184:208)-32)] <- 1
alpha.beta[c(c(55:64,91:94,105:110,119:127,132:137,148:157,165:176)-32)] <- 2

### add annotation layer: relative SASA
SASA <- read.table("../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/SASA_PyMol.txt", header = F)[,1:2]
SASA[,2] <- as.numeric(gsub("%", "", SASA[,2]))/100
SASA <- SASA[-c(1:2),]
SASA <- SASA[grep("/A/", SASA[,1]),]

### annotation layer: PYL1-PYL1 interface (get_contacts - based on Melcher or AF2?! double-check)
PYL1.PYL1.interface <- read.table("../../data/Predictions/GetContacts/PYL1_PYL1_GetContacts.tsv", header = F)[,2:4]
PYL1.PYL1.interface <- str_split_fixed(PYL1.PYL1.interface$V3, ":", 4)[,2:3]
PYL1.PYL1.interface <- unique(paste0(PYL1.PYL1.interface[,1], ":", PYL1.PYL1.interface[,2]))
PYL1.PYL1.interface <- PYL1.PYL1.interface[order(as.numeric(str_split_fixed(PYL1.PYL1.interface, ":", 2)[,2]))]

### summarise annotation layers
meta <- matrix(0, ncol = 3, nrow = 177)
colnames(meta) <- c("PYL1-PYL1", "SASA", "Alpha.Beta   ")
rownames(meta) <- rownames(PYL1.PYL1.mat)
meta[as.numeric(str_split_fixed(PYL1.PYL1.interface, ":", 2)[,2])-1,"PYL1-PYL1"] <- 1
meta[,"Alpha.Beta   "] <- alpha.beta
meta[,"SASA"] <- SASA[,2]
class(meta) <- "numeric"
meta <- as.data.frame(meta)

pheatmap(t(PYL1.PYL1.mat),
         color = magma(n = 1000, begin = 0, end = 1),
         breaks = seq(f = 0, to = 130, length.out = 1000),
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         legend = T,
         legend_breaks = seq(0, 100, by = 20), 
         legend_labels = seq(0, 100, by = 20),
         labels_col = substr(rownames(PYL1.PYL1.mat), 1, 1),
         na_col = "grey",
         display_numbers = WT.mat,
         annotation_col = meta,
         annotation_legend = F,
         annotation_colors = list(`SASA` = colorRampPalette(c("black", "grey95"))(n = 20),
                                  `Alpha.Beta   ` = colorRampPalette(c("white", "grey30", "grey60"))(n = 3),
                                  `PYL1-PYL1` = colorRampPalette(c("white", "black"))(n = 2)),
         border_color = NA, fontsize_row = 35, fontsize_col = 17,
         angle_col = 0, treeheight_row = 100, treeheight_col = 100,
         fontsize = 20, height = 18, width = 45, cellwidth = 17,
         main ="", filename = "../../results/Figure4/Figure4C_aPCA_heatmap.pdf")


## 4. Mean abundance per residue on structure ##
################################################

# ## Re-build matrix
# PYL1.PYL1.mat <- matrix(NA, nrow = 177, ncol = 21)
# rownames(PYL1.PYL1.mat) <- paste0(WT.seq, 33:209)
# colnames(PYL1.PYL1.mat) <- c("G", "A", "V", "L", "M", "I", "F",
#                              "Y", "W", "K", "R", "H", "D", "E",
#                              "S", "T", "C", "N", "Q", "P","*")
# 
# ### fill-in aPCA matrix
# for(i in 1:nrow(PYL1.PYL1.0uM.ABA)){
# 
#   if(PYL1.PYL1.0uM.ABA[i,"WT_AA"] == "WT"){
#     next
#   }
# 
#   tmp.WT <- paste0(PYL1.PYL1.0uM.ABA[i,"WT_AA"], PYL1.PYL1.0uM.ABA[i,"Pos"])
#   tmp.MUT <- PYL1.PYL1.0uM.ABA[i,"Mut"]
#   PYL1.PYL1.mat[tmp.WT,tmp.MUT] <- PYL1.PYL1.0uM.ABA[i,"gr_normalised_WTscaled"]
# 
# }
# 
# ## calculate mean abundance per position
# PYL1.PYL1.mean <- apply(PYL1.PYL1.mat[,-21], 1, mean, na.rm = T)
# 
# ## import PDB
# PYL1_PYL1.pdb <- read.pdb(file = "../../data/PDB/original/Melcher_2009_ApoPYL1/pdb3kay.ent")
# 
# ## overwrite PDB file B-factor values
# PYL1_PYL1.pdb$atom$b <- 0
# for (i in 4:length(PYL1.PYL1.mean)){
#   tmp.res <- substr(names(PYL1.PYL1.mean)[i], 2, nchar(names(PYL1.PYL1.mean)[i]))
#   PYL1_PYL1.pdb$atom[which(PYL1_PYL1.pdb$atom[,"resno"] == tmp.res & PYL1_PYL1.pdb$atom$chain == "A"),"b"] <- PYL1.PYL1.mean[i]
# }
# 
# ### export
# write.pdb(PYL1_PYL1.pdb, file = "../../data/PDB/modified/Figure4C_PYL1_PYL1.ent")


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
# [1] bio3d_2.4-5       pheatmap_1.0.13   viridis_0.6.5     viridisLite_0.4.2 stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4     S7_0.2.0           glue_1.8.0         gridExtra_2.3      scales_1.4.0       grid_4.5.1        
# [11] tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2       R6_2.6.1           tidyselect_1.2.1  
# [21] pillar_1.11.1      parallel_4.5.1     magrittr_2.0.4     tools_4.5.1        gtable_0.3.6       ggplot2_4.0.0