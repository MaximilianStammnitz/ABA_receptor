# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################################
## Figure 3B - min. EC50 fold changes on PYL1-ABI1 structure ##
###############################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "bio3d")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Collect EC50 fold-changes within residues ##
##################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")

## filter: remove stops
Hill.parameters.filtered <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]
Hill.parameters.filtered <- Hill.parameters.filtered[which(Hill.parameters.filtered[,"R^2"] > 0.9),]

## residue-level: summarise EC50s and EC50 fold-changes
ec50_per_res <- matrix(NA, nrow = 177, ncol = 21)
rownames(ec50_per_res) <- paste0(strsplit("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "")[[1]],
                                 33:209)
colnames(ec50_per_res) <- c(c("G", "A", "V", "L", "M", "I", "F", 
                              "Y", "W", "K", "R", "H", "D", "E", 
                              "S", "T", "C", "N", "Q", "P","*"))
for(i in 2:nrow(Hill.parameters.filtered)){
  tmp.id <- rownames(Hill.parameters.filtered)[i]
  aa.id <- substr(tmp.id, 1, nchar(tmp.id) - 1)
  mut.id <- substr(tmp.id, nchar(tmp.id), nchar(tmp.id))
  ec50_per_res[aa.id,mut.id] <- Hill.parameters.filtered[i,"EC50"]
}
ec50_per_res_FC <- ec50_per_res / Hill.parameters.filtered["WT","EC50"]
ec50_per_res_FC_log10 <- log10(ec50_per_res_FC)
ec50_per_res_FC_log10.min <- apply(ec50_per_res_FC_log10[,-21], 1, min, na.rm = T) + 4.5


## 2. Project on complex structure ##
#####################################

PYL1_ABI1.pdb <- read.pdb(file = "../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/pdb3kdj.ent")
PYL1_ABI1_min_EC50_fc <- PYL1_ABI1.pdb

### overwrite B-factor values of PYL1
for (i in 1:length(ec50_per_res_FC_log10.min)){
  
  PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom[,"resno"] == as.numeric(substr(names(ec50_per_res_FC_log10.min), 2, nchar(names(ec50_per_res_FC_log10.min))))[i]),"b"] <- ec50_per_res_FC_log10.min[i]
  
}

PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom$chain == "A" & c(PYL1_ABI1_min_EC50_fc$atom[,"resno"] == 31 | PYL1_ABI1_min_EC50_fc$atom[,"resno"] == 32  | PYL1_ABI1_min_EC50_fc$atom[,"resno"] == 106)),"b"] <- 100

### overwrite the B-factor values of the corresponding ABI1 (chain B)
PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom$chain == "B"),"b"] <- 100

### export
write.pdb(PYL1_ABI1_min_EC50_fc, file = paste0("../../data/PDB/modified/Figure3B_min_EC50_FC.ent"))

### ChimeraX: 
# B-factor colouring from 2.5 (blue) to 4.5 (grey) to 6.5 (red)


## 3. Version ##
################

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