# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

################################################################
## Figure 3H - dendrogram of EC50 fold-changes and chemotypes ##
################################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "wesanderson", "dendextend")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Collect EC50 fold-changes within residues ##
##################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")

## filter like in Figure 2
Hill.parameters.filtered <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.8),]
Hill.parameters.filtered <- Hill.parameters.filtered[-grep("[*]", rownames(Hill.parameters.filtered)),]

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


## 2. Review task: clustering by chemotype ##
#############################################

## match mutation labels and colours
cols <- c("G" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "P" = wes_palette("Darjeeling2", 6, type = "continuous")[6],
          "W" = wes_palette("Darjeeling2", 6, type = "continuous")[5], "Y" = wes_palette("Darjeeling2", 6, type = "continuous")[5], 
          "F" = wes_palette("Darjeeling2", 6, type = "continuous")[5],
          "A" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "V" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "L" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
          "I" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "M" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
          "S" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "T" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "C" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
          "Q" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "N" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
          "K" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "R" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "H" = wes_palette("Darjeeling2", 6, type = "continuous")[2],
          "D" = wes_palette("Darjeeling2", 6, type = "continuous")[1], "E" = wes_palette("Darjeeling2", 6, type = "continuous")[1])

## hierarchical clustering
hc <- hclust(dist(log10(t(ec50_per_res_FC[, 1:20]))), "complete")
d <- as.dendrogram(hc)

## rotating branches
target_order <- c("F", "Y", "W", 
                  "V", "I", "L", "M", "A", 
                  "S", "T", "C", "N", "Q",
                  "H", "K", "R",
                  "D", "E",
                  "G", "P")
d_rot <- rotate(d, order = target_order)

## colouring tip labels
labs <- labels(d_rot)
labs_cols <- cols[match(labs, names(cols))]

## plot
pdf("../../results/Figure3/Figure3H_EC50_hclust_dendrogram.pdf", width = 16, height = 15)
plot(d_rot,
     leaflab = "none", 
     xaxt = "n", 
     yaxt = "n",
     xlab = "", 
     ylab = "", 
     main = "", 
     sub = "",
     bty = "n",
     edgePar = list(lwd = 7))

## add coloured amino acid labels
labs_ord  <- labels(d_rot)
cols_ord  <- labs_cols[labs_ord]
text(x = 1:20,
     y = -0.7, 
     labels = labels(d_rot),
     xpd = NA,
     col = cols_ord,
     cex = 4.5)

dev.off()


## 3. Version ##
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
# [1] dendextend_1.19.1 wesanderson_0.3.7 scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       viridis_0.6.5      magrittr_2.0.4     gtable_0.3.6       glue_1.8.0         gridExtra_2.3      tibble_3.3.0      
# [11] pkgconfig_2.0.3    generics_0.1.4     dplyr_1.1.4        lifecycle_1.0.4    ggplot2_4.0.0      cli_3.6.5          S7_0.2.0           viridisLite_0.4.2  vctrs_0.6.5        grid_4.5.1        
# [21] compiler_4.5.1     tools_4.5.1        pillar_1.11.1      rlang_1.1.6        stringi_1.8.7  