# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

####################################################
## Figure 4B - PCA on full dose response data set ##
####################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "ggplot2", "ggtext")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input and filter dose response curve parameters ##
########################################################

## Input 12-dose measurements
load("../../data/DiMSum/PYL1-ABI1_12conc/PYL1-ABI1_preprocessed_12conc.RData")

## remove synonymous variants
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(x[,"Nham_aa"] == 0 & is.na(x[,"WT"]) == T),]; return(y)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(PYL1.ABI1[[1]][,"Nham_aa"] == 0 & is.na(PYL1.ABI1[[1]][,"WT"]) == F)[-1],]; return(y)})

## summarise in one matrix
PYL1.summary <- matrix(NA, nrow = 177*21, ncol = 12)
colnames(PYL1.summary) <- names(PYL1.ABI1)
ids <- c()
for(i in 1:177){
  tmp.pos <- rep(str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"], "", 177)[i], 21)
  ids <- c(ids, paste0(tmp.pos, i + 32))
}
rownames(PYL1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
                                            "I", "F", "Y", "W", "K",
                                            "R", "H", "D", "E", "S",
                                            "T", "C", "N", "Q", "P",
                                            "*"),177))
WT <- str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"],"",177)[1,]
for(j in 1:length(PYL1.ABI1)){
  
  tmp.vars <- str_split_fixed(PYL1.ABI1[[j]]$aa_seq,"",177)
  
  for(i in 1:nrow(tmp.vars)){
    
    ### only take the "best" WT
    if(all(tmp.vars[i,] == WT)){
      
      PYL1.summary[paste0(WT, 33:209, WT),j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
      next
      
    }
    
    pos.mut <- which(tmp.vars[i,] != WT)
    
    ### skip if higher-order mutant
    if(length(pos.mut) > 1){
      
      next
      
    }
    
    ### otherwise obtain the corresponding fitness
    tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
    PYL1.summary[tmp.mut,j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
    
  }
  
}

## remove WT repetition 
WTs.pos <- match(paste0(WT, 33:209, WT), rownames(PYL1.summary))
PYL1.summary.nonWT <- PYL1.summary[-WTs.pos[-1],]
PYL1.summary.nonWT <- PYL1.summary.nonWT[c(16,1:15,17:nrow(PYL1.summary.nonWT)),]
rownames(PYL1.summary.nonWT)[1] <- "WT"

### remove stops and any mutants with incomplete curves
PYL1.summary.nonWT <- PYL1.summary.nonWT[-sort(unique(do.call(c,apply(PYL1.summary.nonWT[,1:12], 2, function(x){which(is.na(x) == T)})))),1:12]
PYL1.summary.nonWT <- PYL1.summary.nonWT[-grep("[*]", rownames(PYL1.summary.nonWT)),]
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)


## 2. PCA ##
############

PCA <- prcomp(PYL1.summary.nonWT[,1:12], scale. = T, center = T)
sum.PCA <- summary(PCA)
y.PCA.components.vars <- sum.PCA$importance[2,]*100
y.PCA.components.vars <- round(y.PCA.components.vars, 1)
PCA.out <- as.data.frame(PCA$x)

pdf('../../results/Figure4/Figure4B_PCA_scree.pdf', height = 11.41/2, width = 10.26059)
par(mar = c(3,5,2,2))

plot(y.PCA.components.vars, xlim = c(1,6), ylim = c(0, 100), 
     pch = 16, cex = 2,
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n")
lines(y.PCA.components.vars, lwd = 2)
axis(1, at = 1:6, labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"), tick = T)
axis(2, at = seq(f = 0, to = 100, length.out = 6), labels = seq(f = 0, to = 100, length.out = 6), las = 2)
mtext(text = "Percentage of variance explained", side = 2, line = 3, cex = 2)

dev.off()

pdf('../../results/Figure4/Figure4B_PCA_PC1_loadings.pdf', height = 11.41/3, width = 14.616/2)
par(mar = c(6,6.5,2,0))

out <- barplot(rev(PCA$rotation[,1]),
               col = "black", border = "white",
               xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(x = c(0.7, 1.9, 3.1, 4.3, 
           5.5, 6.7, 7.9, 9.1, 
           10.3, 11.5, 12.7, 13.9), 
     y = par("usr")[3] - 0.04,
     labels = c("0.000", "0.009", "0.032", "0.111",
                "0.387", "1.360", "4.760", "16.66",
                "58.31", "204.08", "714.29", "2,500"), 
     srt = 45, xpd = TRUE)
axis(2, at = c(0, -0.05, -0.1, -0.15, -0.2, -0.25, -0.3), 
     labels = c("0.00", "-0.05", "-0.10", "-0.15", "-0.20", "-0.25", "-0.30"), tick = T, las = 2)
mtext(text = "(+)-ABA concentration (µM)", side = 1, line = 4, cex = 2)
mtext(text = "PC1 loading", side = 2, line = 4.5, cex = 2)

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
# [1] ggtext_0.1.2  ggplot2_4.0.0 scales_1.4.0  stringr_1.5.2
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        GGally_2.4.0       stringi_1.8.7      purrr_1.1.0        generics_0.1.4     S7_0.2.0           glue_1.8.0         gridtext_0.1.5    
# [11] grid_4.5.1         ggstats_0.11.0     tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    tidyr_1.3.1       
# [21] farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6       xml2_1.4.0