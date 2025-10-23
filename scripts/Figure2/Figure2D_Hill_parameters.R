# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#################################
## Figure 2D - Hill parameters ##
#################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "beeswarm")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input and filter dose response curve parameters ##
########################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## filter
Hill.parameters.filtered <- parameters.Hill[which(parameters.Hill[,"EC50 P"] < 0.1 & parameters.Hill[,"Hill P"] < 0.1 & parameters.Hill[,"R^2"] > 0.9),]
Hill.parameters.filtered <- Hill.parameters.filtered[-grep("[*]", rownames(Hill.parameters.filtered)),]


## 2. Plot parameter distributions ##
#####################################

cex.all <- 0.25

## B[0]
pdf("../../results/Figure2/Figure2D_B0.pdf", width = 5, height = 15)

par(mar = c(2, 12, 2, 3))

## boxplot frame
boxplot(Hill.parameters.filtered[,"B[0]"],
        ylim = c(-20, 120),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 1.3,
        las = 2,
        cex.main = 2)

## actual data
beeswarm(Hill.parameters.filtered[,"B[0]"], 
         method = 'compactswarm',
         pch = 16,
         cex = cex.all,
         add = T,
         col = "black")

## add median
segments(0.65, median(Hill.parameters.filtered[,"B[0]"]), 
         1.35, median(Hill.parameters.filtered[,"B[0]"]), lwd = 3)

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 100, length.out = 6),
     labels = seq(from = 0, to = 100, length.out = 6),
     las = 2,
     line = 1)
mtext(text = bquote(B[0] ~ "(Binding)"), side = 2, line = 7, cex = 5)

dev.off()


## B[inf]
pdf("../../results/Figure2/Figure2D_Binf.pdf", width = 5, height = 15)

par(mar = c(2, 12, 2, 3))

## boxplot frame
boxplot(Hill.parameters.filtered[,"B[inf]"],
        ylim = c(-20, 120),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 1.3,
        las = 2,
        cex.main = 2)

## actual data
beeswarm(Hill.parameters.filtered[,"B[inf]"], 
         method = 'compactswarm',
         pch = 16,
         cex = cex.all,
         add = T,
         col = "black")

## add median
segments(0.65, median(Hill.parameters.filtered[,"B[inf]"]), 
         1.35, median(Hill.parameters.filtered[,"B[inf]"]), lwd = 3)

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 100, length.out = 6),
     labels = seq(from = 0, to = 100, length.out = 6),
     las = 2,
     line = 1)
mtext(text = bquote(B[infinity] ~ "(Binding)"), side = 2, line = 7, cex = 5)

dev.off()


## EC50
pdf("../../results/Figure2/Figure2D_EC50.pdf", width = 5, height = 15)

par(mar = c(2, 12, 2, 3))

## boxplot frame
boxplot(Hill.parameters.filtered[,"EC50"],
        ylim = c(0.01,1000),
        log = "y",
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 1.3,
        las = 2,
        cex.main = 2)

## actual data
beeswarm(Hill.parameters.filtered[,"EC50"], 
         method = 'compactswarm',
         pch = 16,
         cex = cex.all,
         add = T,
         col = "black")

## add median
segments(0.65, median(Hill.parameters.filtered[,"EC50"]), 
         1.35, median(Hill.parameters.filtered[,"EC50"]), lwd = 3)

## Axis
axis(2, 
     cex.axis = 3,
     at = c(0.01, 0.1, 1, 10, 100, 1000), 
     labels = c(0.01, 0.1, 1, 10, 100, c("1,000")),
     las = 2,
     line = 1)

mtext(text = bquote(EC[50] ~ "((+)-ABA, µM)"), side = 2, line = 7, cex = 5)

dev.off()


## Hill
pdf("../../results/Figure2/Figure2D_Hill.pdf", width = 5, height = 15)

par(mar = c(2, 12, 2, 3))

## boxplot frame
boxplot(Hill.parameters.filtered[,"Hill"],
        ylim = c(0,2.82),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 1.3,
        las = 2,
        cex.main = 2)

## actual data
beeswarm(Hill.parameters.filtered[,"Hill"], 
         method = 'compactswarm',
         pch = 16,
         cex = cex.all,
         add = T,
         col = "black")

## add median
segments(0.65, median(Hill.parameters.filtered[,"Hill"]), 
         1.35, median(Hill.parameters.filtered[,"Hill"]), lwd = 3)

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 2.5, by = 0.5), 
     labels = c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5"),
     las = 2,
     line = 1)
mtext(text = bquote(n ~ "(Hill coefficient)"), side = 2, line = 7, cex = 5)

dev.off()


## 3. Version ##
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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] beeswarm_0.4.0 scales_1.4.0   stringr_1.5.1 
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.3     cli_3.6.5          tools_4.5.1       
# [6] RColorBrewer_1.1-3 glue_1.8.0         rstudioapi_0.17.1  farver_2.1.2       stringi_1.8.7     
# [11] lifecycle_1.0.4    rlang_1.1.6
