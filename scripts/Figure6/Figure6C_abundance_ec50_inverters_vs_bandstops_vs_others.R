# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################################
## Figure 6C - abundance/EC50 of inverters vs. others ##
########################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel",
              "wesanderson", "beeswarm", "readxl")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify hypersensitive and "inverter" variants from hierarchical clustering ##
#####################################################################################

## input hclust results from Table S3
cat.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S3_v2.xlsx"))
colnames(cat.out) <- cat.out[2,]
hypersensitive.clust <- cat.out[which(cat.out[,"Variant label"] == "Hypersensitive binding"
                                      | cat.out[,"Variant label"] == "Inverted Hill shape"
                                      | cat.out[,"Variant label"] == "Bandstop shape"),]

hyper.binders <- hypersensitive.clust[which(hypersensitive.clust[,"Variant label"] == "Hypersensitive binding"), 1]
inverted.binders <- hypersensitive.clust[which(hypersensitive.clust[,"Variant label"] == "Inverted Hill shape"), 1]
bandstop.binders <- hypersensitive.clust[which(hypersensitive.clust[,"Variant label"] == "Bandstop shape"), 1]


## 2. Fetch raw measurements ##
###############################

## abundance
abundance.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S2_v2.xlsx", sheet = 2))
rownames(abundance.out) <- abundance.out[,1]
abundance.out <- abundance.out[,-1]
colnames(abundance.out) <- c(0, 250)
abundance.out <- abundance.out[-c(1:3),]

## EC50
ec50.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S3_v2.xlsx", sheet = 1))
rownames(ec50.out) <- ec50.out[,1]
ec50.out <- ec50.out[,-1]
colnames(ec50.out) <- ec50.out[2,]
ec50.out <- ec50.out[-c(1:2),]


## 3. Plot aPCA beeswarms ##
############################

abundance.comp <- list("Hypersensitive" = as.numeric(abundance.out[match(hyper.binders, rownames(abundance.out)),"0"]),
                       "Inverters" = as.numeric(abundance.out[match(inverted.binders, rownames(abundance.out)),"0"]),
                       "Bandstops" = as.numeric(abundance.out[match(bandstop.binders, rownames(abundance.out)),"0"]),
                       "Others" = as.numeric(abundance.out[-sort(unique(c(which(rownames(abundance.out) %in% hyper.binders),
                                                                          which(rownames(abundance.out) %in% inverted.binders),
                                                                          which(rownames(abundance.out) %in% bandstop.binders)))),"0"]))

wilcox.test(x = abundance.comp$Hypersensitive, y = abundance.comp$Others) # p < 2.2e-16
wilcox.test(x = abundance.comp$Inverters, y = abundance.comp$Others) # p = 2.521e-07
wilcox.test(x = abundance.comp$Bandstops, y = abundance.comp$Others) # p = 2.546e-13

## set up colouring
beeswarm.cols <- c("G" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "P" = wes_palette("Darjeeling2", 6, type = "continuous")[6],
                   "W" = wes_palette("Darjeeling2", 6, type = "continuous")[5], "Y" = wes_palette("Darjeeling2", 6, type = "continuous")[5], 
                   "F" = wes_palette("Darjeeling2", 6, type = "continuous")[5],
                   "A" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "V" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "L" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "I" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "M" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "S" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "T" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "C" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "Q" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "N" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "K" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "R" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "H" = wes_palette("Darjeeling2", 6, type = "continuous")[2],
                   "D" = wes_palette("Darjeeling2", 6, type = "continuous")[1], "E" = wes_palette("Darjeeling2", 6, type = "continuous")[1])

labels.cols <- c(names(abundance.out[match(hyper.binders, rownames(abundance.out)),"0"]),
                 names(abundance.out[match(inverted.binders, rownames(abundance.out)),"0"]),
                 names(abundance.out[match(bandstop.binders, rownames(abundance.out)),"0"]),
                 names(abundance.out[-sort(unique(c(which(rownames(abundance.out) %in% hyper.binders),
                                                    which(rownames(abundance.out) %in% inverted.binders),
                                                    which(rownames(abundance.out) %in% bandstop.binders)))),"0"]))
labels.cols <- substr(labels.cols, nchar(labels.cols), nchar(labels.cols))
labels.cols <- beeswarm.cols[match(labels.cols, names(beeswarm.cols))]
labels.cols[which(is.na(labels.cols) == T)] <- "black"
names(labels.cols[which(is.na(names(labels.cols) == T))]) <- "*"

pdf("../../results/Figure6/Figure6C_aPCA_inverters_bandstops.pdf", height = 13, width = 28)
par(mar = c(8,13,1,0))

boxplot(abundance.comp, 
        cex = 1, 
        ylim = c(-20, 130), 
        outline = F, 
        horizontal = F,
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        border = "white",
        col = "white",
        frame = F,
        xlab = '',
        xaxt = 'n',
        yaxt = 'n')

## set up label positions (via beeswarm)
labels.pos1 <- beeswarm(abundance.comp$Hypersensitive, method = "swarm", do.plot = FALSE, cex = 2.5)
labels.pos2 <- beeswarm(abundance.comp$Inverters, method = "swarm", do.plot = FALSE, cex = 5)
labels.pos3 <- beeswarm(abundance.comp$Bandstops, method = "swarm", do.plot = FALSE, cex = 5)
labels.pos4 <- beeswarm(abundance.comp$Others, method = "swarm", do.plot = FALSE, cex = 0.7)

## add mutations to plot
points(x = labels.pos1$x, 
       y = labels.pos1$y,
       col = "black",
       pch = 16,
       cex = 2.5)

labels.out <- rownames(abundance.out[match(inverted.binders, rownames(abundance.out)),])
labels.out <- substr(labels.out, nchar(labels.out), nchar(labels.out))
points(x = labels.pos2$x + 1, 
       y = labels.pos2$y, 
       col = beeswarm.cols[match(labels.out, names(beeswarm.cols))],
       pch = 16,
       cex = 5)

labels.out <- rownames(abundance.out[match(bandstop.binders, rownames(abundance.out)),])
labels.out <- substr(labels.out, nchar(labels.out), nchar(labels.out))
points(x = labels.pos3$x + 2, 
       y = labels.pos3$y, 
       col = beeswarm.cols[match(labels.out, names(beeswarm.cols))],
       pch = 16,
       cex = 5)

points(x = labels.pos4$x + 3, 
       y = labels.pos4$y,
       col = "black", 
       pch = 16,
       cex = 0.7)

## add X-axis
names(abundance.comp) <- c("Hypersensitive\n(N = 92)", 
                           "Inverters\n(N = 15)", 
                           "Band-stops\n(N = 22)", 
                           "Others\n(N = 3,412)")

axis(1, at = 1, labels = names(abundance.comp)[1], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.6, tick = F)
axis(1, at = 2, labels = names(abundance.comp)[2], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.6, tick = F)
axis(1, at = 3, labels = names(abundance.comp)[3], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.6, tick = F)
axis(1, at = 4, labels = names(abundance.comp)[4], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.6, tick = F)

## add Y-axis
axis(2, at = seq(f = 0, t = 100, length.out = 6), 
     labels = seq(f = 0, t = 100, length.out = 6), las = 2,
     cex.axis = 3.5, lwd = 4)

## Y label
mtext(text = "Relative PYL1 Abundance", side = 2, line = 8, cex = 5.5)

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
# [1] readxl_1.4.5      beeswarm_0.4.0    wesanderson_0.3.7 ggrepel_0.9.6     ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4     S7_0.2.0           glue_1.8.0         plyr_1.8.9         gridtext_0.1.5     cellranger_1.1.0  
# [11] grid_4.5.1         tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2       R6_2.6.1          
# [21] tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6       xml2_1.4.0   