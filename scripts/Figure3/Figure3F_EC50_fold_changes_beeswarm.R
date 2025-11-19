# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################################
## Figure 3F - beeswarm of EC50 fold-changes and chemotypes  ##
###############################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "wesanderson", "beeswarm")

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


## 2. Plot beeswarm of mutations, colour by chemotype ##
########################################################

beeswarm.out <- vector(mode = "list", length = 4)
names(beeswarm.out) <- c("H87", "R164", "Q199", "T118")
for(i in 1:length(beeswarm.out)){
  beeswarm.out[[i]] <- ec50_per_res_FC[match(names(beeswarm.out)[i], rownames(ec50_per_res_FC)),-21]
  beeswarm.out[[i]] <- beeswarm.out[[i]][-which(is.na(beeswarm.out[[i]]))]
}

## extract the mutation labels
beeswarm.out.labels <- do.call(c,lapply(beeswarm.out, names))

## match mutation labels and colours
beeswarm.cols <- c("G" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "P" = wes_palette("Darjeeling2", 6, type = "continuous")[6],
                   "W" = wes_palette("Darjeeling2", 6, type = "continuous")[5], "Y" = wes_palette("Darjeeling2", 6, type = "continuous")[5], 
                   "F" = wes_palette("Darjeeling2", 6, type = "continuous")[5],
                   "A" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "V" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "L" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "I" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "M" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "S" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "T" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "C" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "Q" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "N" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "K" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "R" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "H" = wes_palette("Darjeeling2", 6, type = "continuous")[2],
                   "D" = wes_palette("Darjeeling2", 6, type = "continuous")[1], "E" = wes_palette("Darjeeling2", 6, type = "continuous")[1])

cols.out <- beeswarm.out.labels
cols.out <- beeswarm.cols[match(cols.out, names(beeswarm.cols))]

## plot
pdf("../../results/Figure3/Figure3F_EC50_FC_chemotype.pdf", width = 19, height = 14)
par(mar = c(6,14,1,0))

boxplot(beeswarm.out, 
        cex = 1, 
        log = "y", 
        ylim = c(0.01, 1000), 
        outline = F, 
        horizontal = F,
        medcol = 'darkgrey',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        border = "white",
        col = "white",
        lwd = 3,
        frame = F,
        xlab = '',
        xaxt = 'n',
        yaxt = 'n')

## dotted horizontal line: no fold change in EC50
abline(h = 1, col = "black", lty = "dotted", lwd = 2)

## set up label positions (via beeswarm)
labels.pos <- beeswarm(beeswarm.out, method = "swarm", 
                       do.plot = FALSE, cex = 6)

## add mutation labels to plot
text(x = labels.pos$x, y = labels.pos$y, 
     labels = beeswarm.out.labels, 
     col = cols.out, cex = 3.5)

## add X-axis
names(beeswarm.out) <- c("H87", "R164", "Q199", "T118")
axis(1, at = 1, labels = names(beeswarm.out)[1], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 2, labels = names(beeswarm.out)[2], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 3, labels = names(beeswarm.out)[3], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 4, labels = names(beeswarm.out)[4], las = 1, col = 'black', 
     cex.axis = 5, lwd = 3, col.ticks = F, padj = 0.5, tick = F)

## add Y-axis
axis(2, at = c(0.01, 0.1, 1, 10, 100, 1000),
     labels = c("0.01", "0.1", "1", "10", "100", "1,000"),
     las = 2,
     cex.axis = 3, lwd = 3)

## Y label
mtext(text = bquote(EC[50] * " fold change"), side = 2, line = 7, cex = 5)

## legend
legend("topleft", legend = c("Aromatic",
                             "Aliphatic",
                             "Polar, uncharged",
                             "+ charged",
                             "- charged",
                             "Special"),
       pch = 21,
       col = rep("black", 6),
       pt.bg = rev(wes_palette("Darjeeling2", 6, type = "continuous"))[c(2:6,1)], 
       bty = "n",
       cex = 4)

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
# [1] beeswarm_0.4.0    wesanderson_0.3.7 scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.4     cli_3.6.5          tools_4.5.1        RColorBrewer_1.1-3 glue_1.8.0         farver_2.1.2       stringi_1.8.7      lifecycle_1.0.4   
# [11] rlang_1.1.6  