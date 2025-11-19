# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##############################################
## Figure 3E - EC50 distance of key mutants ##
##############################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "beeswarm")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Calculate log EC50 fold change ##
#######################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")

## remove stops
Hill.parameters.filtered <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.9),]
Hill.parameters.filtered <- Hill.parameters.filtered[-grep("[*]", rownames(Hill.parameters.filtered)),]

## calculate log10-fold change between WT and variant
ec50_logFC <- log10(Hill.parameters.filtered[-1,"EC50"] / Hill.parameters.filtered["WT","EC50"])

## fetch EC50s and values
all_ec50 <- log10(Hill.parameters.filtered[-1,"EC50"])
all_ec50_se <- Hill.parameters.filtered[-1,"EC50 SE"] / c(Hill.parameters.filtered[-1,"EC50"] * log(10))
ref_ec50 <- log10(Hill.parameters.filtered["WT","EC50"])
ref_se <- Hill.parameters.filtered["WT","EC50 SE"] / c(Hill.parameters.filtered["WT","EC50"] * log(10))

## calculate Z-scores
z_scores <- (all_ec50 - ref_ec50) / sqrt(all_ec50_se^2 + ref_se^2)

## calculate p-values (BH adjusted)
p.ec50 <- 2 * (1 - pnorm(abs(z_scores)))
p.adj.ec50 <- p.adjust(p.ec50, method = "BH")

## summarise data
mutants <- as.data.frame(cbind("logFC" = ec50_logFC,
                               "p" = p.ec50,
                               "p.adj" = p.adj.ec50))
mutants$p.adj[mutants$p.adj == 0] <- 1e-18

## annotate known hypersensitive mutants from the literature
## 18 positions
## Mosquna et al., PNAS 2011: 
### H87P/R/G/A/I/W/M/J/V, V110F/P/L, I111Q/H/P/E/K, L114F, A116W, L188T/V/C/I
### F189V/A, T192F, L196Y/F, K200W
## Jones et al., eLife 2014:
### H87P
## Elzinga et al., ACS ChemBio 2019: 
### F88L/M, V108I, I137C/S, E171L, A190C/V
## Rowe et al., Nature Plants 2023: 
### S112A, E141D, R143S, A190V
known.mutant.pos <- c("H87", "F88", "V108", "V110",
                      "I111", "S112", "L114", "A116",
                      "I137", "E141", "R143", "E171",
                      "L188", "F189", "A190", "T192",
                      "L196", "K200")

## annotate new hypersensitive mutants from the data
## 19 positions
new.mutant.pos <- subset(mutants, `p.adj` < 0.1 & `logFC` < 0)
new.mutant.pos <- unique(substr(rownames(new.mutant.pos), 1, nchar(rownames(new.mutant.pos)) - 1))
new.mutant.pos <- new.mutant.pos[!new.mutant.pos %in% known.mutant.pos]

## all other PYL1 positions (137)
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))
WT.seq <- paste0(WT.seq, 33:209)
other.pos <- WT.seq
other.pos <- other.pos[!other.pos %in% known.mutant.pos]
other.pos <- other.pos[!other.pos %in% new.mutant.pos]


## 2. Annotate residue-specific distances/contacts to (+)-ABA and ABI1 ##
#########################################################################

## input residue minimal distances from ABA and ABI1
load("../../data/PYL1_distances_to_interfaces.Rdata")

ABA_dist <- PYL1_ABA_ABI1_dist$ABA$ABA_vs_PYL1_HAmin_ligand[-c(1:2)]
ABI1_dist <- PYL1_ABA_ABI1_dist$ABI1$ABI1_vs_PYL1_HAmin_ligand[-c(1:2)]
names(ABA_dist) <- names(ABI1_dist) <- 33:209

## add to data summary
data.all <- cbind("ABA distance" = ABA_dist,
                  "ABI1 distance" = ABI1_dist)
rownames(data.all) <- WT.seq
data.all <- cbind(data.all,
                  "Min. distance" = apply(data.all, 1, min))
data.all <- cbind(data.all,
                  "Category" = rep(NA, nrow(data.all)))
data.all[which(rownames(data.all) %in% known.mutant.pos),"Category"] <- "Known"
data.all[which(rownames(data.all) %in% new.mutant.pos),"Category"] <- "Novel"
data.all[which(rownames(data.all) %in% other.pos),"Category"] <- "Other"

## convert to list object
beeswarm.out <- vector(mode = "list", length = 3)
names(beeswarm.out) <- c("Known", "Novel", "Other")
for(i in 1:length(beeswarm.out)){
  beeswarm.out[[i]] <- data.all[data.all[,"Category"] == names(beeswarm.out)[i],"Min. distance"]
  beeswarm.out[[i]] <- as.numeric(beeswarm.out[[i]])
  names(beeswarm.out[[i]]) <- names(data.all[data.all[,"Category"] == names(beeswarm.out)[i],"Min. distance"])
}
names(beeswarm.out) <- c(expression("Known\nEC"["50"]*" regulator\n(N = 18)"),
                         expression("Novel\nEC"["50"]*" regulator\n(N = 19)"),
                         expression("\nOther sites\n(N = 140)"))


## 3. Plot ##
#############

pdf("../../results/Figure3/Figure3E_Distance_pos_horizontal.pdf", width = 17, height = 14)
par(mar = c(14,16,1,0))

boxplot(rev(beeswarm.out), 
        cex = 1, 
        ylim = c(0, 32), 
        outline = F, 
        horizontal = T,
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

beeswarm(rev(beeswarm.out),
         horizontal = T,
         method = "compactswarm",
         col = "black", pch = 16,
         cex = 3, add = T)

## add X-axis
text(y = c(3.2, 2.2, 1.1),
     x = -4,
     labels = c("Known", "Novel", "Other sites"),
     xpd = NA, cex = 3.3)

text(y = c(3, 2, 0.9),
     x = -4,
     labels = c(expression("EC"[50]*" regulators"),
                expression("EC"[50]*" regulators"),
                expression("(N = 140)")),
     xpd = NA, cex = 3.3)

text(y = c(2.8, 1.8),
     x = -4,
     labels = c("(N = 18)",
                "(N = 19)"),
     xpd = NA, cex = 3.3)

## vertical line
abline(v = 5, col = "black", lwd = 3)

## add Y-axis
axis(1, at = seq(f = 0, t = 30, length.out = 7),
     labels = seq(f = 0, t = 30, length.out = 7),
     las = 1,
     cex.axis = 3,
     lwd = 5,
     padj = 0.85)

## Y label
mtext(text = "Minimum site distance [Å]\n((+)-ABA and ABI1 interfaces)", 
      side = 1,
      line = 12, 
      cex = 4.5)

dev.off()


## 4. Version ##
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
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.4     cli_3.6.5          tools_4.5.1        RColorBrewer_1.1-3 glue_1.8.0         farver_2.1.2       vctrs_0.6.5        stringi_1.8.7     
# [11] lifecycle_1.0.4    rlang_1.1.6  