# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 6B - aPCA residual volcano plots ##
###########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel", "cowplot")

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

## input Hill parameter distributions from dose-response curve fits, filtering like in Figure 2D/E
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")
parameters.Hill <- parameters.Hill[which(parameters.Hill[,"EC50 P"] < 0.05 & parameters.Hill[,"Hill P"] < 0.05 & parameters.Hill[,"R^2"] > 0.95),]
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


## 3. Calculate residuals & significance ##
###########################################

## LOESS of all four parameters vs. pseudo-abundance
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

## summarise by position
loess.res.pos <- vector(mode = "list", length = 4)
names(loess.res.pos) <- colnames(loess.res)
loess.res.pos <- lapply(loess.res.pos, 
                        function(x){x <- matrix(NA, nrow = 177, ncol = 20); 
                        rownames(x) <- unique(substr(rownames(loess.res), 1, nchar(rownames(loess.res))-1))[-1]; 
                        colnames(x) <- c("G", "A", "V", "L", "M", "I", "F", 
                                         "Y", "W", "K", "R", "H", "D", "E", 
                                         "S", "T", "C", "N", "Q", "P"); return(x)})
for(i in 2:nrow(loess.res)){
  print(i)
  tmp.pos <- substr(rownames(loess.res)[i], 1, nchar(rownames(loess.res)[i])-1)
  tmp.mut <- substr(rownames(loess.res)[i], nchar(rownames(loess.res)[i]), nchar(rownames(loess.res)[i]))
  loess.res.pos$B0[tmp.pos,tmp.mut] <- loess.res[i,"B0"]
  loess.res.pos$Binf[tmp.pos,tmp.mut] <- loess.res[i,"Binf"]
  loess.res.pos$EC50[tmp.pos,tmp.mut] <- loess.res[i,"EC50"]
  loess.res.pos$Hill[tmp.pos,tmp.mut] <- loess.res[i,"Hill"]
}

## position-wise significance test
loess.res.pos.sig <- matrix(NA, nrow = 177, ncol = 4)
rownames(loess.res.pos.sig) <- rownames(loess.res.pos$B0)
colnames(loess.res.pos.sig) <- names(loess.res.pos)
for(i in 1:nrow(loess.res.pos$B0)){
  loess.res.pos.sig[i,"B0"] <- wilcox.test(x = loess.res.pos$B0[-i,], 
                                           y = loess.res.pos$B0[i,], 
                                           alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"Binf"] <- wilcox.test(x = loess.res.pos$Binf[-i,], 
                                             y = loess.res.pos$Binf[i,], 
                                             alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"EC50"] <- wilcox.test(x = loess.res.pos$EC50[-i,], 
                                             y = loess.res.pos$EC50[i,], 
                                             alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"Hill"] <- wilcox.test(x = loess.res.pos$Hill[-i,], 
                                             y = loess.res.pos$Hill[i,], 
                                             alternative = "two.sided")$p.value
}
loess.res.pos.sig.fdr <- loess.res.pos.sig
loess.res.pos.sig.fdr[,"B0"] <- p.adjust(loess.res.pos.sig[,"B0"], method = "BH")
loess.res.pos.sig.fdr[,"Binf"] <- p.adjust(loess.res.pos.sig[,"Binf"], method = "BH")
loess.res.pos.sig.fdr[,"EC50"] <- p.adjust(loess.res.pos.sig[,"EC50"], method = "BH")
loess.res.pos.sig.fdr[,"Hill"] <- p.adjust(loess.res.pos.sig[,"Hill"], method = "BH")

## take mean position-wise residual
loess.res.pos <- as.data.frame(cbind("B0" = apply(loess.res.pos$B0, 1, mean, na.rm = T),
                                     "B0 FDR" = loess.res.pos.sig.fdr[,"B0"],
                                     "Binf" = apply(loess.res.pos$Binf, 1, mean, na.rm = T),
                                     "Binf FDR" = loess.res.pos.sig.fdr[,"Binf"],
                                     "EC50" = exp(apply(loess.res.pos$EC50, 1, mean, na.rm = T)),
                                     "EC50 FDR" = loess.res.pos.sig.fdr[,"EC50"],
                                     "Hill" = exp(apply(loess.res.pos$Hill, 1, mean, na.rm = T)),
                                     "Hill FDR" = loess.res.pos.sig.fdr[,"Hill"]))


## 4. Plots ##
##############

axis.tick.size <- 30
axis.X.label.size <- 27
axis.Y.label.size <- 40
residue.label.size <- 9
point.size <- 7

## B0
p0 <- ggplot(loess.res.pos, aes(x = B0, y = -log10(`B0 FDR`))) +
  geom_point(aes(x = B0, y = -log10(`B0 FDR`), 
                 fill = squish(B0, range = c(-30, 30))),
             color = "black",
             shape = 21,
             size = point.size,
             alpha = 1) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-30,30),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text_repel(data = loess.res.pos[grep("T33|E36|F37|Q39|S43", rownames(loess.res.pos)),],
                  aes(label = rownames(loess.res.pos[grep("T33|E36|F37|Q39|S43", rownames(loess.res.pos)),])),
                  size = residue.label.size, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 2,
                  force_pull = 0.5,
                  box.padding = 0.85,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = seq(f = -15, t = 15, by = 5), 
                     limits = c(-15, 15), oob = scales::squish) +
  scale_y_continuous(breaks = seq(f = 0, t = 10, by = 2), 
                     limits = c(0, 10), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = axis.tick.size),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = axis.X.label.size, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = axis.Y.label.size, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = bquote("Position-wise mean " * B[0] * " % residual"),
       y = bquote(-log[10] * "(adjusted P)"))

## Binf
pinf <- ggplot(loess.res.pos, aes(x = Binf, y = -log10(`Binf FDR`))) +
  geom_point(aes(x = Binf, y = -log10(`Binf FDR`), 
                 fill = squish(Binf, range = c(-30, 30))),
             color = "black",
             shape = 21,
             size = point.size,
             alpha = 1) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-30,30),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text_repel(data = loess.res.pos[grep("D128|D129|R130|R131|E160|E161|E162|E163|R164", rownames(loess.res.pos)),],
                  aes(label = rownames(loess.res.pos[grep("D128|D129|R130|R131|E160|E161|E162|E163|R164", rownames(loess.res.pos)),])),
                  size = residue.label.size, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 5,
                  force_pull = 5,
                  box.padding = 2,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = seq(f = -60, t = 20, by = 20), 
                     limits = c(-65, 20), oob = scales::squish) +
  scale_y_continuous(breaks = seq(f = 0, t = 10, by = 2), 
                     limits = c(0, 10), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = axis.tick.size),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = axis.X.label.size, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = axis.Y.label.size, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = bquote("Position-wise mean " * B[infinity] * " % residual"),
       y = bquote(-log[10] * "(adjusted P)"))

## EC50
pEC50 <- ggplot(loess.res.pos, aes(x = EC50, y = -log10(`EC50 FDR`))) +
  geom_point(aes(x = log10(`EC50`), y = -log10(`EC50 FDR`), 
                 fill = squish(log10(`EC50`), range = c(-5.5, 5.5))),
             color = "black",
             shape = 21,
             size = point.size,
             alpha = 1) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-5.5,5.5),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text_repel(data = loess.res.pos[grep("H87|L196|Y50|E171|D107", rownames(loess.res.pos)),],
                  aes(x = log10(`EC50`), y = -log10(`EC50 FDR`),
                      label = rownames(loess.res.pos[grep("H87|L196|Y50|E171|D107", rownames(loess.res.pos)),])),
                  size = residue.label.size, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 2,
                  force_pull = 0.5,
                  box.padding = 0.85,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = seq(f = -5, t = 5, by = 1), 
                     limits = c(-2, 2),
                     labels = c("A", "B", "C", "0.01", "0.1", "1", "10", "100", "I", "J", "K")) +
  scale_y_continuous(breaks = seq(f = 0, t = 10, by = 2), 
                     limits = c(0, 10), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = axis.tick.size),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = axis.X.label.size, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = axis.Y.label.size, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = bquote("Position-wise mean " * EC[50] * " fold change residual"),
       y = bquote(-log[10] * "(adjusted P)"))

## Hill
pHill <- ggplot(loess.res.pos, aes(x = Hill, y = -log10(`Hill FDR`))) +
  geom_point(aes(x = log10(`Hill`), y = -log10(`Hill FDR`), 
                 fill = squish(log10(`Hill`), range = c(-0.7, 0.7))),
             color = "black",
             shape = 21,
             size = point.size,
             alpha = 1) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-0.7,0.7),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text_repel(data = loess.res.pos[grep("Q42|H154|S41|P68|S43", rownames(loess.res.pos)),],
                  aes(x = log10(`Hill`), y = -log10(`Hill FDR`),
                      label = rownames(loess.res.pos[grep("Q42|H154|S41|P68|S43", rownames(loess.res.pos)),])),
                  size = residue.label.size, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 2,
                  force_pull = 0.5,
                  box.padding = 0.85,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = c(log10(0.5), log10(1), log10(2)),
                     limits = c(log10(0.5), log10(2)),
                     labels = c("0.5", "1", "2")) +
  scale_y_continuous(breaks = seq(f = 0, t = 10, by = 2), 
                     limits = c(0, 10), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = axis.tick.size),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = axis.X.label.size, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = axis.Y.label.size, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = "Position-wise mean Hill coefficient (n) fold change residual",
       y = bquote(-log[10] * "(adjusted P)"))

## plot all
pdf('../../results/FigureS6/FigureS6B_Hill_paramater_residuals_volcanos.pdf', width = 12, height = 42)

plot_grid(p0, pinf, pEC50, pHill, 
          align = "hv", axis = "tblr", ncol = 1)

dev.off()


## 5. Version ##
################

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
# [1] cowplot_1.1.3 ggrepel_0.9.6 ggtext_0.1.2  ggplot2_3.5.1 reshape_0.8.9 scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.4     generics_0.1.3    labeling_0.4.3   
# [8] glue_1.8.0        colorspace_2.1-1  plyr_1.8.9        gridtext_0.1.5    grid_4.4.1        munsell_0.5.1     tibble_3.2.1     
# [15] lifecycle_1.0.4   compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2     
# [22] R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    withr_3.0.2       tools_4.4.1       gtable_0.3.6     
# [29] xml2_1.3.6       