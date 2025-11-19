# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#############################################################
## Supplementary Figure 5B - LOESS fit vs. Hill parameters ##
#############################################################


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


## 3. Calculate residuals ##
############################

## LOESS of all four parameters vs. pseudo-abundance
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric") ## log?
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

## predict full-range LOESS curves
loess.pred <- expand.grid("aPCA" = seq(f = range(parameters.Hill[,"aPCA"])[1], 
                                       to = range(parameters.Hill[,"aPCA"])[2], 
                                       length.out = 1000))
loess.pred$B0 <- as.numeric(predict(loess.B0, newdata = loess.pred, se = TRUE)$fit)
loess.pred$Binf <- as.numeric(predict(loess.Binf, newdata = loess.pred, se = TRUE)$fit)
loess.pred$EC50 <- exp(as.numeric(predict(loess.EC50, newdata = loess.pred, se = TRUE)$fit))
loess.pred$Hill <- exp(as.numeric(predict(loess.Hill, newdata = loess.pred, se = TRUE)$fit)) ## log?

### correlations
r.B0 <- cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"B[0]"], method = "spearman")
r.Binf <- cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"B[inf]"], method = "spearman")
r.EC50 <- cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"EC50"], method = "spearman")
r.n <- cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"Hill"], method = "spearman")

## summarise data
Hill.df <- as.data.frame(cbind(parameters.Hill[,c(2,3,4,1,17)],
                               "B[0] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"B0"],
                               "B[inf] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Binf"],
                               "EC50 res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"EC50"],
                               "Hill res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Hill"]))
Hill.df$names <- rownames(Hill.df)
Hill.df$res <- as.factor(c("WT", substr(Hill.df$names, 1, nchar(Hill.df$names)-1)[-1]))

## B0
p1 <- ggplot(Hill.df) +
  geom_hex(aes(x = `aPCA`, y = `B[0]`), bins = 50, inherit.aes = F) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = B0), linewidth = 2) +
  scale_fill_gradientn(colours = colorRampPalette(c("gray90", "black"))(n = 30), limits = c(0, 30)) +
  annotate("text",
           x = 0,
           y = 115,
           label = bquote(rho == .(format(r.B0, nsmall = 2, digits = 2))),
           hjust = 0, size = 15, color = "black") +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 115)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 30),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Abundance",
       y = bquote(B[0]))

# Binf
p2 <- ggplot(Hill.df) +
  geom_hex(aes(x = `aPCA`, y = `B[inf]`), bins = 50, inherit.aes = F) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = Binf), linewidth = 2) +
  scale_fill_gradientn(colours = colorRampPalette(c("gray90", "black"))(n = 30), limits = c(0, 30)) +
  annotate("text",
           x = 0,
           y = 115,
           label = bquote(rho == .(format(r.Binf, nsmall = 2, digits = 2))),
           hjust = 0, size = 15, color = "black") +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 115)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 30),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Abundance",
       y = bquote(B[infinity]))

# EC50
p3 <- ggplot(Hill.df) +
  geom_hex(aes(x = `aPCA`, y = `EC50`), bins = 50, inherit.aes = F) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = EC50), linewidth = 2) +
  scale_fill_gradientn(colours = colorRampPalette(c("gray90", "black"))(n = 30), limits = c(0, 30)) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000), 
                limits = c(0.01, 20000),
                labels = c("0.01", "0.1", "1", "10", "100", "1,000", "10,000")) +
  annotate("text", x = 0, y = 10000, label = bquote(rho == .(format(r.EC50, nsmall = 2, digits = 2))), 
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 30),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Abundance",
       y = bquote(EC[50] ~ "(µM (+)-ABA)"))

# Hill
p4 <- ggplot(Hill.df) +
  geom_hex(aes(x = `aPCA`, y = `Hill`), bins = 50, inherit.aes = F) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = Hill), linewidth = 2) +
  scale_fill_gradientn(colours = colorRampPalette(c("gray90", "black"))(n = 30), limits = c(0, 30)) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3), 
                limits = c(0.3, 3.5)) +
  annotate("text", x = 0, y = 3, label = bquote(rho == .(format(r.n, nsmall = 2, digits = 2))), 
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 30),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Abundance",
       y = "Hill coefficient (n)")

pdf('../../results/FigureS5/FigureS5B_abundancePCA_vs_Hill_parameters.pdf', width = 45, height = 13)

plot_grid(p1, p2, p3, p4, align = "hv", axis = "tblr", ncol = 4)

dev.off()


## 4. Version ##
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
# [1] vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.4     generics_0.1.3    labeling_0.4.3    glue_1.8.0       
# [8] colorspace_2.1-1  plyr_1.8.9        gridtext_0.1.5    grid_4.4.1        munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4  
# [15] compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 lattice_0.22-6    farver_2.1.2     
# [22] R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    withr_3.0.2       tools_4.4.1       gtable_0.3.6     
# [29] hexbin_1.28.5     xml2_1.3.6    