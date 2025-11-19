# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################
## Figure 5A - aPCA residual outliers ##
########################################


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

## input Hill parameter distributions from dose-response curve fits, filtering like in Figure 2
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
parameters.Hill <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.9),]
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
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

## assess raw correlations
cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"B[0]"], method = "spearman")
cor(x = parameters.Hill[,"aPCA"], y = parameters.Hill[,"B[inf]"], method = "spearman")
cor(x = parameters.Hill[,"aPCA"], y = log(parameters.Hill[,"EC50"]), method = "spearman")
cor(x = parameters.Hill[,"aPCA"], y = log(parameters.Hill[,"Hill"]), method = "spearman")

## predict full-range LOESS curves
loess.pred <- expand.grid("aPCA" = seq(f = range(parameters.Hill[,"aPCA"])[1], 
                                       to = range(parameters.Hill[,"aPCA"])[2], 
                                       length.out = 1000))
loess.pred$B0 <- as.numeric(predict(loess.B0, newdata = loess.pred, se = TRUE)$fit)
loess.pred$Binf <- as.numeric(predict(loess.Binf, newdata = loess.pred, se = TRUE)$fit)
loess.pred$EC50 <- exp(as.numeric(predict(loess.EC50, newdata = loess.pred, se = TRUE)$fit))
loess.pred$Hill <- exp(as.numeric(predict(loess.Hill, newdata = loess.pred, se = TRUE)$fit)) ## log?
rm(loess.B0, loess.Binf, loess.EC50, loess.Hill)

## summarise data
Hill.df <- as.data.frame(cbind(parameters.Hill[,c(2,3,4,1,17)],
                               "B[0] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"B0"],
                               "B[inf] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Binf"],
                               "EC50 res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"EC50"],
                               "Hill res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Hill"]))
Hill.df$res <- as.factor(c("WT", substr(Hill.df$names, 1, nchar(Hill.df$names)-1)[-1]))


## 4. Plot ##
#############

out.5A.B0 <- ggplot(Hill.df) +
  geom_point(data = Hill.df, 
             aes(x = `aPCA`, y = `B[0]`, color = `B[0] res`), size = 2, stroke = 0) +
  geom_point(data = Hill.df[grep("^E36", rownames(Hill.df)),],
             aes(x = `aPCA`, y = `B[0]`, fill = `B[0] res`), size = 4,
             shape = 21) +  geom_point(data = Hill.df[grep("^E36", rownames(Hill.df)),],
                                       aes(x = `aPCA`, y = `B[0]`, fill = `B[0] res`), size = 4,
                                       shape = 21) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = B0), linewidth = 2) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-30,30),
                        oob = squish) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-30,30),
                       oob = squish) +
  geom_segment(data = Hill.df[grep("^E36", rownames(Hill.df)),],
               aes(x = `aPCA`, y = `B[0]`,
                   xend = rep(120, length(grep("^E36", rownames(Hill.df)))),
                   yend = rep(10, length(grep("^E36", rownames(Hill.df))))),
               linewidth = 0.5,
               alpha = 0.2) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 130)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 120)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 120)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1 Abundance",
       y = bquote(B[0] ~ "(Relative PYL1/ABI1 Binding)"))

out.5A.Binf <- ggplot(Hill.df) +
  geom_point(data = Hill.df, 
             aes(x = `aPCA`, y = `B[inf]`, color = `B[inf] res`), size = 2, stroke = 0) +
  geom_point(data = Hill.df[grep("^D129", rownames(Hill.df)),], 
             aes(x = `aPCA`, y = `B[inf]`, fill = `B[inf] res`), size = 4,
             shape = 21) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = `Binf`), linewidth = 2) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-30,30),
                        oob = squish) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-30,30),
                       oob = squish) +
  geom_segment(data = Hill.df[grep("^D129", rownames(Hill.df)),],
               aes(x = `aPCA`, y = `B[inf]`,
                   xend = rep(100, length(grep("^D129", rownames(Hill.df)))),
                   yend = rep(115, length(grep("^D129", rownames(Hill.df))))),
               linewidth = 0.5,
               alpha = 0.2) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 130)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 120)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 120)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1 Abundance",
       y = bquote(B[infinity] ~ "(Relative PYL1/ABI1 Binding)"))

out.5A.EC50 <- ggplot(Hill.df) +
  geom_point(data = Hill.df, 
             aes(x = `aPCA`, y = `EC50`, color = `EC50 res`), size = 2, stroke = 0) +
  geom_point(data = Hill.df[grep("^H87|^E171", rownames(Hill.df)),],
             aes(x = `aPCA`, y = `EC50`, fill = `EC50 res`), size = 4,
             shape = 21) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = `EC50`), linewidth = 2) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-5.5, 5.5),
                        oob = squish) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-5.5, 5.5),
                       oob = squish) +
  geom_segment(data = Hill.df[grep("^H87", rownames(Hill.df)),],
               aes(x = `aPCA`, y = `EC50`,
                   xend = rep(median(Hill.df[grep("^H87", rownames(Hill.df)),"aPCA"]), length(grep("^H87", rownames(Hill.df)))),
                   yend = rep(0.05, length(grep("^H87", rownames(Hill.df))))),
               linewidth = 0.5,
               alpha = 0.2) +
  geom_segment(data = Hill.df[grep("^E171", rownames(Hill.df)),],
               aes(x = `aPCA`, y = `EC50`,
                   xend = rep(13.84615, length(grep("^E171", rownames(Hill.df)))),
                   yend = rep(0.05, length(grep("^E171", rownames(Hill.df))))),
               linewidth = 0.5,
               alpha = 0.2) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 130)) +
  coord_cartesian(xlim = c(-20, 130)) +
  scale_y_log10(breaks = c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000), 
                limits = c(0.01, 20000)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1 Abundance",
       y = bquote(EC[50] ~ "((+)-ABA, µM)"))

out.5A.n <- ggplot(Hill.df) +
  geom_point(data = Hill.df, 
             aes(x = `aPCA`, y = `Hill`, color = `Hill res`), size = 2, stroke = 0) +
  geom_point(data = Hill.df[grep("^Q42", rownames(Hill.df)),], 
             aes(x = `aPCA`, y = `Hill`, fill = `Hill res`), size = 4,
             shape = 21) +
  geom_line(data = loess.pred, aes(x = `aPCA`, y = `Hill`), linewidth = 2) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-0.7, 0.7),
                        oob = squish) +
  scale_fill_gradient2(low = "red", mid = "gray90", high = "blue", 
                       midpoint = 0, 
                       limits = c(-0.7, 0.7),
                       oob = squish) +
  geom_segment(data = Hill.df[grep("^Q42", rownames(Hill.df)),],
               aes(x = `aPCA`, y = `Hill`,
                   xend = rep(73.84615, length(grep("^Q42", rownames(Hill.df)))),
                   yend = rep(0.4, length(grep("^Q42", rownames(Hill.df))))),
               linewidth = 0.5,
               alpha = 0.2) +
  scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-20, 130)) +
  coord_cartesian(xlim = c(-20, 130)) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3), 
                limits = c(0.3, 3.5)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1 Abundance",
       y = "Hill coefficient (n)")

pdf('../../results/Figure5/Figure5A_Hill_vs_aPCA_individual_outliers.pdf', width = 9, height = 35)
plot_grid(out.5A.B0, out.5A.Binf, out.5A.EC50, out.5A.n, align = "hv", axis = "tblr", ncol = 1)
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
# [1] cowplot_1.2.0  ggrepel_0.9.6  ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 scales_1.4.0   stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4     S7_0.2.0           labeling_0.4.3     glue_1.8.0         plyr_1.8.9         gridtext_0.1.5    
# [11] grid_4.5.1         tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2       R6_2.6.1          
# [21] tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6       xml2_1.4.0