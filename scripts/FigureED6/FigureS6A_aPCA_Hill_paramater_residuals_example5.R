# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##################################################################
## Supplementary Figure 6A - aPCA residual outliers (example 5) ##
##################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel", "cowplot")

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
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
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
loess.pred$Hill <- exp(as.numeric(predict(loess.Hill, newdata = loess.pred, se = TRUE)$fit))

## summarise data
Hill.df <- as.data.frame(cbind(parameters.Hill[,c(2,3,4,1,17)],
                               "B[0] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"B0"],
                               "B[inf] res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Binf"],
                               "EC50 res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"EC50"],
                               "Hill res" = loess.res[match(rownames(parameters.Hill), rownames(loess.res)),"Hill"]))
Hill.df$res <- as.factor(c("WT", substr(Hill.df$names, 1, nchar(Hill.df$names)-1)[-1]))


## 4. Plots ##
##############

## B0
B0.res <- c("E36", "T33", "S43", "Q39", "F37")
p.out.B0 <- vector(mode = "list", length = 5)

for(i in 1:5){
  
  p.out.B0[[i]] <- ggplot(Hill.df) +
    geom_point(data = Hill.df, 
               aes(x = `aPCA`, y = `B[0]`, color = `B[0] res`), size = 1) +
    geom_point(data = Hill.df[grep(B0.res[i], rownames(Hill.df)),],
               aes(x = `aPCA`, y = `B[0]`, fill = `B[0] res`), size = 6,
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
    scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 120)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 35),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "Abundance",
         y = bquote(B[0] ~ "(Binding)"),
         title = paste0(B0.res[i], " (mean ", format(round(mean(Hill.df[grep(B0.res[i], rownames(Hill.df)),"B[0] res"]),1), nsmall = 1), "% residual)"))
  
}

## Binf
Binf.res <- c("R130", "R164", "D129", "E162", "E160")
p.out.Binf <- vector(mode = "list", length = 5)

for(i in 1:5){
  
  p.out.Binf[[i]] <- ggplot(Hill.df) +
    geom_point(data = Hill.df, 
               aes(x = `aPCA`, y = `B[inf]`, color = `B[inf] res`), size = 1) +
    geom_point(data = Hill.df[grep(Binf.res[i], rownames(Hill.df)),], 
               aes(x = `aPCA`, y = `B[inf]`, fill = `B[inf] res`), size = 6,
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
    scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 120)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 35),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "Abundance",
         y = bquote(B[infinity] ~ "(Binding)"),
         title = paste0(Binf.res[i], " (mean ", format(round(mean(Hill.df[grep(Binf.res[i], rownames(Hill.df)),"B[inf] res"]),1), nsmall = 1), "% residual)"))
  
}

## EC50
EC50.res <- c("L196", "H87", "E171", "K200", "Y50")
p.out.EC50 <- vector(mode = "list", length = 5)

for(i in 1:5){
  
  p.out.EC50[[i]] <- ggplot(Hill.df) +
    geom_point(data = Hill.df, 
               aes(x = `aPCA`, y = `EC50`, color = `EC50 res`), size = 1) +
    geom_point(data = Hill.df[grep(EC50.res[i], rownames(Hill.df)),],
               aes(x = `aPCA`, y = `EC50`, fill = `EC50 res`), size = 6,
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
    scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    coord_cartesian(xlim = c(-5, 120)) +
    scale_y_log10(breaks = c(0.01, 0.1, 0, 1, 10, 100, 1000, 10000), 
                  limits = c(0.01, 20000)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 35),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "Abundance",
         y = bquote(EC[50] ~ "((+)-ABA, µM)"),
         title = paste0(EC50.res[i], " (mean ", format(round(exp(mean(Hill.df[grep(EC50.res[i], rownames(Hill.df)),"EC50 res"])),3), nsmall = 3), "× residual)"))
  
}

## n
n.res <- c("Q42", "H154", "S41", "P68", "S43")
p.out.n <- vector(mode = "list", length = 5)

for(i in 1:5){
  
  p.out.n[[i]] <- ggplot(Hill.df) +
    geom_point(data = Hill.df, 
               aes(x = `aPCA`, y = `Hill`, color = `Hill res`), size = 1) +
    geom_point(data = Hill.df[grep(n.res[i], rownames(Hill.df)),], 
               aes(x = `aPCA`, y = `Hill`, fill = `Hill res`), size = 6,
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
    scale_x_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
    coord_cartesian(xlim = c(-5, 130)) +
    scale_y_log10(breaks = c(0.5, 1, 2, 3), 
                  limits = c(0.3, 3.5)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 35),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "Abundance",
         y = "Hill coefficient (n)",
         title = paste0(n.res[i], " (mean ", format(round(exp(mean(Hill.df[grep(n.res[i], rownames(Hill.df)),"Hill res"])),3), nsmall = 3), "× residual)"))
  
}

## plot all
pdf('../../results/FigureS6/FigureS6A_Hill_paramater_residuals_example5.pdf', width = 60, height = 42)

plot_grid(p.out.B0[[1]], p.out.B0[[2]], p.out.B0[[3]], p.out.B0[[4]], p.out.B0[[5]], 
          p.out.Binf[[1]], p.out.Binf[[2]], p.out.Binf[[3]], p.out.Binf[[4]], p.out.Binf[[5]], 
          p.out.EC50[[1]], p.out.EC50[[2]], p.out.EC50[[3]], p.out.EC50[[4]], p.out.EC50[[5]], 
          p.out.n[[1]], p.out.n[[2]], p.out.n[[3]], p.out.n[[4]], p.out.n[[5]], 
          align = "hv", axis = "tblr", ncol = 5)

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
# [1] cowplot_1.1.3 ggrepel_0.9.6 ggtext_0.1.2  ggplot2_3.5.1 reshape_0.8.9 scales_1.3.0  stringr_1.5.1 readxl_1.4.3 
# 
# loaded via a namespace (and not attached):
# [1] gtable_0.3.6      dplyr_1.1.4       compiler_4.4.1    tidyselect_1.2.1  Rcpp_1.0.14       xml2_1.3.6        R6_2.6.1         
# [8] plyr_1.8.9        commonmark_1.9.2  labeling_0.4.3    generics_0.1.3    tibble_3.2.1      munsell_0.5.1     pillar_1.10.1    
# [15] rlang_1.1.5       stringi_1.8.4     xfun_0.51         cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.4.1       
# [22] gridtext_0.1.5    rstudioapi_0.17.1 markdown_1.13     lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2     
# [29] cellranger_1.1.0  colorspace_2.1-1  tools_4.4.1       pkgconfig_2.0.3  