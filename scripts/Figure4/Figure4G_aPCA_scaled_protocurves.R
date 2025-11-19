# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#########################################
## Figure 4G - aPCA scaled protocurves ##
#########################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "viridis")

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

## input Hill parameter distributions from dose-response curve fits, filter like in Figure 2
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


## 3. Calibrate protocurves ##
##############################

## LOESS of all four parameters vs. pseudo-abundance
loess.pred <- expand.grid("aPCA" = seq(f = min(parameters.Hill[,"aPCA"]), 
                                       to = max(parameters.Hill[,"aPCA"]), 
                                       length.out = 50))
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric") ## log?
loess.pred$B0 <- as.numeric(predict(loess.B0, newdata = loess.pred, se = TRUE)$fit)
loess.pred$Binf <- as.numeric(predict(loess.Binf, newdata = loess.pred, se = TRUE)$fit)
loess.pred$EC50 <- exp(as.numeric(predict(loess.EC50, newdata = loess.pred, se = TRUE)$fit))
loess.pred$Hill <- exp(as.numeric(predict(loess.Hill, newdata = loess.pred, se = TRUE)$fit)) ## log?
rm(loess.B0, loess.Binf, loess.EC50, loess.Hill)

## calculate "proto"-curves for an equally distributed range of 100 abundances 

### Define the four-parameter logistic function
four_pl_function <- function(x, Hill, B0, Binf, EC50) {
  B0 + (Binf - B0) / (1 + (EC50 / x)^c(Hill))
}

### Calculate the predicted response over the dose range

### plot these "proto"-curves
protocurves.all <- matrix(NA, nrow = 50, ncol = 1000)
rownames(protocurves.all) <- loess.pred$aPCA
for(i in 1:nrow(loess.pred)){
  
  print(i)
  protocurves.all[i,] <- four_pl_function(exp(seq(log(0.001), log(5000), length = 1000)), 
                                          Hill = loess.pred[i,"Hill"], 
                                          B0 = loess.pred[i,"B0"],
                                          Binf = loess.pred[i,"Binf"],
                                          EC50 = loess.pred[i,"EC50"])
  
}
protocurves.all <- rbind(protocurves.all,
                         "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
protocurves.all <- t(protocurves.all)
protocurves.all <- as.data.frame(protocurves.all)
protocurves.lines <- melt(protocurves.all, id.vars = "conc")
protocurves.lines$variable <- as.numeric(as.character(protocurves.lines$variable))


## 4. Plot ##
#############

### assign the right colour for each curve
cols <- magma(n = 1000)
names(cols) <- seq(f = 0, to = 130, length.out = 1000)
closest.col <- sapply(as.numeric(colnames(protocurves.all)[-101]), function(xi) {cols[which.min(abs(as.numeric(names(cols)) - xi))]})
names(closest.col) <- colnames(protocurves.all)[-101]
protocurves.lines$col <- closest.col[match(as.character(protocurves.lines$variable), names(closest.col))]

pdf('../../results/Figure4/Figure4G_aPCA_scaled_protocurves.pdf', height = 15, width = 23)

out.4G <- ggplot(protocurves.lines, aes(x = conc, y = value)) +
  geom_line(aes(x = conc, y = value, group = variable, colour = variable), size = 1.5) +
  scale_colour_gradientn(colours = magma(n = 1000), 
                         limits = c(0, 130),
                         breaks = seq(from = 0, to = 100, length.out = 6),
                         oob = scales::squish) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.1, 10, 1000),
                labels = c(0, 0.1, 10, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        axis.text = element_text(size = 45),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
        legend.title = element_text(size = 35, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 30, family = "Helvetica"),
        legend.key.size = unit(3, "lines"),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "(+)-ABA concentration (µM)",
       y = "Relative PYL1/ABI1 Binding (expected)",
       color = "Abundance")

print(out.4G)

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
# [1] viridis_0.6.5     viridisLite_0.4.2 ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4     S7_0.2.0           glue_1.8.0         plyr_1.8.9         gridExtra_2.3     
# [11] gridtext_0.1.5     grid_4.5.1         tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2      
# [21] R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6       xml2_1.4.0  