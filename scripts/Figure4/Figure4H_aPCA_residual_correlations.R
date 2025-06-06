# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

############################################
## Figure 4H - aPCA residual correlations ##
############################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "ggplot2", "GGally")

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
                                 "log(EC50)" = loess.EC50$residuals,
                                 "log(Hill)" = loess.Hill$residuals))

rm(loess.B0, loess.Binf, loess.EC50, loess.Hill)


## 4. Plot ##
#############

# Define custom correlation display 
cor_func <- function(data, mapping, method, digits, display_grid, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method = method, use = 'complete.obs')
  ggally_text(
    label = as.expression(bquote(rho == .(formatC(corr, format = "f", digits = digits)))),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black', ...)
  
}

# Plotting function
pdf("../../results/Figure4/Figure4H_Hill_parameter_aPCA_res_corr.pdf", width = 10, height = 11)

out.4H <- ggpairs(data = loess.res, 
                  lower = list(continuous = wrap("points", size = 0.5, col = alpha("black", 0.05))),
                  upper = list(continuous = wrap(cor_func, size = 12, method = "spearman",
                                                 display_grid = F, digits = 2)),
                  diag = list(continuous = wrap("blank"))) +
  theme_test(base_size = 30) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank())

print(out.4H)

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
# [1] GGally_2.2.1  ggplot2_3.5.1 scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.4          rlang_1.1.5        stringi_1.8.4      purrr_1.0.4       
# [7] generics_0.1.3     labeling_0.4.3     glue_1.8.0         prettyunits_1.2.0  colorspace_2.1-1   plyr_1.8.9        
# [13] hms_1.1.3          grid_4.4.1         ggstats_0.8.0      munsell_0.5.1      tibble_3.2.1       progress_1.2.3    
# [19] lifecycle_1.0.4    compiler_4.4.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.0.14        pkgconfig_2.0.3   
# [25] tidyr_1.3.1        rstudioapi_0.17.1  farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.10.1     
# [31] magrittr_2.0.3     withr_3.0.2        tools_4.4.1        gtable_0.3.6      