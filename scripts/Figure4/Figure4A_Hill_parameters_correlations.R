# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################################
## Figure 4A - Hill parameters, pairwise correlations ##
########################################################


## 0. Environment ##
####################

# # Libraries
packages <- c("stringr", "scales", "ggplot2", "GGally", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input and filter dose response curve parameters ##
########################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## filter like in Figure 2
parameters.Hill <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.9),]
parameters.Hill <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]

## subset, log(EC50), log(n)
corr.in <- parameters.Hill[,c(2:4,1)]
corr.in[,"EC50"] <- log(corr.in[,"EC50"])
corr.in[,"Hill"] <- log(corr.in[,"Hill"])
colnames(corr.in)[1:4] <- c("B[0]", "B[inf]", "log(EC50)", "log(n)")
corr.in <- as.data.frame(corr.in)


## 2. Plot parameter correlations ##
####################################

# Define custom correlation display 
cor_func <- function(data, mapping, method, digits, display_grid, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  corr <- cor(x, y, method = method, use = 'complete.obs')
  ggally_text(
    label = expr_text(bquote(rho == .(formatC(corr, format = "f", digits = digits)))),
    parse = T,
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black', ...)
  
}

# Plotting function
pdf("../../results/Figure4/Figure4A_Hill_parameter_corr.pdf", width = 10, height = 11)

out.4A <- ggpairs(data = corr.in, 
                  lower = list(continuous = wrap("points", size = 0.01, col = alpha("black", 0.1), shape = 16)),
                  upper = list(continuous = wrap(cor_func, size = 12, method = "spearman",
                                                 display_grid = F, digits = 2)),
                  diag = list(continuous = wrap("blank"))) +
  theme_test(base_size = 30) + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_blank())

print(out.4A)

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
# [1] rlang_1.1.6   GGally_2.4.0  ggplot2_4.0.0 scales_1.4.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          stringi_1.8.7      purrr_1.1.0       
# [6] generics_0.1.4     S7_0.2.0           labeling_0.4.3     glue_1.8.0         prettyunits_1.2.0 
# [11] hms_1.1.3          grid_4.5.1         ggstats_0.11.0     tibble_3.3.0       progress_1.2.3    
# [16] lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 pkgconfig_2.0.3   
# [21] tidyr_1.3.1        rstudioapi_0.17.1  farver_2.1.2       R6_2.6.1           tidyselect_1.2.1  
# [26] pillar_1.11.0      magrittr_2.0.3     tools_4.5.1        withr_3.0.2        gtable_0.3.6  
