# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##################################################################################
## Extended Data Figure 3D - correlation of Elzinga PYR1 mutants vs. PYL1 EC50s ##
##################################################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "readxl", "growthrates", "drc", 
              "scales", "reshape", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify sites with significantly lower EC50 ##
#####################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")


## 2. Load Elzinga et al. single mutations (ACS ChemBio 2019) ##
################################################################

## note that the exact system here is:
## PYR1 (not PYL1)
## HAB1 (not ABI1)
## Y2M (not DHFR complementation)
pyr1.muts.elzinga <- matrix(NA, nrow = 9, ncol = 2)
rownames(pyr1.muts.elzinga) <- c("WT", "A160C", "A160V", "F61L",
                                 "F61M", "I110C", "I110S", "V81I",
                                 "E141L")
colnames(pyr1.muts.elzinga) <- c("EC50", "SD")
pyr1.muts.elzinga[,"EC50"] <- c(0.864, 0.167, 0.0982, 0.613,
                                2.98, 0.873, 1.15, 0.242,
                                0.238)
pyr1.muts.elzinga[,"SD"] <- c(0.061, 0.022, 0.011, 0.079,
                              1.1, 0.17, 0.28, 0.242,
                              0.238)

## translate the mutants to PYL1 positions
rownames(pyr1.muts.elzinga)[-1] <- paste0(substr(rownames(pyr1.muts.elzinga)[-1], 1, 1),
                                          as.numeric(substr(rownames(pyr1.muts.elzinga)[-1], 2, nchar(rownames(pyr1.muts.elzinga)[-1]) - 1)) + 30,
                                          substr(rownames(pyr1.muts.elzinga)[-1], nchar(rownames(pyr1.muts.elzinga)[-1]), nchar(rownames(pyr1.muts.elzinga)[-1])))
rownames(pyr1.muts.elzinga)[4] <- "F88L"
rownames(pyr1.muts.elzinga)[5] <- "F88M"
rownames(pyr1.muts.elzinga)[6] <- "I137C"
rownames(pyr1.muts.elzinga)[7] <- "I137S"
rownames(pyr1.muts.elzinga)[8] <- "V108I"


## 3. Combine and plot ##
#########################

data.out <- pyr1.muts.elzinga
data.out <- cbind(data.out,
                  parameters.Hill[match(rownames(data.out), rownames(parameters.Hill)),c(4,8)])
colnames(data.out) <- c("Elzinga EC50", "Elzinga EC50_SD", "DMS EC50", "DMS EC50_SE")
data.out <- as.data.frame(data.out)

## correlation
r.EC50 <- cor(x = log(data.out$`Elzinga EC50`),  y = log(data.out$`DMS EC50`), method = "pearson") ## 0.61
p.EC50 <- summary(lm(log10(data.out$`Elzinga EC50`) ~ log10(data.out$`DMS EC50`)))$coefficients[2,4] ## P = 0.07926899


pdf('../../results/FigureED3/FigureED3D_Elzinga_correlation.pdf', height = 15, width = 15)

ggplot(data.out, aes(x = `Elzinga EC50`, y = `DMS EC50`)) +
  scale_x_log10(breaks = c(0.1, 0.1, 1, 10), 
                labels = c(0.1, 0.1, 1, 10),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10), 
                labels = c(0.01, 0.1, 1, 10),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.1, 75), ylim = c(0.1, 75), expand = T) +
  geom_smooth(data = data.out,
              mapping = aes(x = `Elzinga EC50`, 
                            y = `DMS EC50`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_point(aes(x = `Elzinga EC50`, y = `DMS EC50`), 
             pch = 16, 
             col = "black", size = 9) +
  geom_errorbar(aes(x = `Elzinga EC50`,
                    ymin = `DMS EC50` - `DMS EC50_SE`, ymax = `DMS EC50` + `DMS EC50_SE`),
                color = "black", 
                linewidth = 0.75, width = 0) +
  geom_errorbar(aes(y = `DMS EC50`,
                    xmin = `Elzinga EC50` - `Elzinga EC50_SD`, xmax = `Elzinga EC50` + `Elzinga EC50_SD`),
                color = "black", 
                linewidth = 0.75, height = 0) +
  annotate("text",
           x = 0.1,
           y = 60,
           label = expr_text(bquote(italic(r) == .(format(r.EC50, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.EC50, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 2),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = expression(atop("PYR1-HAB1 "*EC[50]* " (µM)", "(Elzinga et al., 2019)")),
       y = expression(atop("PYL1-ABI1 "*EC[50]* " (µM)", "(this study)")))

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
# [1] rlang_1.1.6       ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    scales_1.4.0      drc_3.0-1         MASS_7.3-65       growthrates_0.8.5 deSolve_1.40      lattice_0.22-7   
# [11] readxl_1.4.5      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      FME_1.3.6.4        magrittr_2.0.4     grid_4.5.1         RColorBrewer_1.1-3 mvtnorm_1.3-3     
# [11] cellranger_1.1.0   plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    mgcv_1.9-3         TH.data_1.1-4      codetools_0.2-20   abind_1.4-8       
# [21] cli_3.6.5          crayon_1.5.3       splines_4.5.1      withr_3.0.2        plotrix_3.8-4      rootSolve_1.8.2.4  tools_4.5.1        parallel_4.5.1     coda_0.19-4.1      minpack.lm_1.2-4  
# [31] minqa_1.2.8        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4    car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
# [41] glue_1.8.0         Rcpp_1.1.0         tidyselect_1.2.1   tibble_3.3.0       farver_2.1.2       nlme_3.1-168       carData_3.0-5      compiler_4.5.1     S7_0.2.0           gridtext_0.1.5 