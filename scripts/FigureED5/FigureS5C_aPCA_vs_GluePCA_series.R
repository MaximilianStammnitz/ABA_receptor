# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##############################################################################################
## Supplementary Figure 5C - LOESS fit of aPCA vs. all 12 GluePCA concentrations of (+)-ABA ##
##############################################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel", "cowplot")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input aPCA data and GluePCA Hill parameters ##
####################################################

## aPCA
aPCA <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S2.xlsx", sheet = 2))[-c(1:3),1:2]
colnames(aPCA) <- c("Variant", "aPCA")
rownames(aPCA) <- aPCA[,1]
aPCA <- aPCA[,-1,drop=F]
class(aPCA) <- "numeric"

## GluePCA
GluePCA <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S2.xlsx", sheet = 1))[-c(1:2),]
colnames(GluePCA) <- GluePCA[1,]
colnames(GluePCA)[1] <- "Variant"
GluePCA <- GluePCA[-1,]
rownames(GluePCA) <- GluePCA[,1]
GluePCA <- GluePCA[,-1]
class(GluePCA) <- "numeric"
colnames(GluePCA) <- 1:12


## 2. Combine data ##
#####################

all.out <- as.data.frame(cbind(aPCA, GluePCA[match(rownames(aPCA), rownames(GluePCA)),]))

## clean up
rm(aPCA, GluePCA, packages)


## 3. Fit LOESS residuals ##
############################

## ABA dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
dosages <- rev(dosages)

## LOESS of all twelve comparisons (using "alive" aPCA values >25)
all.out.alive <- all.out[which(all.out$aPCA >25),]
loess.1 <- loess(`1` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.2 <- loess(`2` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.3 <- loess(`3` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.4 <- loess(`4` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.5 <- loess(`5` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.6 <- loess(`6` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.7 <- loess(`7` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.8 <- loess(`8` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.9 <- loess(`9` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.10 <- loess(`10` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.11 <- loess(`11` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")
loess.12 <- loess(`12` ~ aPCA, all.out.alive, span = 0.5, family = "symmetric")

## predict full-range LOESS curves
loess.pred <- expand.grid("aPCA" = seq(f = 25, 
                                       to = range(all.out.alive$aPCA, na.rm = T)[2], 
                                       length.out = 1000))
loess.pred$`1` <- as.numeric(predict(loess.1, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`2` <- as.numeric(predict(loess.2, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`3` <- as.numeric(predict(loess.3, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`4` <- as.numeric(predict(loess.4, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`5` <- as.numeric(predict(loess.5, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`6` <- as.numeric(predict(loess.6, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`7` <- as.numeric(predict(loess.7, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`8` <- as.numeric(predict(loess.8, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`9` <- as.numeric(predict(loess.9, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`10` <- as.numeric(predict(loess.10, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`11` <- as.numeric(predict(loess.11, newdata = loess.pred$aPCA, se = TRUE)$fit)
loess.pred$`12` <- as.numeric(predict(loess.12, newdata = loess.pred$aPCA, se = TRUE)$fit)

## LOESS R^2
R_2.out <- vector(mode = "list", length = 12)
R_2.out[[1]] <- 1 - (sum((loess.1$residuals)^2) / sum((all.out.alive$`1` - mean(all.out.alive$`1`, na.rm = T))^2, na.rm = T))
R_2.out[[2]] <- 1 - (sum((loess.2$residuals)^2) / sum((all.out.alive$`2` - mean(all.out.alive$`2`, na.rm = T))^2, na.rm = T))
R_2.out[[3]] <- 1 - (sum((loess.3$residuals)^2) / sum((all.out.alive$`3` - mean(all.out.alive$`3`, na.rm = T))^2, na.rm = T))
R_2.out[[4]] <- 1 - (sum((loess.4$residuals)^2) / sum((all.out.alive$`4` - mean(all.out.alive$`4`, na.rm = T))^2, na.rm = T))
R_2.out[[5]] <- 1 - (sum((loess.5$residuals)^2) / sum((all.out.alive$`5` - mean(all.out.alive$`5`, na.rm = T))^2, na.rm = T))
R_2.out[[6]] <- 1 - (sum((loess.6$residuals)^2) / sum((all.out.alive$`6` - mean(all.out.alive$`6`, na.rm = T))^2, na.rm = T))
R_2.out[[7]] <- 1 - (sum((loess.7$residuals)^2) / sum((all.out.alive$`7` - mean(all.out.alive$`7`, na.rm = T))^2, na.rm = T))
R_2.out[[8]] <- 1 - (sum((loess.8$residuals)^2) / sum((all.out.alive$`8` - mean(all.out.alive$`8`, na.rm = T))^2, na.rm = T))
R_2.out[[9]] <- 1 - (sum((loess.9$residuals)^2) / sum((all.out.alive$`9` - mean(all.out.alive$`9`, na.rm = T))^2, na.rm = T))
R_2.out[[10]] <- 1 - (sum((loess.10$residuals)^2) / sum((all.out.alive$`10` - mean(all.out.alive$`10`, na.rm = T))^2, na.rm = T))
R_2.out[[11]] <- 1 - (sum((loess.11$residuals)^2) / sum((all.out.alive$`11` - mean(all.out.alive$`11`, na.rm = T))^2, na.rm = T))
R_2.out[[12]] <- 1 - (sum((loess.12$residuals)^2) / sum((all.out.alive$`12` - mean(all.out.alive$`12`, na.rm = T))^2, na.rm = T))

## define colours
grad_colors <- alpha(colorRampPalette(c("grey75", "darkgreen"))(12), 0.8)

## define vector of plots
plots.out <- vector(mode = "list", length = 12)
for (i in 1:12){
  
  name.tmp <- colnames(all.out)[i+1]
  
  plots.out[[i]] <- ggplot(all.out, mapping = aes(x = `aPCA`, y = .data[[name.tmp]])) +
    geom_point(shape = 16,
               color = grad_colors[i], 
               size = 1.5) +
    geom_line(data = loess.pred, aes(x = `aPCA`, y = .data[[name.tmp]]), linewidth = 2) +
    scale_x_continuous(breaks = seq(from = 0, to = 120, length.out = 7), limits = c(-5, 120)) +
    scale_y_continuous(breaks = seq(from = 0, to = 120, length.out = 7), limits = c(-5, 120)) +
    coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(hjust = 0.5),
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
         y = "Binding",
         title = paste0(format(round(dosages[i], 3), nsmall = 3), " µM (+)-ABA"))
  
}

pdf('../../results/FigureS5/FigureS5C_abundancePCA_vs_binding_series.pdf', width = 58, height = 20)

plot_grid(plots.out[[1]], plots.out[[2]], plots.out[[3]], plots.out[[4]], plots.out[[5]], plots.out[[6]], 
          plots.out[[7]], plots.out[[8]], plots.out[[9]], plots.out[[10]], plots.out[[11]], plots.out[[12]], align = "hv", axis = "tblr", ncol = 6)

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
# [1] cowplot_1.1.3 ggrepel_0.9.6 ggtext_0.1.2  ggplot2_3.5.1 reshape_0.8.9 scales_1.3.0  stringr_1.5.1 readxl_1.4.3 
# 
# loaded via a namespace (and not attached):
# [1] gtable_0.3.6      dplyr_1.1.4       compiler_4.4.1    tidyselect_1.2.1  Rcpp_1.0.14       xml2_1.3.6        R6_2.6.1         
# [8] plyr_1.8.9        commonmark_1.9.2  generics_0.1.3    tibble_3.2.1      munsell_0.5.1     pillar_1.10.1     rlang_1.1.5      
# [15] stringi_1.8.4     xfun_0.51         cli_3.6.4         withr_3.0.2       magrittr_2.0.3    grid_4.4.1        gridtext_0.1.5   
# [22] rstudioapi_0.17.1 markdown_1.13     lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        farver_2.1.2      cellranger_1.1.0 
# [29] colorspace_2.1-1  tools_4.4.1       pkgconfig_2.0.3  