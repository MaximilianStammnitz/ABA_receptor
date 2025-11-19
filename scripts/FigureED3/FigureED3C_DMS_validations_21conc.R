# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

####################################################################################
## Extended Data Figure 3C - PYL1-ABI1 Hill parameter validations (more DMS data) ##
####################################################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "drc", "ggplot2", "ggtext", 
              "cowplot", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)



## Libraries
packages <- c("stringr", "readxl", "growthcurver", "drc", 
              "scales", "reshape", "ggplot2", "ggtext", "cowplot",
              "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input original PYL1-ABI1 Hill parameters (12-concentration series) ##
###########################################################################

## input Hill parameter distributions from dose-response curve fits
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
main.Hill.params <- parameters.Hill[!is.na(parameters.Hill[,1]),]
main.Hill.params <- main.Hill.params[,c(2:4,1,10:12,9,ncol(main.Hill.params))]


## 2. Input new PYL1-ABI1 Hill parameters (21-concentration series) ##
######################################################################

## load data
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_21conc.RData")
extra.Hill.params <- parameters.Hill
extra.Hill.params <- extra.Hill.params[,c(2:4,1,10:12,9,ncol(extra.Hill.params))]
rm(parameters.Hill)


## 3. Combine ##
################

## build shared table
colnames(main.Hill.params) <- paste0("12-conc ", colnames(main.Hill.params))
colnames(extra.Hill.params) <- paste0("21-conc ", colnames(extra.Hill.params))
Hill.params <- extra.Hill.params
Hill.params <- cbind(Hill.params,
                     main.Hill.params[match(rownames(Hill.params), rownames(main.Hill.params)),])
class(Hill.params) <- "numeric"
Hill.params <- as.data.frame(Hill.params)

## regressions and plots
p.thresh <- 0.1
r2.thresh <- 0.5
point.size <- 3

### B0
out.B0 <- Hill.params[which(Hill.params$`21-conc B[0] P` < p.thresh & 
                            Hill.params$`21-conc R^2` > r2.thresh & 
                            Hill.params$`12-conc B[0] P` < p.thresh & 
                            Hill.params$`12-conc R^2` > r2.thresh),]
p.B0 <- summary(lm(out.B0$`21-conc B[0]` ~ out.B0$`12-conc B[0]`))$coefficients[2,4]

### Binf
out.Binf <- Hill.params[which(Hill.params$`21-conc B[inf] P` < p.thresh & 
                              Hill.params$`21-conc R^2` > r2.thresh & 
                              Hill.params$`12-conc B[inf] P` < p.thresh & 
                              Hill.params$`12-conc R^2` > r2.thresh),]
p.Binf <- summary(lm(out.Binf$`21-conc B[inf]` ~ out.Binf$`12-conc B[inf]`))$coefficients[2,4]

### EC50
out.EC50 <- Hill.params[which(Hill.params$`21-conc EC50 P` < p.thresh & 
                              Hill.params$`21-conc R^2` > r2.thresh & 
                              Hill.params$`12-conc EC50 P` < p.thresh & 
                              Hill.params$`12-conc R^2` > r2.thresh),]
p.EC50 <- summary(lm(log10(out.EC50$`21-conc EC50`) ~ log10(out.EC50$`12-conc EC50`)))$coefficients[2,4]

### Hill coefficient
out.Hill <- Hill.params[which(Hill.params$`21-conc Hill P` < p.thresh & 
                              Hill.params$`21-conc R^2` > r2.thresh & 
                              Hill.params$`12-conc Hill P` < p.thresh & 
                              Hill.params$`12-conc R^2` > r2.thresh),]
p.Hill <- summary(lm(log10(out.Hill$`21-conc Hill`) ~ log10(out.Hill$`12-conc Hill`)))$coefficients[2,4]

## B0
r.B0 <- cor(x = out.B0$`21-conc B[0]`, y = out.B0$`12-conc B[0]`, 
            method = "pearson", use = "complete.obs")

out.S3B_B0 <- ggplot(out.B0, aes(x = `21-conc B[0]`, y = `12-conc B[0]`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 130)) +
  geom_point(data = out.B0,
             mapping = aes(x = `21-conc B[0]`, y = `12-conc B[0]`),
             color = "black", size = point.size, shape = 16) +
  # geom_smooth(data = out.B0,
  #             mapping = aes(x = `21-conc B[0]`, y = `12-conc B[0]`),
  #             method = 'lm',
  #             color = "black", fullrange = T,
  #             linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = -20,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.B0, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.B0, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  # annotate("text",
  #          x = -20,
  #          y = 110,
  #          label = expr_text(bquote(italic(N) == .(format(nrow(out.B0))))),
  #          parse = T,
  #          hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[0] ~ "(21-conc. DMS series)"),
       y = bquote(B[0] ~ "(12-conc. DMS series)"))

## Binf
r.Binf <- cor(x = out.Binf$`21-conc B[inf]`, y = out.Binf$`12-conc B[inf]`, 
              method = "pearson", use = "complete.obs")

out.S3B_Binf <- ggplot(out.Binf, aes(x = `21-conc B[inf]`, y = `bulk_Binf`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 130)) +
  geom_point(data = out.Binf,
             mapping = aes(x = `21-conc B[inf]`, y = `12-conc B[inf]`),
             color = "black", size = point.size, shape = 16) +
  # geom_smooth(data = out.Binf,
  #             mapping = aes(x = `21-conc B[inf]`, y = `12-conc B[inf]`),
  #             method = 'lm',
  #             color = "black", fullrange = T,
  #             linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = -20,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.Binf, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.Binf, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  # annotate("text",
  #          x = -20,
  #          y = 110,
  #          label = expr_text(bquote(italic(N) == .(format(nrow(out.Binf))))),
  #          parse = T,
  #          hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[infinity] ~ "(21-conc. DMS series)"),
       y = bquote(B[infinity] ~ "(12-conc. DMS series)"))

## EC50
r.EC50 <- cor(x = log(out.EC50$`21-conc EC50`), y = log(out.EC50$`12-conc EC50`), 
              method = "pearson", use = "complete.obs")

out.S3B_EC50 <- ggplot(out.EC50, aes(x = `21-conc EC50`, y = `12-conc EC50`)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.001, 10000), ylim = c(0.001, 10000), expand = T) +
  geom_point(data = out.EC50,
             mapping = aes(x = `21-conc EC50`, y = `12-conc EC50`),
             color = "black", size = point.size, shape = 16) +
  # geom_smooth(data = out.EC50,
  #             mapping = aes(x = `21-conc EC50`, y = `12-conc EC50`),
  #             method = 'lm',
  #             color = "black", fullrange = T,
  #             linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = 0.0018,
           y = 3025,
           label = expr_text(bquote(italic(r) == .(format(r.EC50, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.EC50, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  # annotate("text",
  #          x = 0.1,
  #          y = 294,
  #          label = expr_text(bquote(italic(N) == .(format(nrow(out.EC50))))),
  #          parse = T,
  #          hjust = 0, size = 14, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(EC[50] ~ "(21-conc. DMS series)"),
       y = bquote(EC[50] ~ "(12-conc. DMS series)"))

## n
r.Hill <- cor(x = log(out.Hill$`21-conc Hill`), y = log(out.Hill$`12-conc Hill`), 
              method = "pearson", use = "complete.obs")

out.S3B_Hill <- ggplot(out.Hill, aes(x = `21-conc Hill`, y = `12-conc Hill`)) +
  scale_x_log10(breaks = c(0.4, 0.6, 1, 2, 4), 
                labels = c(0.4, 0.6, 1, 2, 4),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.4, 0.6, 1, 2, 4), 
                labels = c(0.4, 0.6, 1, 2, 4),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.4, 5), ylim = c(0.4, 5), expand = T) +
  geom_point(data = out.Hill,
             mapping = aes(x = `21-conc Hill`, y = `12-conc Hill`),
             color = "black", size = point.size, shape = 16) +
  # geom_smooth(data = out.Hill,
  #             mapping = aes(x = `21-conc Hill`, y = `12-conc Hill`),
  #             method = 'lm',
  #             color = "black", fullrange = T,
  #             linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = 0.4,
           y = 4.22,
           label = expr_text(bquote(italic(r) == .(format(r.Hill, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.Hill, digits = 4, nsmall = 4)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  # annotate("text",
  #          x = 0.4,
  #          y = 2.95,
  #          label = expr_text(bquote(italic(N) == .(format(nrow(out.Hill))))),
  #          parse = T,
  #          hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(n ~ "(21-conc. DMS series)"),
       y = bquote(n ~ "(12-conc. DMS series)"))

## combined plot
pdf("../../results/FigureED3/FigureED3C_Hill_parameter_DMS_validation.pdf", height = 15, width = 55)
plot_grid(out.S3B_B0, out.S3B_Binf, out.S3B_EC50, out.S3B_Hill, align = "hv", axis = "tblr", ncol = 4)
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
# [1] reshape_0.8.10     growthcurver_0.3.1 readxl_1.4.5       rlang_1.1.6        cowplot_1.2.0      ggtext_0.1.2       ggplot2_4.0.0      drc_3.0-1          MASS_7.3-65        scales_1.4.0      
# [11] stringr_1.5.2     
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      lattice_0.22-7     magrittr_2.0.4     grid_4.5.1         RColorBrewer_1.1-3 mvtnorm_1.3-3     
# [11] cellranger_1.1.0   plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    TH.data_1.1-4      codetools_0.2-20   abind_1.4-8        cli_3.6.5         
# [21] crayon_1.5.3       splines_4.5.1      withr_3.0.2        plotrix_3.8-4      tools_4.5.1        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4   
# [31] car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         tibble_3.3.0       tidyselect_1.2.1   farver_2.1.2       carData_3.0-5     
# [41] compiler_4.5.1     S7_0.2.0           gridtext_0.1.5   