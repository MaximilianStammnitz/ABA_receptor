# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 22.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#############################################
## Figure 1A - Example dose-response curve ##
#############################################


## 0. Environment ##
####################

## Libraries
packages <- c("reshape", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Calculate ##
##################

## generic DRC function, given the four Hill parameters
four_pl_function <- function(x, Hill, G0, Ginf, EC50) {
  G0 + (Ginf - G0) / (1 + (EC50 / x)^c(Hill))
}

## Calculate the predicted response over the dose range
out <- four_pl_function(exp(seq(log(0.001), log(5000), length = 1000)), 
                        G0 = 2, Ginf = 98, EC50 = 1, Hill = 1)

curve <- cbind(out, "conc" = exp(seq(log(0.001), log(5000), length = 1000)))
curve <- as.data.frame(curve)


## 2. Plot ##
#############

## Calculate the WT dose-response curve
pdf("../../results/Figure1/Figure1A_example_dose_response_curve.pdf",
    height = 11, width = 25.5)

out.1A <- ggplot(curve[1:958,], aes(x = conc, y = out)) +
  geom_line(data = curve[1:958,], aes(x = conc, y = out), linewidth = 3) +
  geom_segment(aes(x = 0.001, xend = 2800, y = 0, yend = 0), 
              linetype = "dashed", size = 1.5, color = "black") +
  geom_segment(aes(x = 0.001, xend = 2800, y = 100, yend = 100), 
               linetype = "dashed", size = 1.5, color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = 0, yend = 100), 
               linetype = "dashed", size = 1.5, color = "black") +
  annotate("text", 
           x = 4900, 
           y = 0, 
           label = expr_text(bquote(italic(B[0]))),
           parse = T,
           color = "black",
           size = 23) +
  annotate("text", 
         x = 4900, 
         y = 100, 
         label = expr_text(bquote(italic(B[infinity]))),
         parse = T,
         color = "black",
         size = 23) +
  annotate("text", 
           x = 1, 
           y = 120, 
           label = expr_text(bquote(italic(EC[50]))),
           parse = T,
           vjust = 1, 
           color = "black",
           size = 23) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.001, 0.01, 0.1, 1, 10, 100, "1,000")) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 120)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 45),
        axis.line.x = element_line(size = 2, color = 'black'),
        axis.line.y = element_line(size = 2, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 75, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 75, vjust = 3),
        panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
        plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot background
        legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend background
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Small molecule conc. (µM)",
       y = "% Binding")

print(out.1A)

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
# [1] rlang_1.1.6    ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10
# 
# loaded via a namespace (and not attached):
# [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.3     gtable_0.3.6      
# [7] glue_1.8.0         tibble_3.3.0       pkgconfig_2.0.3    generics_0.1.4     dplyr_1.1.4        lifecycle_1.0.4   
# [13] xml2_1.4.0         cli_3.6.5          S7_0.2.0           scales_1.4.0       grid_4.5.1         vctrs_0.6.5       
# [19] withr_3.0.2        compiler_4.5.1     plyr_1.8.9         rstudioapi_0.17.1  tools_4.5.1        pillar_1.11.0     
# [25] Rcpp_1.1.0         crayon_1.5.3       gridtext_0.1.5
