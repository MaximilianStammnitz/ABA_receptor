# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 21.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################
## Figure 1G - bulk growth rate summmary ##
###########################################


## 0. Environment ##
####################

# ## Libraries
packages <- c("stringr", "scales", "ggplot2", "ggtext")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")


## 2. Plot ##
############

## Data frame to store density values
plot_data <- data.frame()
for (i in 1:12) {
  
  tmp.dist <- density(PYL1.ABI1[[i]]$gr_normalised_WTscaled, 
                      n = 1000, kernel = 'gaussian', 
                      from = -30, to = 130)
  if(i == 12){
    
    tmp_df <- data.frame(x = tmp.dist$x, 
                         y = tmp.dist$y + 0.01 * ((13 - i) / 12), 
                         group = factor(i, levels = 1:12), 
                         dosage = 9.062741e-03/3.5/3.5/3.5)
    
  }else{
    
    tmp_df <- data.frame(x = tmp.dist$x, 
                         y = tmp.dist$y + 0.01 * ((13 - i) / 12), 
                         group = factor(i, levels = 1:12), 
                         dosage = as.numeric(names(PYL1.ABI1)[i]))
    
  }

  plot_data <- rbind(plot_data, tmp_df)
  
}

## Define colors palette
grad_colors <- rev(alpha(colorRampPalette(c("darkgreen", "white"))(12), 0.8))

## Plot
pdf("../../results/Figure1/Figure1G_bulk_growth_rates.pdf", height = 15, width = 18)

out.1G <- ggplot(plot_data, aes(x = x, y = y, group = group, fill = dosage)) +
  geom_polygon(aes(group = group),
               color = "black",
               size = 0.5,
               alpha = 0.9,
               linewidth = 1) +
  scale_fill_gradientn(colors = grad_colors,
                       name = "(+)-ABA (µM)",
                       trans = "log",
                       breaks = c(0.001, 0.1, 10, 1000),
                       labels = c(0.001, 0.1, 10, "1,000"),
                       guide = guide_colorbar(reverse = F,
                                              barwidth = 40,
                                              barheight = 4,
                                              frame.colour = "black",
                                              frame.linewidth = 1.5,
                                              ticks.colour = "black",
                                              ticks.linewidth = 1.5)) +
  scale_x_continuous(limits = c(-30, 130),
                     breaks = seq(0, 100, length.out = 6),
                     labels = seq(0, 100, length.out = 6)) +
  coord_cartesian(xlim = c(-25, 125), ylim = c(0,0.05)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        legend.position = "top",
        legend.text = element_text(size = 60),
        legend.title = element_text(size = 70, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Binding (library sequencing)", y = "")

print(out.1G)

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
# [1] ggtext_0.1.2  ggplot2_4.0.0 scales_1.4.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4    
# [7] S7_0.2.0           labeling_0.4.3     glue_1.8.0         gridtext_0.1.5     grid_4.5.1         tibble_3.3.0      
# [13] lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3   
# [19] rstudioapi_0.17.1  farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.11.0      magrittr_2.0.3    
# [25] tools_4.5.1        withr_3.0.2        gtable_0.3.6       xml2_1.4.0  
