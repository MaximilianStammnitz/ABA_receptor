# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###################################
## Figure 3A - EC50 Volcano plot ##
###################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify sites with significantly lower EC50 ##
#####################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")

## remove stops
Hill.parameters.filtered <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]
Hill.parameters.filtered <- Hill.parameters.filtered[which(Hill.parameters.filtered[,"R^2"] > 0.9),]

## calculate log10-fold change between WT and variant
ec50_logFC <- log10(Hill.parameters.filtered[-1,"EC50"] / Hill.parameters.filtered["WT","EC50"])

## fetch B0, calculate z-tests
# B0.z_scores <- c(Hill.parameters.filtered[-1,"B[0]"] - Hill.parameters.filtered[1,"B[0]"]) / 
#   sqrt(Hill.parameters.filtered[-1,"B[0] SE"]^2 + Hill.parameters.filtered[1,"B[0] SE"]^2)
# p.B0 <- 2 * (1 - pnorm(abs(B0.z_scores)))
# p.adj.B0 <- p.adjust(p.B0, method = "BH")
# 
# ## fetch Binf, calculate z-tests
# Binf.z_scores <- c(Hill.parameters.filtered[-1,"B[inf]"] - Hill.parameters.filtered[1,"B[inf]"]) / 
#   sqrt(Hill.parameters.filtered[-1,"B[inf] SE"]^2 + Hill.parameters.filtered[1,"B[inf] SE"]^2)
# p.Binf <- 2 * (1 - pnorm(abs(Binf.z_scores)))
# p.adj.Binf <- p.adjust(p.Binf, method = "BH")

## fetch EC50, calculate z-tests
all_ec50 <- log10(Hill.parameters.filtered[-1,"EC50"])
all_ec50_se <- Hill.parameters.filtered[-1,"EC50 SE"] / c(Hill.parameters.filtered[-1,"EC50"] * log(10))
ref_ec50 <- log10(Hill.parameters.filtered["WT","EC50"])
ref_se <- Hill.parameters.filtered["WT","EC50 SE"] / c(Hill.parameters.filtered["WT","EC50"] * log(10))
ec50.z_scores <- (all_ec50 - ref_ec50) / sqrt(all_ec50_se^2 + ref_se^2)
p.ec50 <- 2 * (1 - pnorm(abs(ec50.z_scores)))
p.adj.ec50 <- p.adjust(p.ec50, method = "BH")

## fetch n, calculate z-tests
# all_n <- log10(Hill.parameters.filtered[-1,"Hill"])
# all_n_se <- Hill.parameters.filtered[-1,"Hill SE"] / c(Hill.parameters.filtered[-1,"Hill"] * log(10))
# ref_n <- log10(Hill.parameters.filtered["WT","Hill"])
# ref_n_se <- Hill.parameters.filtered["WT","Hill SE"] / c(Hill.parameters.filtered["WT","Hill"] * log(10))
# n.z_scores <- (all_n - ref_n) / sqrt(all_n_se^2 + ref_n_se^2)
# p.n <- 2 * (1 - pnorm(abs(n.z_scores)))
# p.adj.n <- p.adjust(p.n, method = "BH")

## summarise data
data.out <- as.data.frame(cbind("logFC" = ec50_logFC,
                                "p" = p.ec50,
                                "p.adj" = p.adj.ec50))
data.out$label <- names(ec50_logFC)
data.out$p.adj[data.out$p.adj == 0] <- 1e-15

# data.out <- as.data.frame(cbind("p.adj.B0" = p.adj.B0,
#                                 "p.adj.Binf" = p.adj.B0,
#                                 "p.adj.ec50" = p.adj.ec50,
#                                 "p.adj.n" = p.adj.n))
# length(which(data.out$p.adj.B0 < 0.01 | data.out$p.adj.Binf < 0.01 | data.out$p.adj.ec50 < 0.01 | data.out$p.adj.n < 0.01)) / nrow(data.out)

## 2. Display ##
################

## highlight selected hits
labs.out <- subset(data.out, `p.adj` < 0.1 & `logFC` < 0)
labs.out <- labs.out[order(labs.out$logFC),]
labs.out <- data.out[grep("T33V|T33E|
                           |H87P|H87V|H87R|H87Q|H87C|H87A|
                           |D107V|D107I|D107T|
                           |N117K|
                           |T118N|
                           |L125I|
                           |E141D|
                           |R164S|R164N|R164D|R164E|
                           |I165E|I165P|
                           |W166M|W166T|W166C|
                           |E171T|
                           |D185E|
                           |L188T|L188C|
                           |A190V|A190T|
                           |V193I|
                           |L196T|L196C|L196M|
                           |Q199R|Q199K|
                           |K200S|K200T|
                           |S203E", rownames(data.out)),]

pdf("../../results/Figure3/Figure3A_Volcano_EC50.pdf", width = 7, height = 8)

out.3A <- ggplot(data.out, aes(x = logFC, y = -log10(p.adj))) +
  geom_point(data = subset(data.out, !(p.adj < 0.1 & logFC < 0)),
             aes(x = logFC, y = -log10(p.adj), color = squish(logFC, range = c(-4.5, 4.5))),
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = subset(data.out, p.adj < 0.1 & logFC < 0),
             aes(x = logFC, y = -log10(p.adj), fill = squish(logFC, range = c(-4.5, 4.5))),
             shape = 21,
             color = "black",
             size = 3,
             stroke = 0.5,
             alpha = 1) +
  scale_colour_gradient2(low = "blue", mid = "gray90", high = "red", 
                         midpoint = 0, 
                         limits = c(-2,2),
                         oob = squish) +
  scale_fill_gradient2(low = "blue", mid = "gray90", high = "red", 
                       midpoint = 0, 
                       limits = c(-2,2),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.1), 
             linetype = "dashed", 
             color = "black", alpha = 0.2, linewidth = 0.5) +
  geom_text_repel(data = labs.out,
                  aes(label = label),
                  size = 4, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 2,
                  force_pull = 0.5,
                  box.padding = 0.85,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = seq(f = -5, t = 5, by = 1), 
                     limits = c(-2, 1),
                     labels = c("A", "B", "C", "0.01", "0.1", "1", "10", "100", "I", "J", "K")) +
  scale_y_continuous(breaks = seq(f = 0, t = 15, by = 3), 
                     limits = c(0, 15), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = 18),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 25, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 25, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = bquote(EC[50] * " fold change"),
       y = bquote(-log[10] * "(adjusted P)"))

print(out.3A)

dev.off()


## 3. Version ##
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
# [1] ggrepel_0.9.6  ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 scales_1.4.0   stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        stringi_1.8.7      generics_0.1.4     S7_0.2.0           labeling_0.4.3     glue_1.8.0         plyr_1.8.9        
# [11] gridtext_0.1.5     grid_4.5.1         tibble_3.3.0       lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2      
# [21] R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.2        tools_4.5.1        gtable_0.3.6       xml2_1.4.0  