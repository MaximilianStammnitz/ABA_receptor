# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 3D - EC50 distance decay to ABI1 ##
###########################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "cowplot")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Calculate log EC50 fold change ##
#######################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## remove stops
Hill.parameters.filtered <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]

## calculate log10-fold change between WT and variant
ec50_log10FC <- log10(Hill.parameters.filtered[-1,"EC50"] / Hill.parameters.filtered["WT","EC50"])


## 2. Annotate residue-specific distances/contacts to (+)-ABA and ABI1 ##
#########################################################################

## input residue minimal distances from ABA and ABI1
load("../../data/PYL1_distances_to_interfaces.Rdata")

ABA_dist <- PYL1_ABA_ABI1_dist$ABA$ABA_vs_PYL1_HAmin_ligand[-c(1:2)]
ABI1_dist <- PYL1_ABA_ABI1_dist$ABI1$ABI1_vs_PYL1_HAmin_ligand[-c(1:2)]
names(ABA_dist) <- names(ABI1_dist) <- 33:209

## add to data summary
data.all <- as.data.frame(cbind("EC50_logFC" = ec50_log10FC))
data.all$ABA_dist <- ABA_dist[match(as.numeric(substr(rownames(data.all), 2, nchar(rownames(data.all)) - 1)),
                                    as.numeric(names(ABA_dist)))]
data.all$ABI1_dist <- ABI1_dist[match(as.numeric(substr(rownames(data.all), 2, nchar(rownames(data.all)) - 1)),
                                      as.numeric(names(ABI1_dist)))]

## annotate (+)-ABA and ABI1 contacts (Miyazono et al., Nature 2009)
PYL1.ABA.contact <- PYL1.ABI1.contact <- rep(0, 177)
names(PYL1.ABA.contact) <- names(PYL1.ABI1.contact) <- 33:209
PYL1.ABA.contact[match(c(86,88,110,116,121,135,137,143,144,147,149,171,189,193,197), names(PYL1.ABA.contact))] <- 1
PYL1.ABI1.contact[match(c(87,88,90,111:117,142:144,178,180,181,185,186,188,189,192,193,196,199), names(PYL1.ABI1.contact))] <- 1
names(PYL1.ABA.contact) <- names(PYL1.ABI1.contact) <- paste0(c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177)) ,33:209)
data.all$ABA_contact <- PYL1.ABA.contact[match(substr(rownames(data.all), 1, nchar(rownames(data.all)) - 1),names(PYL1.ABA.contact))]
data.all$ABI1_contact <- PYL1.ABI1.contact[match(substr(rownames(data.all), 1, nchar(rownames(data.all)) - 1),names(PYL1.ABI1.contact))]


## 3. Calculate decay parameters ##
###################################

## partition data into positive and negative ∆∆Gb
data.all.posFC <- data.all[which(data.all$EC50_logFC >= 0),]
data.all.negFC <- data.all[which(data.all$EC50_logFC < 0),]

## exponential decay function
exp_decay <- function(x, a, b) {
  a * exp(-b * x)
}

## calculate exponential distance decay regressions
PYL1.posEC50_logFC.ABI1.exp.decay <- nls(abs(EC50_logFC) ~ a * exp(-b * ABI1_dist),
                                         data = data.all.posFC,
                                         start = list(a = 4, b = 0.5))

PYL1.negEC50_logFC.ABI1.exp.decay <- nls(abs(EC50_logFC) ~ a * exp(-b * ABI1_dist),
                                         data = data.all.negFC,
                                         start = list(a = 4, b = 0.5))

## calculate exponential distance decay regression confidence intervals
PYL1.posEC50_logFC.ABI1.exp.decay.full <- data.frame(ABI1_dist = seq(0, 36, length.out = 1000))
PYL1.negEC50_logFC.ABI1.exp.decay.full <- data.frame(ABI1_dist = seq(0, 36, length.out = 1000))

#### predict fitted values
PYL1.posEC50_logFC.ABI1.exp.decay.full$fit <- predict(PYL1.posEC50_logFC.ABI1.exp.decay, newdata = PYL1.posEC50_logFC.ABI1.exp.decay.full)
PYL1.negEC50_logFC.ABI1.exp.decay.full$fit <- predict(PYL1.negEC50_logFC.ABI1.exp.decay, newdata = PYL1.negEC50_logFC.ABI1.exp.decay.full)

#### bootstrap the model parameters
set.seed(123)
PYL1.posEC50_logFC.ABI1.boot_params <- replicate(1000, {
  samp <- data.all.posFC[sample(nrow(data.all.posFC), replace = TRUE), ]
  fit <- try(nls(abs(EC50_logFC) ~ a * exp(-b * ABI1_dist), data = samp, start = list(a = 4, b = 0.5)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    c(NA, NA)
  } else {
    coef(fit)
  }
})
PYL1.posEC50_logFC.ABI1.boot_params <- t(PYL1.posEC50_logFC.ABI1.boot_params)
PYL1.posEC50_logFC.ABI1.boot_params <- PYL1.posEC50_logFC.ABI1.boot_params[complete.cases(PYL1.posEC50_logFC.ABI1.boot_params), ]
PYL1.posEC50_logFC.ABI1.preds <- apply(PYL1.posEC50_logFC.ABI1.boot_params, 1, function(p) {
  p[1] * exp(-p[2] * PYL1.posEC50_logFC.ABI1.exp.decay.full$ABI1_dist)
})
PYL1.posEC50_logFC.ABI1.exp.decay.full$lower <- apply(PYL1.posEC50_logFC.ABI1.preds, 1, quantile, probs = 0.025)
PYL1.posEC50_logFC.ABI1.exp.decay.full$upper <- apply(PYL1.posEC50_logFC.ABI1.preds, 1, quantile, probs = 0.975)

set.seed(123)
PYL1.negEC50_logFC.ABI1.boot_params <- replicate(1000, {
  samp <- data.all.negFC[sample(nrow(data.all.negFC), replace = TRUE), ]
  fit <- try(nls(abs(EC50_logFC) ~ a * exp(-b * ABI1_dist), data = samp, start = list(a = 4, b = 0.5)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    c(NA, NA)
  } else {
    coef(fit)
  }
})
PYL1.negEC50_logFC.ABI1.boot_params <- t(PYL1.negEC50_logFC.ABI1.boot_params)
PYL1.negEC50_logFC.ABI1.boot_params <- PYL1.negEC50_logFC.ABI1.boot_params[complete.cases(PYL1.negEC50_logFC.ABI1.boot_params), ]
PYL1.negEC50_logFC.ABI1.preds <- apply(PYL1.negEC50_logFC.ABI1.boot_params, 1, function(p) {
  p[1] * exp(-p[2] * PYL1.negEC50_logFC.ABI1.exp.decay.full$ABI1_dist)
})
PYL1.negEC50_logFC.ABI1.exp.decay.full$lower <- apply(PYL1.negEC50_logFC.ABI1.preds, 1, quantile, probs = 0.025)
PYL1.negEC50_logFC.ABI1.exp.decay.full$upper <- apply(PYL1.negEC50_logFC.ABI1.preds, 1, quantile, probs = 0.975)

#### display main parameters of the fit:
PYL1.posEC50_logFC.ABI1.Gb0 <- 10^c(coef(PYL1.posEC50_logFC.ABI1.exp.decay)["a"],
                                    coef(PYL1.posEC50_logFC.ABI1.exp.decay)["a"] - 1.96*summary(PYL1.posEC50_logFC.ABI1.exp.decay)[[10]][1,2],
                                    coef(PYL1.posEC50_logFC.ABI1.exp.decay)["a"] + 1.96*summary(PYL1.posEC50_logFC.ABI1.exp.decay)[[10]][1,2])
PYL1.posEC50_logFC.ABI1.rate <- c(coef(PYL1.posEC50_logFC.ABI1.exp.decay)["b"],
                                  coef(PYL1.posEC50_logFC.ABI1.exp.decay)["b"] - 1.96*summary(PYL1.posEC50_logFC.ABI1.exp.decay)[[10]][2,2],
                                  coef(PYL1.posEC50_logFC.ABI1.exp.decay)["b"] + 1.96*summary(PYL1.posEC50_logFC.ABI1.exp.decay)[[10]][2,2])
PYL1.posEC50_logFC.ABI1.Gb_half <- c(log(2)/coef(PYL1.posEC50_logFC.ABI1.exp.decay)["b"],
                                     quantile(log(2)/PYL1.posEC50_logFC.ABI1.boot_params[,"b"], probs = c(0.025, 0.975)))

PYL1.negEC50_logFC.ABI1.Gb0 <- 10^-c(coef(PYL1.negEC50_logFC.ABI1.exp.decay)["a"],
                                     coef(PYL1.negEC50_logFC.ABI1.exp.decay)["a"] - 1.96*summary(PYL1.negEC50_logFC.ABI1.exp.decay)[[10]][1,2],
                                     coef(PYL1.negEC50_logFC.ABI1.exp.decay)["a"] + 1.96*summary(PYL1.negEC50_logFC.ABI1.exp.decay)[[10]][1,2])
PYL1.negEC50_logFC.ABI1.rate <- c(coef(PYL1.negEC50_logFC.ABI1.exp.decay)["b"],
                                  coef(PYL1.negEC50_logFC.ABI1.exp.decay)["b"] - 1.96*summary(PYL1.negEC50_logFC.ABI1.exp.decay)[[10]][2,2],
                                  coef(PYL1.negEC50_logFC.ABI1.exp.decay)["b"] + 1.96*summary(PYL1.negEC50_logFC.ABI1.exp.decay)[[10]][2,2])
PYL1.negEC50_logFC.ABI1.Gb_half <- c(log(2)/coef(PYL1.negEC50_logFC.ABI1.exp.decay)["b"],
                                     quantile(log(2)/PYL1.negEC50_logFC.ABI1.boot_params[,"b"], probs = c(0.025, 0.975)))


## 4. Plot ##
#############

out.1 <- ggplot(data.all.posFC, aes(x = ABI1_dist, y = abs(EC50_logFC))) +
  geom_point(data = data.all.posFC, 
             aes(x = ABI1_dist, y = abs(EC50_logFC), color = squish(abs(EC50_logFC), range = c(-4.5, 4.5))), 
             shape = 16, size = 1.5, alpha = 1) +  
  scale_colour_gradient2(low = "blue", mid = "gray90", high = "red", 
                         midpoint = 0, 
                         limits = c(-2,2),
                         oob = squish) +
  geom_ribbon(data = PYL1.posEC50_logFC.ABI1.exp.decay.full, aes(x = ABI1_dist, y = fit, ymin = lower, ymax = upper), alpha = 0.2, fill = "grey") + 
  geom_line(data = PYL1.posEC50_logFC.ABI1.exp.decay.full, aes(x = ABI1_dist, y = fit), color = "black", linewidth = 1.5) +
  geom_vline(xintercept = 5, color = "black", linewidth = 0.75) +
  scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5), limits = c(0, 37)) +
  scale_y_continuous(breaks = seq(from = 0, to = 4, by = 1), limits = c(-10000,10000),
                     labels = c("1", "10", "100", "1,000", "10,000")) +
  coord_cartesian(xlim = c(0, 37), ylim = c(0, 4), expand = F) +
  annotate("text",
           x = 16.5,
           y = 3.8,
           label = bquote(log[10]("k") == .(format(-round(PYL1.posEC50_logFC.ABI1.rate[1], 2), nsmall = 2)) * " Å"^{-1}),
           hjust = 0, size = 10, color = "black") +
  annotate("text",
           x = 16.5,
           y = 3.5,
           label = bquote(log[10](d[50]) == .(format(round(PYL1.posEC50_logFC.ABI1.Gb_half[1], 2), nsmall = 2)) * " Å"),
           hjust = 0, size = 10, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 30),
        legend.position = "none",
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = "Distance to ABI1 [Å]",
       y = bquote(EC[50] * " fold change"))

out.2 <- ggplot(data.all.negFC, aes(x = ABI1_dist, y = abs(EC50_logFC))) +
  geom_point(data = data.all.negFC, 
             aes(x = ABI1_dist, y = abs(EC50_logFC), color = squish(abs(EC50_logFC), range = c(-4.5, 4.5))), 
             shape = 16, size = 1.5, alpha = 1) +  
  scale_colour_gradient2(low = "red", mid = "gray90", high = "blue", 
                         midpoint = 0, 
                         limits = c(-2,2),
                         oob = squish) +
  geom_point(data = data.all.negFC[grep("87|164|165|166", rownames(data.all.negFC)),], 
             aes(x = ABI1_dist, y = abs(EC50_logFC), color = squish(abs(EC50_logFC), range = c(-4.5, 4.5))), 
             shape = 16, color = "black", 
             size = 4, alpha = 1) +  
  geom_ribbon(data = PYL1.negEC50_logFC.ABI1.exp.decay.full, aes(x = ABI1_dist, y = fit, ymin = lower, ymax = upper), alpha = 0.2, fill = "grey") + 
  geom_line(data = PYL1.negEC50_logFC.ABI1.exp.decay.full, aes(x = ABI1_dist, y = fit), color = "black", linewidth = 1.5) +
  geom_vline(xintercept = 5, color = "black", linewidth = 0.75) +
  scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5), limits = c(0, 37)) +
  scale_y_continuous(breaks = seq(from = 0, to = 2, by = 1), limits = c(-10000,10000),
                     labels = c("1", "0.1", "0.01")) +
  coord_cartesian(xlim = c(0, 37), ylim = c(0, 2), expand = F) +
  annotate("text",
           x = 16.5,
           y = 1.9,
           label = bquote(log[10]("k") == .(format(-round(PYL1.negEC50_logFC.ABI1.rate[1], 2), nsmall = 2)) * " Å"^{-1}),
           hjust = 0, size = 10, color = "black") +
  annotate("text",
           x = 16.5,
           y = 1.75,
           label = bquote(log[10](d[50]) == .(format(round(PYL1.negEC50_logFC.ABI1.Gb_half[1], 2), nsmall = 2)) * " Å"),
           hjust = 0, size = 10, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 30),
        legend.position = "none",
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = "Distance to ABI1 [Å]",
       y = bquote(EC[50] * " fold change"))

pdf("../../results/FigureS3/FigureS3D_PYL1_EC50_FC_decay_ABI1.pdf", width = 20, height = 9)

plot_grid(out.1, out.2, align = "hv", axis = "tblr", ncol = 2)

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
# [1] cowplot_1.1.3 ggtext_0.1.2  ggplot2_3.5.1 reshape_0.8.9 scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.4     generics_0.1.3    labeling_0.4.3   
# [8] glue_1.8.0        colorspace_2.1-1  plyr_1.8.9        gridtext_0.1.5    grid_4.4.1        munsell_0.5.1     tibble_3.2.1     
# [15] lifecycle_1.0.4   compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.14       pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2     
# [22] R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    withr_3.0.2       tools_4.4.1       gtable_0.3.6     
# [29] xml2_1.3.6    