# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

############################################################################
## Supplementary Figure 2D - Potency shift of (+)-ABA GluePCA vs in vitro ##
############################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("drc", "ggrepel", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## generate mutant's dose response curves
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[-which(x[,"Nham_aa"] == 0)[-1],]; return(x)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x$id <- paste0(x[,"WT_AA"], x[,"Pos"], x[,"Mut"]); return(x)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[grep("WT|H87A|F88A|I111A|S112A|L114A|P115A|H142A|R143A|L144A|P178A|N181A|F189A", x$id),]; return(x)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[order(as.numeric(x$Pos)),]; return(x)})

### Build a dose-response matrix for each nucleotide sequence
PYL1.ABI1.mat <- matrix(NA, ncol = 12, nrow = 13)
colnames(PYL1.ABI1.mat) <- names(PYL1.ABI1)
rownames(PYL1.ABI1.mat) <- PYL1.ABI1[[1]]$id
for(i in 1:12){
  PYL1.ABI1.mat[,i] <- PYL1.ABI1[[i]][match(rownames(PYL1.ABI1.mat), PYL1.ABI1[[i]]$id),"gr_normalised"]
}

### Calculate dose-response curves

### Create a list of dose-response curves
PYL1.ABI1.ls <- vector(mode = "list", length = nrow(PYL1.ABI1.mat))
names(PYL1.ABI1.ls) <- rownames(PYL1.ABI1.mat)

### Fit curves using the DRC package
for (i in 1:length(PYL1.ABI1.ls)){
  
  ## Temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.mat[i,],colnames(PYL1.ABI1.mat))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("GR", "concentration")
  
  ## Curve fitting & parameter extraction
  tmp.drc <- drm(tmp.drc$GR ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]
  
  ## Predict binding across 1000 concentrations
  tmp.drc.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence")
  PYL1.ABI1.ls[[i]] <- tmp.drc.predict[,1]
  
  ## Clean up environment
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict, tmp.drc.predict.newdata)
  
}

### Stratify & set relative to WT
PYL1.ABI1.ls <- do.call(rbind, PYL1.ABI1.ls)
colnames(PYL1.ABI1.ls) <- c(0, exp(seq(log(0.001), log(5000), length = 999)))
PYL1.ABI1.ls <- apply(PYL1.ABI1.ls, 2, function(x){y <- 100*x[1:12]/x[13]; return(y)})


## 2. Miyazono et al. Nature 2009 alanine mutants ##
####################################################

## input ABI1 pull-down % values from Miyazono et al., Nature 2009
miyazono <- 100*c(0.94, 0.52, 0.90, 0.46, 0.52,
                  0.74, 0.13, 0.30, 0.06, 0.89,
                  0.66, 0.30)
names(miyazono) <- c("H87A", "F88A", "I111A", "S112A", "L114A",
                     "P115A", "H142A", "R143A", "L144A", "P178A",
                     "N181A", "F189A")


## 3. Calculate correlations across 1,000 concentrations ##
###########################################################

## correlate
out.cor <- cbind("conc" = c(9.062741e-03/3.5/3.5/3.5, exp(seq(log(0.001), log(5000), length = 999))), 
                 "corr" = rep(NA,1000))
for(i in 1:1000){
  out.cor[i,2] <- cor(x = miyazono,
                      y = PYL1.ABI1.ls[,i],
                      method = "pearson")
}
out.cor <- as.data.frame(out.cor)

## ggPlots
pdf("../../results/FigureS2/FigureS2D_potency_shift.pdf", height = 15, width = 1.130682*18)

out.S2D <- ggplot(out.cor[-1,], aes(x = `conc`, y = `corr`)) +
  geom_line(data = out.cor[-1,],
             mapping = aes(x = `conc`, y = `corr`),
             color = "black", size = 3) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, length.out = 6), limits = c(0, 1.1)) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_segment(aes(x = 10,
                   xend = 10,
                   y = out.cor[which.min(abs(out.cor$conc - 10)),2] + 0.015, 
                   yend = out.cor[which.min(abs(out.cor$conc - 10)),2] + 0.035),
               size = 1, color = "black") +
  annotate("text", 
           x = 10, 
           y = out.cor[which.min(abs(out.cor$conc - 10)),2] + 0.1, 
           label = expr_text(bquote(italic(r[obs.]))), 
           parse = T,
           vjust = 1, 
           color = "black",
           size = 20) +
  geom_segment(aes(x = out.cor[which.max(out.cor[,2]),1],
                   xend = out.cor[which.max(out.cor[,2]),1],
                   y = out.cor[which.max(out.cor[,2]),2] + 0.015, 
                   yend = out.cor[which.max(out.cor[,2]),2] + 0.035),
               size = 1, color = "black") +
  annotate("text",
           x = out.cor[which.max(out.cor[,2]),1],
           y = out.cor[which.max(out.cor[,2]),2] + 0.1,
           label = expr_text(bquote(italic(r[max.]))),
           parse = T,
           vjust = 1,
           color = "black",
           size = 20) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "(+)-ABA conc. (µM)",
       y = expression("Pearson's r (GluePCA vs. 10 µM" ~ italic("in vitro") ~ ")"))

print(out.S2D)

dev.off()


## 4. Small-scale plots for robs and rmax ##
############################################

### Re-calculate two dose-response curves: for 10 µM (observed) and for 136.45266 µM (theoretical max. r)
PYL1.ABI1.two <- vector(mode = "list", length = nrow(PYL1.ABI1.mat))
names(PYL1.ABI1.two) <- rownames(PYL1.ABI1.mat)

### Fit curves using the DRC package & predict WT-relative binding
for (i in 1:length((PYL1.ABI1.two))){
  
  ## Temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.mat[i,],colnames(PYL1.ABI1.mat))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("GR", "concentration")
  
  ## Curve fitting & parameter extraction
  tmp.drc <- drm(tmp.drc$GR ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]
  
  ## Predict binding at 10 µM
  tmp.drc.predict.newdata <- expand.grid(conc = c(10, 136.45266))
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence")
  PYL1.ABI1.two[[i]] <- tmp.drc.predict[,1]
  
  ## Clean up environment
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict, tmp.drc.predict.newdata)
  
}

### Stratify & set relative to WT
PYL1.ABI1.two <- do.call(rbind, PYL1.ABI1.two)
colnames(PYL1.ABI1.two) <- c(10, 136.45266)
PYL1.ABI1.two <- apply(PYL1.ABI1.two, 2, function(x){y <- 100*x[1:12]/x[13]; return(y)})

out.miya.10uM.df <- as.data.frame(cbind("miyazono" = miyazono, "10 µM" = PYL1.ABI1.two[,"10"]))
out.miya.136uM.df <- as.data.frame(cbind("miyazono" = miyazono, "136 µM" = PYL1.ABI1.two[,"136.45266"]))
r.10 <- cor(x = out.miya.10uM.df$miyazono, y = out.miya.10uM.df$`10 µM`, method = "pearson")
r.136 <- cor(x = out.miya.136uM.df$miyazono, y = out.miya.136uM.df$`136 µM`, method = "pearson")

out.S2D_obs <- ggplot(out.miya.10uM.df, aes(x = `miyazono`, y = `10 µM`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  coord_cartesian(xlim = c(-5, 121), ylim = c(-5, 115)) +
  geom_point(data = out.miya.10uM.df,
             mapping = aes(x = `miyazono`, y = `10 µM`),
             color = "black", size = 16, shape = 16) +
  geom_smooth(out.miya.10uM.df,
              mapping = aes(x = `miyazono`, y = `10 µM`),
              method = 'lm',
              color = "black",
              fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = "grey50") +
  annotate("text",
           x = -5,
           y = 105,
           label = expr_text(bquote(italic(r) == .(format(r.10, digits = 2)))),
           parse = T,
           hjust = 0, size = 35, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1/ABI1 Binding (Miyazono et al., 2009)",
       y = "Relative PYL1/ABI1 Binding (this study, 10 µM)")

pdf("../../results/FigureS2/FigureS2D_potency_shift_obs_v2.pdf", height = 15, width = 18)
print(out.S2D_obs)
dev.off()

out.S2D_theormax <- ggplot(out.miya.136uM.df, aes(x = `miyazono`, y = `136 µM`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-100, 200)) +
  scale_y_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-100, 200)) +
  coord_cartesian(xlim = c(-5, 121), ylim = c(-5, 115)) +
  geom_point(data = out.miya.136uM.df,
             mapping = aes(x = `miyazono`, y = `136 µM`),
             color = "black", size = 16, shape = 16) +
  geom_smooth(out.miya.136uM.df,
              mapping = aes(x = `miyazono`, y = `136 µM`),
              method = 'lm',
              color = "black",
              fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = "grey50") +
  annotate("text",
           x = -5,
           y = 105,
           label = expr_text(bquote(italic(r) == .(format(r.136, digits = 2)))),
           parse = T,
           hjust = 0, size = 35, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1/ABI1 Binding (Miyazono et al., 2009)",
       y = "Relative PYL1/ABI1 Binding (this study, 136.5 µM)")

pdf("../../results/FigureS2/FigureS2D_potency_shift_theor_max_v2.pdf", height = 15, width = 18)
print(out.S2D_theormax)
dev.off()


## 5. Version ##
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
# [1] rlang_1.1.6   ggtext_0.1.2  ggrepel_0.9.6 ggplot2_4.0.0 drc_3.0-1     MASS_7.3-65  
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       crayon_1.5.3       dplyr_1.1.4        compiler_4.5.1     gtools_3.9.5      
# [7] tidyselect_1.2.1   plotrix_3.8-4      Rcpp_1.1.0         xml2_1.4.0         splines_4.5.1      scales_1.4.0      
# [13] lattice_0.22-7     TH.data_1.1-4      R6_2.6.1           generics_0.1.4     Formula_1.2-5      tibble_3.3.0      
# [19] car_3.1-3          pillar_1.11.0      RColorBrewer_1.1-3 multcomp_1.4-28    S7_0.2.0           cli_3.6.5         
# [25] mgcv_1.9-3         withr_3.0.2        magrittr_2.0.3     gridtext_0.1.5     grid_4.5.1         rstudioapi_0.17.1 
# [31] mvtnorm_1.3-3      sandwich_3.1-1     nlme_3.1-168       lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0        
# [37] farver_2.1.2       codetools_0.2-20   zoo_1.8-14         survival_3.8-3     abind_1.4-8        carData_3.0-5     
# [43] pkgconfig_2.0.3    tools_4.5.1  
