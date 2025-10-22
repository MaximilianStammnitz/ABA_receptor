# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 22.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#################################################
## Figure 1F - external Ala variant validation ##
#################################################


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
PYL1.ABI1.mat <- PYL1.ABI1.sigma.mat <- matrix(NA, ncol = 12, nrow = 13)
colnames(PYL1.ABI1.mat) <- colnames(PYL1.ABI1.sigma.mat) <- names(PYL1.ABI1)
rownames(PYL1.ABI1.mat) <- rownames(PYL1.ABI1.sigma.mat) <- PYL1.ABI1[[1]]$id
for(i in 1:12){
  PYL1.ABI1.mat[,i] <- PYL1.ABI1[[i]][match(rownames(PYL1.ABI1.mat), PYL1.ABI1[[i]]$id),"gr_normalised"]
  PYL1.ABI1.sigma.mat[,i] <- PYL1.ABI1[[i]][match(rownames(PYL1.ABI1.mat), PYL1.ABI1[[i]]$id),"gr_sigma_normalised"]
}

### Calculate dose-response curves

### Create a list of dose-response curves
PYL1.ABI1.ls <- vector(mode = "list", length = nrow(PYL1.ABI1.mat))
names(PYL1.ABI1.ls) <- rownames(PYL1.ABI1.mat)

### Fit curves using the DRC package, predict WT-relative binding (+ SE) at 10 µM
for (i in 1:length(PYL1.ABI1.ls)){

  ## Temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.mat[i,],PYL1.ABI1.sigma.mat[i,],colnames(PYL1.ABI1.mat))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("GR", "GR sigma", "concentration")

  ## Curve fitting & parameter extraction
  tmp.drc <- drm(tmp.drc$GR ~ tmp.drc$concentration,
                 weights = 1/tmp.drc$`GR sigma`,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]

  ## Predict binding at 10 µM
  tmp.drc.predict.newdata <- expand.grid(conc = 10)
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence")
  PYL1.ABI1.ls[[i]] <- tmp.drc.predict

  ## Clean up environment
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict, tmp.drc.predict.newdata)

}

### Stratify
PYL1.ABI1.10uM <- do.call(rbind, PYL1.ABI1.ls)
rownames(PYL1.ABI1.10uM)[nrow(PYL1.ABI1.10uM)] <- "WT"

### Normalise other values to the WT at 10 µM (not at WT B[inf]!)
PYL1.ABI1.10uM <- 100*PYL1.ABI1.10uM[1:12,]/PYL1.ABI1.10uM[13,1]


## 2. Miyazono et al. Nature 2009 alanine mutants ##
####################################################

## input ABI1 pull-down % values from Miyazono et al., Nature 2009
miyazono <- c(0.94, 0.52, 0.90, 0.46, 0.52,
              0.74, 0.13, 0.30, 0.06, 0.89,
              0.66, 0.30)
names(miyazono) <- c("H87A", "F88A", "I111A", "S112A", "L114A",
                     "P115A", "H142A", "R143A", "L144A", "P178A",
                     "N181A", "F189A")


## 3. Combine two the data sets ##
##################################

out.miya.10uM.df <- cbind("miyazono" = 100*miyazono, "10 µM mean" = as.numeric(PYL1.ABI1.10uM[,"Prediction"]), 
                          "10 µM lower" = as.numeric(PYL1.ABI1.10uM[,"Lower"]), "10 µM upper" = as.numeric(PYL1.ABI1.10uM[,"Upper"]))
out.miya.10uM.df <- as.data.frame(out.miya.10uM.df)

## correlation
r <- cor(x = out.miya.10uM.df$miyazono,
         y = out.miya.10uM.df$`10 µM mean`,
         method = "pearson")


## ggPlots
pdf("../../results/Figure1/Figure1F_Miyazono_correlation.pdf", height = 15, width = 18)
 
out.1F <- ggplot(out.miya.10uM.df, aes(x = `miyazono`, y = `10 µM mean`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  coord_cartesian(xlim = c(-5, 121), ylim = c(-5, 115)) +
  geom_smooth(out.miya.10uM.df,
              mapping = aes(x = `miyazono`, y = `10 µM mean`),
              method = 'lm',
              color = "black",
              fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = "grey50") +
  geom_errorbar(aes(ymin = `10 µM lower`, ymax = `10 µM upper`),
                 color = "black", linewidth = 1) +
  geom_point(data = out.miya.10uM.df,
             mapping = aes(x = `miyazono`, y = `10 µM mean`),
             color = "black", size = 8, shape = 21, fill = "#569956") +
  geom_text_repel(data = out.miya.10uM.df,
                  aes(x = `miyazono`, y = `10 µM mean`,
                      label = rownames(out.miya.10uM.df)),
                  max.overlaps = 20,
                  force = 3,
                  min.segment.length = 0,
                  box.padding = 1,
                  point.padding = 1,
                  color = "black",
                  size = 15) +
  annotate("text",
           x = 0,
           y = 105,
           label = expr_text(bquote(italic(r) == .(format(r, digits = 2)))),
           parse = T,
           hjust = 0, size = 25, color = "black") +
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
  labs(x = "Relative PYL1/ABI1 Binding (Miyazono et al., 2009)",
       y = "Relative PYL1/ABI1 Binding (this study)")

print(out.1F)

dev.off()


## 4. Version ##
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
