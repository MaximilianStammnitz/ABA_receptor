# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 21.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#############################################
## Figure 1D - PYL1 WT dose-response curve ##
#############################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "drc", "reshape", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## Extract WT aa variants (incl. synonymous ones)
PYL1.ABI1.WT <- lapply(PYL1.ABI1, function(x){y <- x[which(x[,"Nham_aa"] == 0),]; return(y)})

## Build a dose-response matrix for each nucleotide sequence
PYL1.ABI1.WT.mat <- PYL1.ABI1.WT.sigma.mat <- matrix(NA, ncol = 12, nrow = length(unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))))
colnames(PYL1.ABI1.WT.mat) <- colnames(PYL1.ABI1.WT.sigma.mat) <- names(PYL1.ABI1.WT)
rownames(PYL1.ABI1.WT.mat) <- rownames(PYL1.ABI1.WT.sigma.mat) <- unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))
for(i in 1:12){
  PYL1.ABI1.WT.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_normalised_WTscaled"]
  PYL1.ABI1.WT.sigma.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.sigma.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_sigma_normalised_WTscaled"]
}

## Remove synonymous variants not fully covered across (+)-ABA concentrations
PYL1.ABI1.WT.mat <- na.omit(PYL1.ABI1.WT.mat)
PYL1.ABI1.WT.sigma.mat <- na.omit(PYL1.ABI1.WT.sigma.mat)

## Clean up environment
rm(packages, install_if_missing, i)


## 2. Calculate dose-response curves ##
#######################################

## Fit curves using the DRC packages
WT.PYL1.drc <- cbind(PYL1.ABI1.WT.mat[1,],PYL1.ABI1.WT.sigma.mat[1,],colnames(PYL1.ABI1.WT.mat))
class(WT.PYL1.drc) <- "numeric"
WT.PYL1.drc <- as.data.frame(WT.PYL1.drc)
colnames(WT.PYL1.drc) <- c("B", "B sigma", "concentration")

WT.PYL1.drc <- drm(WT.PYL1.drc$B ~ WT.PYL1.drc$concentration,
                   weights = 1/WT.PYL1.drc$`B sigma`,
                   fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                   type = 'continuous',
                   lowerl = c(NA, 0, NA, NA))
WT.PYL1.drc.par <- WT.PYL1.drc$fit$par
names(WT.PYL1.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.par <- WT.PYL1.drc.par[c(2:4,1)]
WT.PYL1.drc.par[4] <- -WT.PYL1.drc.par[4]

## Predict the full, smoothened curve (using 1000 data points) and confidence interval
WT.PYL1.drc$concentration[12] <- 9.062741e-03/3.5/3.5/3.5 ## "0-conc." positioning for log scale
WT.drc.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
WT.PYL1.drc.predict <- predict(WT.PYL1.drc,
                               newdata = WT.drc.predict.newdata,
                               interval = "confidence")
WT.drc.predict.newdata$p <- WT.PYL1.drc.predict[,1]
WT.drc.predict.newdata$pmin <- WT.PYL1.drc.predict[,2]
WT.drc.predict.newdata$pmax <- WT.PYL1.drc.predict[,3]

## Summarise all curve fits in a single data frame
PYL1.ABI1.WT.curves <- cbind("B" = WT.drc.predict.newdata$p, 
                             "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
PYL1.ABI1.WT.curves <- as.data.frame(PYL1.ABI1.WT.curves)


## 3. Visualise ##
##################

## Display data points only for the main (exact) WT sequence
PYL1.ABI1.WT.mat.exact <- cbind("conc" = as.numeric(colnames(PYL1.ABI1.WT.mat)),
                                "WT" = as.numeric(PYL1.ABI1.WT.mat[PYL1.ABI1.WT$`2500`[which(PYL1.ABI1.WT$`2500`$WT == T)[1],"nt_seq"],]),
                                "sigma" = as.numeric(PYL1.ABI1.WT.sigma.mat[PYL1.ABI1.WT$`2500`[which(PYL1.ABI1.WT$`2500`$WT == T)[1],"nt_seq"],]))
PYL1.ABI1.WT.mat.exact[12,1] <- 9.062741e-03/3.5/3.5/3.5 ## "0-conc." positioning for log scale
PYL1.ABI1.WT.mat.exact <- as.data.frame(PYL1.ABI1.WT.mat.exact)

## R^2 of fit
predicted_values <- predict(WT.PYL1.drc)
rss <- sum((WT.PYL1.drc$data$`WT.PYL1.drc$B` - predicted_values)^2)
tss <- sum((WT.PYL1.drc$data$`WT.PYL1.drc$B` - mean(WT.PYL1.drc$data$`WT.PYL1.drc$B`))^2)
r_squared <- 1 - (rss/tss)

## Visualise the main WT dose-response curve, including the confidence interval ribbon
pdf("../../results/Figure1/Figure1D_WT_main_dose_response_curve.pdf",
    height = 15, width = 18)

out.1D <- ggplot(data = PYL1.ABI1.WT.curves) +
  geom_ribbon(data = WT.drc.predict.newdata, aes(x = conc, y = p, ymin = pmin, ymax = pmax),
              alpha = 0.2, fill = "grey50") +
  geom_point(data = PYL1.ABI1.WT.mat.exact, aes(x = conc, y = WT),
             color = "black", size = 10, shape = 16) +
  geom_line(data = PYL1.ABI1.WT.curves, aes(x = conc, y = B), linewidth = 1.5) +
#  geom_errorbar(data = PYL1.ABI1.WT.mat.exact, 
#                aes(x = conc, 
#                    ymin = PYL1.ABI1.WT.mat.exact$WT - PYL1.ABI1.WT.mat.exact$sigma * 1.96, 
#                    ymax = PYL1.ABI1.WT.mat.exact$WT + PYL1.ABI1.WT.mat.exact$sigma * 1.96),
#                color = "black", linewidth = 0.5) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
  coord_cartesian(ylim = c(-5, 115)) +
  annotate("text",
           x = 9.062741e-03/3.5/3.5/3.5,
           y = 105,
           label = expr_text(bquote(italic(R)^2 == .(format(r_squared, digits = 2)))),
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
  labs(x = "(+)-ABA conc. (µM)",
       y = "Relative PYL1/ABI1 Binding")

print(out.1D)

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
# [1] rlang_1.1.6    ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 drc_3.0-1      MASS_7.3-65    stringr_1.5.1 
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       crayon_1.5.3       dplyr_1.1.4        compiler_4.5.1     gtools_3.9.5      
# [7] tidyselect_1.2.1   plotrix_3.8-4      Rcpp_1.1.0         xml2_1.4.0         splines_4.5.1      scales_1.4.0      
# [13] lattice_0.22-7     TH.data_1.1-4      R6_2.6.1           plyr_1.8.9         generics_0.1.4     Formula_1.2-5     
# [19] tibble_3.3.0       car_3.1-3          pillar_1.11.0      RColorBrewer_1.1-3 multcomp_1.4-28    stringi_1.8.7     
# [25] S7_0.2.0           cli_3.6.5          withr_3.0.2        magrittr_2.0.3     gridtext_0.1.5     grid_4.5.1        
# [31] rstudioapi_0.17.1  mvtnorm_1.3-3      sandwich_3.1-1     lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0        
# [37] farver_2.1.2       codetools_0.2-20   zoo_1.8-14         survival_3.8-3     abind_1.4-8        carData_3.0-5     
# [43] pkgconfig_2.0.3    tools_4.5.1
