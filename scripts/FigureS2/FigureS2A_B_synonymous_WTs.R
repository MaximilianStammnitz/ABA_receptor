# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################################
## Supplementary Figure 2A & 2B - PYL1 synonymous WTs ##
########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "drc", "reshape", "ggplot2", "ggtext", "scales", "beeswarm", "rlang")

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
PYL1.ABI1.WT.mat <- matrix(NA, ncol = 12, nrow = length(unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))))
colnames(PYL1.ABI1.WT.mat) <- names(PYL1.ABI1.WT)
rownames(PYL1.ABI1.WT.mat) <- unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))
for(i in 1:12){
  PYL1.ABI1.WT.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_normalised_WTscaled"]
}

## Remove synonymous variants not fully covered across (+)-ABA concentrations
PYL1.ABI1.WT.mat <- na.omit(PYL1.ABI1.WT.mat)

## Clean up environment
rm(packages, install_if_missing, i)


## 2. Calculate dose-response curves ##
#######################################

## Create a list of dose-response curves
PYL1.ABI1.WT.ls <- PYL1.ABI1.WT.params <- vector(mode = "list", length = nrow(PYL1.ABI1.WT.mat))

## Fit curves using the DRC packages
for (i in 1:length(PYL1.ABI1.WT.ls)){
  
  ## Temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.WT.mat[i,],colnames(PYL1.ABI1.WT.mat))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")
  
  ## Curve fitting & parameter extraction
  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]
  
  ## Predict the full, smoothened curve (using 1000 data points) and confidence interval
  tmp.drc.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence")
  tmp.drc.predict.newdata$p <- tmp.drc.predict[,1]
  tmp.drc.predict.newdata$pmin <- tmp.drc.predict[,2]
  tmp.drc.predict.newdata$pmax <- tmp.drc.predict[,3]
  
  ## Save
  PYL1.ABI1.WT.params[[i]] <- tmp.drc.par
  PYL1.ABI1.WT.ls[[i]] <- tmp.drc.predict.newdata
  
  ## Clean up environment
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict, tmp.drc.predict.newdata)
  
}

## Summarise all curve fits in a single data frame
PYL1.ABI1.WT.curves <- cbind(sapply(PYL1.ABI1.WT.ls, function(x) x$p), 
                             "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
PYL1.ABI1.WT.curves <- as.data.frame(PYL1.ABI1.WT.curves)
PYL1.ABI1.WT.curves <- melt(PYL1.ABI1.WT.curves, id.vars = "conc")


## 3. Visualise DRCs ##
#######################

## Visualise all WT dose-response curves
pdf("../../results/FigureS2/FigureS2A_synonymous_WT_dose_response_curves.pdf", 
    height = 15, width = 18)

out.S2A <- ggplot(data = PYL1.ABI1.WT.curves) +
  geom_line(data = PYL1.ABI1.WT.curves, aes(group = variable, x = conc, y = value), linewidth = 1, alpha = 0.05) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
  coord_cartesian(ylim = c(-5, 115)) +
  annotate("text",
           x = 9.062741e-03/3.5/3.5/3.5,
           y = 105,
           label = expr_text(bquote(italic(N) == .(format(length(PYL1.ABI1.WT.ls), digits = 2)))),
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
print(out.S2A)

dev.off()


## 4. Visualise parameters ##
#############################

## Fetch chosen "main" WT parameters
PYL1.ABI1.WT.params.main <- matrix(NA, ncol = 3, nrow = 4)
colnames(PYL1.ABI1.WT.params.main) <- c("Mean", "0.05", "0.95")
rownames(PYL1.ABI1.WT.params.main) <- c("B[0]", "B[inf]", "EC50", "Hill")
tmp.drc <- cbind(PYL1.ABI1.WT.mat[1,],colnames(PYL1.ABI1.WT.mat))
class(tmp.drc) <- "numeric"
tmp.drc <- as.data.frame(tmp.drc)
colnames(tmp.drc) <- c("B", "concentration")
tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
               fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
               type = 'continuous')
tmp.drc.summary <- summary(tmp.drc)
PYL1.ABI1.WT.params.main["B[0]",] <- c(tmp.drc$parmMat[2], tmp.drc$parmMat[2] - tmp.drc.summary$coefficients[2,2]*1.96, tmp.drc$parmMat[2] + tmp.drc.summary$coefficients[2,2]*1.96)
PYL1.ABI1.WT.params.main["B[inf]",] <- c(tmp.drc$parmMat[3], tmp.drc$parmMat[3] - tmp.drc.summary$coefficients[3,2]*1.96, tmp.drc$parmMat[3] + tmp.drc.summary$coefficients[3,2]*1.96)
PYL1.ABI1.WT.params.main["EC50",] <- exp(c(log(tmp.drc$parmMat[4]), 
                                           log(tmp.drc$parmMat[4]) - c(c(tmp.drc.summary$coefficients[4,2]/tmp.drc$parmMat[4])*1.96), 
                                           log(tmp.drc$parmMat[4]) + c(c(tmp.drc.summary$coefficients[4,2]/tmp.drc$parmMat[4])*1.96)))
PYL1.ABI1.WT.params.main["Hill",] <- exp(c(log(-tmp.drc$parmMat[1]), 
                                           log(-tmp.drc$parmMat[1]) - c(c(tmp.drc.summary$coefficients[1,2]/c(-tmp.drc$parmMat[1]))*1.96), 
                                           log(-tmp.drc$parmMat[1]) + c(c(tmp.drc.summary$coefficients[1,2]/c(-tmp.drc$parmMat[1]))*1.96)))

## Look at parameter output stats for all other curves
PYL1.ABI1.WT.params <- do.call(cbind, PYL1.ABI1.WT.params[-1])
### quantile(PYL1.ABI1.WT.params["B[0]",-1], c(0.05, 0.95))
### quantile(PYL1.ABI1.WT.params["B[inf]",-1], c(0.05, 0.95))
### exp(quantile(log(PYL1.ABI1.WT.params["EC50",-1]), c(0.05, 0.95)))
### exp(quantile(log(PYL1.ABI1.WT.params["Hill",-1]), c(0.05, 0.95)))


## plot

### B[0]
pdf('../../results/FigureS2/FigureS2B_B0.pdf', 
    width = 6, height = 12)

par(mar = c(22, 18, 1, 1))

## boxplot frame
boxplot(list("WT" = c(),
             "synonymous WT" = PYL1.ABI1.WT.params["B[0]",]),
        ylim = c(-5,112),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 3,
        las = 2,
        cex.main = 2)

## chosen WT 95% conf. interval
## error bars
arrows(x0 = 1,
       x1 = 1,
       y1 = PYL1.ABI1.WT.params.main["B[0]","0.95"], 
       y0 = PYL1.ABI1.WT.params.main["B[0]","0.05"],
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'grey50')
points(x = 1, y = PYL1.ABI1.WT.params.main["B[0]","Mean"], pch = 16, cex = 2)


## synonymous WT data (as beeswarm)
beeswarm(list(c(),
              PYL1.ABI1.WT.params["B[0]",]), 
         method = 'compactswarm',
         pch = 16,
         cex = 0.45,
         add = T,
         col = "black")

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 100, length.out = 6), 
     labels = seq(from = 0, to = 100, length.out = 6),
     las = 2,
     line = 1)
mtext(text = bquote(B[0] ~ "(Binding)"), side = 2, line = 10, cex = 4)

dev.off()

### B[inf]
pdf('../../results/FigureS2/FigureS2B_Binf.pdf', 
    width = 6, height = 12)

par(mar = c(22, 18, 1, 1))

## boxplot frame
boxplot(list("WT" = c(),
             "synonymous WT" = PYL1.ABI1.WT.params["B[inf]",]),
        ylim = c(-5,112),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 3,
        las = 2,
        cex.main = 2)

## chosen WT 95% conf. interval
## error bars
arrows(x0 = 1,
       x1 = 1,
       y1 = PYL1.ABI1.WT.params.main["B[inf]","0.95"], 
       y0 = PYL1.ABI1.WT.params.main["B[inf]","0.05"],
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'grey50')
points(x = 1, y = PYL1.ABI1.WT.params.main["B[inf]","Mean"], pch = 16, cex = 2)

## synonymous WT data (as beeswarm)
beeswarm(list(c(),
              PYL1.ABI1.WT.params["B[inf]",]), 
         method = 'compactswarm',
         pch = 16,
         cex = 0.45,
         add = T,
         col = "black")

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 100, length.out = 6), 
     labels = seq(from = 0, to = 100, length.out = 6),
     las = 2,
     line = 1)
mtext(text = bquote(B[infinity] ~ "(Binding)"), side = 2, line = 10, cex = 4)

dev.off()

## EC50
pdf('../../results/FigureS2/FigureS2B_EC50.pdf', 
    width = 6, height = 12)

par(mar = c(22, 18, 1, 1))

## boxplot frame
boxplot(list("WT" = c(),
             "synonymous WT" = PYL1.ABI1.WT.params["EC50",]),
        ylim = c(0.01,1000),
        log = "y",
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 3,
        las = 2,
        cex.main = 2)

## chosen WT 95% conf. interval
## error bars
arrows(x0 = 1,
       x1 = 1,
       y1 = PYL1.ABI1.WT.params.main["EC50","0.95"], 
       y0 = PYL1.ABI1.WT.params.main["EC50","0.05"],
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'grey50')
points(x = 1, y = PYL1.ABI1.WT.params.main["EC50","Mean"], pch = 16, cex = 2)

## synonymous WT data (as beeswarm)
beeswarm(list(c(),
              PYL1.ABI1.WT.params["EC50",]), 
         method = 'compactswarm',
         pch = 16,
         cex = 0.45,
         add = T,
         col = "black")

## Axis
axis(2, 
     cex.axis = 3,
     at = c(0.01, 0.1, 1, 10, 100, 1000), 
     labels = c(0.01, 0.1, 1, 10, 100, c("1,000")),
     las = 2,
     line = 1)

mtext(text = bquote(EC[50] ~ "((+)-ABA, µM)"), side = 2, line = 10, cex = 4)

dev.off()

## Hill
pdf('../../results/FigureS2/FigureS2B_Hill.pdf', 
    width = 6, height = 12)

par(mar = c(22, 18, 1, 1))

## boxplot frame
boxplot(list("WT" = c(),
             "synonymous WT" = PYL1.ABI1.WT.params["Hill",]),
        ylim = c(0,2.82),
        yaxt = "n",
        main = '',
        col = 'white',
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        staplecol = 'white',
        frame = F,
        ylab = '',
        cex.lab = 3,
        cex.axis = 3,
        las = 2,
        cex.main = 2)

## chosen WT 95% conf. interval
## error bars
arrows(x0 = 1,
       x1 = 1,
       y1 = PYL1.ABI1.WT.params.main["Hill","0.95"], 
       y0 = PYL1.ABI1.WT.params.main["Hill","0.05"],
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'grey50')
points(x = 1, y = PYL1.ABI1.WT.params.main["Hill","Mean"], pch = 16, cex = 2)

## synonymous WT data (as beeswarm)
beeswarm(list(c(),
              PYL1.ABI1.WT.params["Hill",]), 
         method = 'compactswarm',
         pch = 16,
         cex = 0.45,
         add = T,
         col = "black")

## Axis
axis(2, 
     cex.axis = 3,
     at = seq(from = 0, to = 2.5, by = 0.5), 
     labels = c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5"),
     las = 2,
     line = 1)
mtext(text = "Hill coefficient", side = 2, line = 10, cex = 4)

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
# [1] rlang_1.1.6    beeswarm_0.4.0 scales_1.4.0   ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 drc_3.0-1     
# [8] MASS_7.3-65    stringr_1.5.1 
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       crayon_1.5.3       dplyr_1.1.4        compiler_4.5.1     gtools_3.9.5      
# [7] tidyselect_1.2.1   plotrix_3.8-4      Rcpp_1.1.0         xml2_1.4.0         splines_4.5.1      lattice_0.22-7    
# [13] TH.data_1.1-4      R6_2.6.1           plyr_1.8.9         generics_0.1.4     Formula_1.2-5      tibble_3.3.0      
# [19] car_3.1-3          pillar_1.11.0      RColorBrewer_1.1-3 multcomp_1.4-28    stringi_1.8.7      S7_0.2.0          
# [25] cli_3.6.5          withr_3.0.2        magrittr_2.0.3     gridtext_0.1.5     grid_4.5.1         rstudioapi_0.17.1 
# [31] mvtnorm_1.3-3      sandwich_3.1-1     lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2      
# [37] codetools_0.2-20   zoo_1.8-14         survival_3.8-3     abind_1.4-8        carData_3.0-5      pkgconfig_2.0.3   
# [43] tools_4.5.1   
