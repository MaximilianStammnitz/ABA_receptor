# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 28.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###################################################################################
## Supplementary Figure 2C - DMS vs. microtiter plate based dose-response curves ##
###################################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "growthcurver", "reshape", "drc", 
              "scales", "ggplot2", "ggtext", "rlang",
              "grid", "gridExtra")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. TECAN data pre-processing ##
##################################

## import results
PYL1.ABI1.TECAN.plate1 <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[20:381,1:301]
rownames(PYL1.ABI1.TECAN.plate1) <- PYL1.ABI1.TECAN.plate1[,1]
colnames(PYL1.ABI1.TECAN.plate1) <- PYL1.ABI1.TECAN.plate1[1,]
PYL1.ABI1.TECAN.plate1 <- PYL1.ABI1.TECAN.plate1[,-1]
PYL1.ABI1.TECAN.plate1 <- PYL1.ABI1.TECAN.plate1[-c(1:2),]
colnames(PYL1.ABI1.TECAN.plate1) <- as.numeric(colnames(PYL1.ABI1.TECAN.plate1))/3600 ## convert to hours
class(PYL1.ABI1.TECAN.plate1) <- "numeric"
PYL1.ABI1.TECAN.plate1 <- as.data.frame(PYL1.ABI1.TECAN.plate1)

## import yeast colony matrix
setup1 <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[1:18,1:25]
rownames(setup1) <- setup1[,1]
setup1 <- setup1[,-1]
colnames(setup1) <- 1:24

## shape format
PYL1.ABI1.TECAN.plate1 <- cbind("time" = colnames(PYL1.ABI1.TECAN.plate1), 
                                t(PYL1.ABI1.TECAN.plate1))
class(PYL1.ABI1.TECAN.plate1) <- "numeric"
PYL1.ABI1.TECAN.plate1 <- as.data.frame(PYL1.ABI1.TECAN.plate1)

## run
PYL1.ABI1.TECAN.plate1.curves <- SummarizeGrowthByPlate(PYL1.ABI1.TECAN.plate1[1:262,], 
                                                        t_trim = 70,
                                                        bg_correct = "min",
                                                        plot_fit = F)

## label wells by genotype
for(i in 1:nrow(PYL1.ABI1.TECAN.plate1.curves)){
  PYL1.ABI1.TECAN.plate1.curves[i,1] <- setup1[substr(PYL1.ABI1.TECAN.plate1.curves[i,1], 1, 1),
                                               substr(PYL1.ABI1.TECAN.plate1.curves[i,1], 2, 
                                                      nchar(PYL1.ABI1.TECAN.plate1.curves[i,1]))]
}

## categorise
PYL1.ABI1.TECAN.plate1.curves.auc <- vector(mode = "list", length = length(unique(c(setup1[-c(1,17,18),]))))
names(PYL1.ABI1.TECAN.plate1.curves.auc) <- unique(c(setup1[-c(1,17,18),]))
names(PYL1.ABI1.TECAN.plate1.curves.auc) <- c(names(PYL1.ABI1.TECAN.plate1.curves.auc)[1],
                                              names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1][order(as.numeric(substr(names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1], 2, nchar(names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1]) - 1)))])
PYL1.ABI1.TECAN.plate1.curves.auc <- PYL1.ABI1.TECAN.plate1.curves.auc[-which(names(PYL1.ABI1.TECAN.plate1.curves.auc) %in% c("T118N", "S119N"))]
for(i in 1:length(PYL1.ABI1.TECAN.plate1.curves.auc)){
  PYL1.ABI1.TECAN.plate1.curves.auc[[i]] <- which(PYL1.ABI1.TECAN.plate1.curves$sample == names(PYL1.ABI1.TECAN.plate1.curves.auc)[i])
  PYL1.ABI1.TECAN.plate1.curves.auc[[i]] <- PYL1.ABI1.TECAN.plate1.curves$auc_e[PYL1.ABI1.TECAN.plate1.curves.auc[[i]]]
}
PYL1.ABI1.TECAN.plate1.curves.auc <- lapply(PYL1.ABI1.TECAN.plate1.curves.auc, rev)

### dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
dosages <- rev(dosages)

## take mean value of TECAN replicates
for(i in 1:length(PYL1.ABI1.TECAN.plate1.curves.auc)){
  PYL1.ABI1.TECAN.plate1.curves.auc[[i]] <- c(mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][1:2]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][3:4]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][5:6]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][7:8]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][9:10]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][11:12]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][13:14]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][15:16]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][17:18]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][19:20]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][21:22]),
                                              mean(PYL1.ABI1.TECAN.plate1.curves.auc[[i]][23:24]))
}

## fit plate-specific WT curves
WT.PYL1.drc.1 <- cbind(PYL1.ABI1.TECAN.plate1.curves.auc$WT, dosages)
class(WT.PYL1.drc.1) <- "numeric"
WT.PYL1.drc.1 <- as.data.frame(WT.PYL1.drc.1)
colnames(WT.PYL1.drc.1) <- c("GR", "concentration")
WT.PYL1.drc.1 <- drm(WT.PYL1.drc.1$GR ~ WT.PYL1.drc.1$concentration,
                     fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                     type = 'continuous')
WT.PYL1.drc.1.par <- WT.PYL1.drc.1$fit$par
names(WT.PYL1.drc.1.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.1.par <- WT.PYL1.drc.1.par[c(2:4,1)]
WT.PYL1.drc.1.par[4] <- -WT.PYL1.drc.1.par[4]

## linearly rescale all curves to WT Bmax (100%)
PYL1.ABI1.TECAN.plate1.curves.auc <- lapply(PYL1.ABI1.TECAN.plate1.curves.auc, function(x){y <- 100*x/WT.PYL1.drc.1.par["B[inf]"]; return(y)})

## Calculate dose-response curves using duplicate measurements
PYL1.TECAN.curves.auc.predict <- PYL1.TECAN.curves.auc.params <- vector(mode = "list", length = length(PYL1.ABI1.TECAN.plate1.curves.auc))
names(PYL1.TECAN.curves.auc.predict) <- names(PYL1.TECAN.curves.auc.params) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
for(i in 1:length(PYL1.ABI1.TECAN.plate1.curves.auc)){

  ## temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.TECAN.plate1.curves.auc[[i]],
                   dosages)
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")
  
  ## temporary Hill fit & parameters 
  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]
  PYL1.TECAN.curves.auc.params[[i]] <- tmp.drc.par
  
  ## temporary prediction
  tmp.drc.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence",
                             level = 0.95)
  rownames(tmp.drc.predict) <- c(0, exp(seq(log(0.001), log(5000), length = 999)))

  ## export
  PYL1.TECAN.curves.auc.predict[[i]] <- tmp.drc.predict
  
  ## clean up
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict.newdata, tmp.drc.predict)
  
}

## clean up
rm(packages,results,setup1,i,k,tmp.out,WT.PYL1.drc.1,WT.PYL1.drc.1.par)


## 2. Pre-processed DiMSum data ##
##################################

## Load raw measurements
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## Load dose-response curve metrics
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## Which TECAN variants have high conf. library fits?
parameters.Hill.lib <- parameters.Hill[match(names(PYL1.ABI1.TECAN.plate1.curves.auc), rownames(parameters.Hill)),]
parameters.Hill.lib.hq <- parameters.Hill.lib[which(parameters.Hill.lib[,"R^2"] > 0.9),]

## Filter the TECAN fits accordingly
PYL1.ABI1.TECAN.plate1.curves.auc <- PYL1.ABI1.TECAN.plate1.curves.auc[match(rownames(parameters.Hill.lib.hq), names(PYL1.ABI1.TECAN.plate1.curves.auc))]
PYL1.TECAN.curves.auc.predict <- PYL1.TECAN.curves.auc.predict[match(rownames(parameters.Hill.lib.hq), names(PYL1.TECAN.curves.auc.predict))]

## Summarise the key data
PYL1.ABI1.keyvars <- matrix(NA, ncol = 12, nrow = length(PYL1.ABI1.TECAN.plate1.curves.auc))
rownames(PYL1.ABI1.keyvars) <- names(PYL1.TECAN.curves.auc.predict)
colnames(PYL1.ABI1.keyvars) <- names(PYL1.ABI1)
for (i in 1:nrow(PYL1.ABI1.keyvars)){
  if(i == 1){
    PYL1.ABI1.keyvars[i,] <- sapply(PYL1.ABI1, function(x){x <- x[which(x$WT == T)[1],"gr_normalised_WTscaled"]; return(x)})
  }else{
    tmp.wt <- substr(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i], 1, 1)
    tmp.pos <- as.numeric(substr(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i], 2, nchar(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i]) - 1))
    tmp.mut <- substr(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i], nchar(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i]), nchar(names(PYL1.ABI1.TECAN.plate1.curves.auc)[i]))
    PYL1.ABI1.keyvars[i,] <- sapply(PYL1.ABI1, function(x){x <- x[which(x$WT_AA == tmp.wt & x$Pos == tmp.pos & x$Mut == tmp.mut),"gr_normalised_WTscaled"]; return(x)})
  }
}
 
## Remake the curves
PYL1.ABI1.DMS.curves <- PYL1.ABI1.DMS.curves.params <- vector(mode = "list", length = nrow(PYL1.ABI1.keyvars))
names(PYL1.ABI1.DMS.curves) <- names(PYL1.ABI1.DMS.curves.params) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
for(i in 1:length(PYL1.ABI1.DMS.curves)){

  ## temporary input data frame
  tmp.drc <- cbind(PYL1.ABI1.keyvars[i,],colnames(PYL1.ABI1.keyvars))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")
  
  ## temporary Hill fit & parameters 
  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')
  tmp.drc.par <- tmp.drc$fit$par
  names(tmp.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
  tmp.drc.par <- tmp.drc.par[c(2:4,1)]
  tmp.drc.par[4] <- -tmp.drc.par[4]
  PYL1.ABI1.DMS.curves.params[[i]] <- tmp.drc.par
  
  ## predict the full, smoothened curve (using 1000 data points) and confidence interval
  tmp.drc.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  PYL1.ABI1.DMS.curves[[i]] <- predict(tmp.drc,
                                       newdata = tmp.drc.predict.newdata,
                                       interval = "confidence",
                                       level = 0.95)
  rownames(PYL1.ABI1.DMS.curves[[i]]) <- c(0, exp(seq(log(0.001), log(5000), length = 999)))
  
  ## clean up
  rm(tmp.drc, tmp.drc.par, tmp.drc.predict.newdata)
  
}

## clean up
rm(i, tmp.mut, tmp.pos, tmp.wt, parameters.Hill, parameters.Hill.lib, parameters.Hill.lib.hq, PYL1.ABI1)


## 3. Plot the two curve sets ##
################################

drc.plot.microtiter <- function(curves, datapoints, variant){
  
  ## smooth curve prediction
  curves <- curves[[which(names(curves) == variant)]]
  class(curves) <- "numeric"
  curves <- as.data.frame(curves)
  curves$concentration <- as.numeric(rownames(curves))
  
  ## raw data points
  datapoints <- datapoints[[which(names(datapoints) == variant)]]
  datapoints <- cbind(dosages, datapoints)
  colnames(datapoints) <- c("concentration", "binding")
  class(datapoints) <- "numeric"
  datapoints <- as.data.frame(datapoints)
  datapoints[1,"concentration"] <- 9.062741e-03/3.5/3.5/3.5
  
  ## plot
  ggplot(curves, aes(x = concentration, y = Prediction)) +
    geom_ribbon(data = curves,
                aes(x = concentration, y = Prediction, ymin = Lower, ymax = Upper),
                alpha = 0.2, fill = alpha("darkgreen", alpha = 0.5)) +
    geom_point(data = datapoints, aes(x = concentration, y = binding),
               color = "black", size = 8, shape = 16) +
    geom_line(data = curves, aes(x = concentration, y = Prediction), linewidth = 1.5) +
    scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                  limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-1000, 1000)) +
    coord_cartesian(ylim = c(-15, 120)) +
    annotate("text",
             x = 9.062741e-03/3.5/3.5/3.5,
             y = 110,
             label = variant,
             hjust = 0, size = 25, color = "black") +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 40),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(size = 1, color = 'black'),
          axis.line.y = element_line(size = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 45, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 45, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "(+)-ABA conc. (µM)",
         y = "Relative PYL1/ABI1 Binding")
  
}
drc.plot.DMS <- function(curves, datapoints, variant){
  
  ## smooth curve prediction
  curves <- curves[[which(names(curves) == variant)]]
  class(curves) <- "numeric"
  curves <- as.data.frame(curves)
  curves$concentration <- as.numeric(rownames(curves))
  
  ## raw data points
  datapoints <- rev(datapoints[variant,])
  datapoints <- cbind(dosages, datapoints)
  colnames(datapoints) <- c("concentration", "binding")
  class(datapoints) <- "numeric"
  datapoints <- as.data.frame(datapoints)
  datapoints[1,"concentration"] <- 9.062741e-03/3.5/3.5/3.5
  
  ## plot
  ggplot(curves, aes(x = concentration, y = Prediction)) +
    geom_ribbon(data = curves,
                aes(x = concentration, y = Prediction, ymin = Lower, ymax = Upper),
                alpha = 0.2, fill = alpha("cornflowerblue", alpha = 0.5)) +
    geom_point(data = datapoints, aes(x = concentration, y = binding),
               color = "black", size = 8, shape = 16) +
    geom_line(data = curves, aes(x = concentration, y = Prediction), linewidth = 1.5) +
    scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                  limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-1000, 1000)) +
    coord_cartesian(ylim = c(-15, 120)) +
    annotate("text",
             x = 9.062741e-03/3.5/3.5/3.5,
             y = 110,
             label = variant,
             hjust = 0, size = 25, color = "black") +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(),
          plot.subtitle = element_markdown(size = 20),
          title = element_text(size = 40),
          axis.text = element_text(size = 30),
          axis.line.x = element_line(size = 1, color = 'black'),
          axis.line.y = element_line(size = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 45, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 45, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    labs(x = "(+)-ABA conc. (µM)",
         y = "Relative PYL1/ABI1 Binding")
  
}

pdf("../../results/FigureS2/FigureS2C_curves_microtiter_vs_DMS.pdf", height = 25, width = 100)
grid.arrange(drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "WT"), 
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "Q34R"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "Q34I"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "E36R"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "H87A"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "H87V"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "H87P"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "A116S"),
             drc.plot.microtiter(PYL1.TECAN.curves.auc.predict, PYL1.ABI1.TECAN.plate1.curves.auc, "S119A"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "WT"), 
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "Q34R"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "Q34I"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "E36R"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "H87A"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "H87V"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "H87P"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "A116S"),
             drc.plot.DMS(PYL1.ABI1.DMS.curves, PYL1.ABI1.keyvars, "S119A"),
             ncol = 9, nrow = 2, respect = F)
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
# [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] gridExtra_2.3      rlang_1.1.6        ggtext_0.1.2       ggplot2_4.0.0      scales_1.4.0       drc_3.0-1          MASS_7.3-65        reshape_0.8.10     growthcurver_0.3.1
# [10] readxl_1.4.5      
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       lattice_0.22-7     magrittr_2.0.4     RColorBrewer_1.1-3 mvtnorm_1.3-3      cellranger_1.1.0  
# [10] plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    TH.data_1.1-4      codetools_0.2-20   abind_1.4-8        cli_3.6.5         
# [19] crayon_1.5.3       splines_4.5.1      withr_3.0.2        plotrix_3.8-4      tools_4.5.1        minpack.lm_1.2-4   dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1          
# [28] zoo_1.8-14         lifecycle_1.0.4    car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         tibble_3.3.0      
# [37] tidyselect_1.2.1   farver_2.1.2       carData_3.0-5      compiler_4.5.1     S7_0.2.0           gridtext_0.1.5  
