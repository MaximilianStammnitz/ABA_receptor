# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 29.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################
## Figure 1E - individual variant validation ##
###############################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "growthcurver", "reshape", "drc", "ggplot2", "ggtext", "rlang")

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

## linearly rescale all curves to WT Bmax (100%) and 0%
coefs <- lm(c(0, 100) ~ c(mean(PYL1.ABI1.TECAN.plate1.curves.auc$L144A[1:10]), WT.PYL1.drc.1.par["B[inf]"]))
PYL1.ABI1.TECAN.plate1.curves.auc <- lapply(PYL1.ABI1.TECAN.plate1.curves.auc, function(x){y <- x*coefs$coefficients[[2]] + coefs$coefficients[[1]]; return(y)})

## Calculate dose-response curves using duplicate measurements, generate the confidence intervals at the key concentrations
PYL1.ABI1.TECAN.curves.auc.predict <- vector(mode = "list", length = length(PYL1.ABI1.TECAN.plate1.curves.auc))
names(PYL1.ABI1.TECAN.curves.auc.predict) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
for(i in 1:length(PYL1.ABI1.TECAN.plate1.curves.auc)){

  print(names(PYL1.ABI1.TECAN.curves.auc.predict)[i])
    
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

  ## temporary prediction
  tmp.drc.predict.newdata <- expand.grid(conc = dosages)
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence",
                             level = 0.95)
  rownames(tmp.drc.predict) <- dosages

  ## export
  PYL1.ABI1.TECAN.curves.auc.predict[[i]] <- tmp.drc.predict
  
  ## clean up
  rm(tmp.drc, tmp.drc.predict.newdata, tmp.drc.predict)
  
}

## clean up
rm(results,setup1,tecan.30,i,k,tmp.out,WT.PYL1.drc,WT.PYL1.drc.par)


## 2. Pre-processed DiMSum data ##
##################################

## load pre-processed GluePCA DMS data 
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## load dose-response curve metrics
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## Which TECAN variants have high conf. library fits?
parameters.Hill.lib <- parameters.Hill[match(names(PYL1.ABI1.TECAN.plate1.curves.auc), rownames(parameters.Hill)),]
parameters.Hill.lib.hq <- parameters.Hill.lib[which(parameters.Hill.lib[,"R^2"] > 0.9),]

## Filter the TECAN fits accordingly
PYL1.ABI1.TECAN.plate1.curves.auc <- PYL1.ABI1.TECAN.plate1.curves.auc[match(rownames(parameters.Hill.lib.hq), names(PYL1.ABI1.TECAN.plate1.curves.auc))]
PYL1.ABI1.TECAN.curves.auc.predict <- PYL1.ABI1.TECAN.curves.auc.predict[match(rownames(parameters.Hill.lib.hq), names(PYL1.ABI1.TECAN.curves.auc.predict))]

## Summarise the key data
PYL1.ABI1.keyvars <- matrix(NA, ncol = 12, nrow = length(PYL1.ABI1.TECAN.plate1.curves.auc))
rownames(PYL1.ABI1.keyvars) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
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
 
## Remake the curves, generate the confidence intervals at the key concentrations
PYL1.ABI1.keyvars.out <- vector(mode = "list", length = nrow(PYL1.ABI1.keyvars))
names(PYL1.ABI1.keyvars.out) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
for(i in 1:length(PYL1.ABI1.keyvars.out)){

  tmp.drc <- cbind(PYL1.ABI1.keyvars[i,],colnames(PYL1.ABI1.keyvars))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")

  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')

  ## Predict the full, smoothened curve (using 1000 data points) and confidence interval
  tmp.predict.newdata <- expand.grid(conc = as.numeric(colnames(PYL1.ABI1.keyvars)))
  PYL1.ABI1.keyvars.out[[i]] <- predict(tmp.drc,
                                        newdata = tmp.predict.newdata,
                                        interval = "confidence",
                                        level = 0.95)
  rownames(PYL1.ABI1.keyvars.out[[i]]) <- as.numeric(colnames(PYL1.ABI1.keyvars))

}


## 3. Combine two the data sets ##
##################################

## integrate data
PYL1.mutant.summary <- vector(mode = "list", length = length(PYL1.ABI1.TECAN.plate1.curves.auc))
names(PYL1.mutant.summary) <- names(PYL1.ABI1.TECAN.plate1.curves.auc)
for(i in 1:length(PYL1.mutant.summary)){
  
  #### matrix setup
  PYL1.mutant.summary[[i]] <- matrix(NA, ncol = 8, nrow = 12)
  colnames(PYL1.mutant.summary[[i]]) <- c("TECAN",
                                          "TECAN_predict",
                                          "TECAN_predict_min",
                                          "TECAN_predict_max",
                                          "Competition_gr",
                                          "Competition_predict",
                                          "Competition_predict_min",
                                          "Competition_predict_max")
  rownames(PYL1.mutant.summary[[i]]) <- sort(dosages)
  
  #### fill in TECAN data
  PYL1.mutant.summary[[i]][,"TECAN"] <- PYL1.ABI1.TECAN.plate1.curves.auc[[i]]
  PYL1.mutant.summary[[i]][,"TECAN_predict"] <- PYL1.ABI1.TECAN.curves.auc.predict[[i]][,"Prediction"]
  PYL1.mutant.summary[[i]][,"TECAN_predict_min"] <- PYL1.ABI1.TECAN.curves.auc.predict[[i]][,"Lower"]
  PYL1.mutant.summary[[i]][,"TECAN_predict_max"] <- PYL1.ABI1.TECAN.curves.auc.predict[[i]][,"Upper"]
  
  #### fill in competition data
  if(i == 1){
    PYL1.mutant.summary[[i]][,"Competition_gr"] <- rev(apply(sapply(PYL1.ABI1, function(x){out <- x[which(x$WT == T),"gr_normalised_WTscaled"]; return(out)}), 2, median))
    PYL1.mutant.summary[[i]][,"Competition_predict"] <- rev(PYL1.ABI1.keyvars.out[[i]][,1])
    PYL1.mutant.summary[[i]][,"Competition_predict_min"] <- rev(PYL1.ABI1.keyvars.out[[i]][,2])
    PYL1.mutant.summary[[i]][,"Competition_predict_max"] <- rev(PYL1.ABI1.keyvars.out[[i]][,3])
  }else{
    var.tmp <- names(PYL1.mutant.summary)[i]
    var.tmp.pos <- as.numeric(substring(var.tmp, 2, nchar(var.tmp)-1))
    var.tmp.mut <- substring(var.tmp, nchar(var.tmp))
    PYL1.mutant.summary[[i]][,"Competition_gr"] <- rev(sapply(sapply(PYL1.ABI1, function(x){out <- x[which(x$Pos == var.tmp.pos & x$Mut == var.tmp.mut),"gr_normalised_WTscaled"]; return(out)}), median))
    PYL1.mutant.summary[[i]][,"Competition_predict"] <- rev(PYL1.ABI1.keyvars.out[[i]][,1])
    PYL1.mutant.summary[[i]][,"Competition_predict_min"] <- rev(PYL1.ABI1.keyvars.out[[i]][,2])
    PYL1.mutant.summary[[i]][,"Competition_predict_max"] <- rev(PYL1.ABI1.keyvars.out[[i]][,3])
  }
  
}

## summarise
TECAN <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,1]}))
TECAN.pred <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,2]}))
TECAN.pred.min <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,3]}))
TECAN.pred.max <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,4]}))

DMS <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,5]}))
DMS.pred <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,6]}))
DMS.pred.min <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,7]}))
DMS.pred.max <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,8]}))

out.tecan.df <- cbind("conc" = rep(1:12, 9),
                      "tecan mean" = TECAN, 
                      "tecan pred" = TECAN.pred, 
                      "tecan pred min" = TECAN.pred.min, 
                      "tecan pred max" = TECAN.pred.max, 
                      "competition" = DMS, 
                      "competition pred" = DMS.pred, 
                      "competition pred min" = DMS.pred.min, 
                      "competition pred max" = DMS.pred.max)
out.tecan.df <- as.data.frame(out.tecan.df)

## correlation
r <- cor(x = out.tecan.df$`tecan pred`,
         y = out.tecan.df$`competition pred`,
         method = "pearson")

## plot
pdf("../../results/Figure1/Figure1E_individual_mutant_correlation.pdf",
    height = 15, width = 18)

out.1E <- ggplot(out.tecan.df, aes(x = `tecan pred`, y = `competition pred`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 100, length.out = 6), limits = c(-1000, 1000)) +
  coord_cartesian(xlim = c(-5, 121), ylim = c(-5, 115)) +
  geom_smooth(out.tecan.df,
              mapping = aes(x = `tecan pred`, y = `competition pred`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = "grey50") +
  geom_errorbar(data = out.tecan.df,
                aes(x = `tecan pred`,
                    ymin = `competition pred min`,
                    ymax = `competition pred max`),
                linewidth = 0.2,
                color = "black") +
  geom_errorbar(data = out.tecan.df,
                 aes(x = `competition pred`, 
                     xmin = `tecan pred min`, 
                     xmax = `tecan pred max`),
                 linewidth = 0.2,
                 color = "black") +
  geom_point(data = out.tecan.df,
             mapping = aes(x = `tecan pred`, y = `competition pred`, fill = `conc`),
             color = "black", size = 5, shape = 21, stroke = 0.5) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
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
  labs(x = "Relative PYL1/ABI1 Binding (individual mutants)",
       y = "Relative PYL1/ABI1 Binding (library)")

print(out.1E)

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
# [1] rlang_1.1.6        ggtext_0.1.2       ggplot2_4.0.0      drc_3.0-1          MASS_7.3-65        reshape_0.8.10     growthcurver_0.3.1
# [8] readxl_1.4.5      
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       lattice_0.22-7     magrittr_2.0.4     grid_4.5.1        
# [8] RColorBrewer_1.1-3 mvtnorm_1.3-3      cellranger_1.1.0   plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3    
# [15] multcomp_1.4-28    mgcv_1.9-3         scales_1.4.0       TH.data_1.1-4      codetools_0.2-20   abind_1.4-8        cli_3.6.5         
# [22] crayon_1.5.3       splines_4.5.1      withr_3.0.2        plotrix_3.8-4      tools_4.5.1        minpack.lm_1.2-4   dplyr_1.1.4       
# [29] vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4    car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1     
# [36] gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         tibble_3.3.0       tidyselect_1.2.1   farver_2.1.2       nlme_3.1-168      
# [43] labeling_0.4.3     carData_3.0-5      compiler_4.5.1     S7_0.2.0           gridtext_0.1.5   
