# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 29.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

####################################################################
## Supplementary Figure 3B - PYL1-ABI1 Hill parameter validations ##
####################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthcurver", "drc", 
              "scales", "reshape", "ggplot2", "ggtext", "cowplot",
              "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input PYL1-ABI1 data ##
#############################

## input Hill parameter distributions from dose-response curve fits
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")
parameters.Hill <- parameters.Hill[!is.na(parameters.Hill[,1]),]


## 2. Import first TECAN plate ##
#################################

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


## 3. Import second TECAN plate ##
##################################

## raw data import
PYL1.ABI1.TECAN.plate2 <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[20:405,1:263]
rownames(PYL1.ABI1.TECAN.plate2) <- PYL1.ABI1.TECAN.plate2[,1]
colnames(PYL1.ABI1.TECAN.plate2) <- PYL1.ABI1.TECAN.plate2[1,]
PYL1.ABI1.TECAN.plate2 <- PYL1.ABI1.TECAN.plate2[,-1]
PYL1.ABI1.TECAN.plate2 <- PYL1.ABI1.TECAN.plate2[-c(1:2),]
colnames(PYL1.ABI1.TECAN.plate2) <- as.numeric(colnames(PYL1.ABI1.TECAN.plate2))/3600 ## convert to hours
class(PYL1.ABI1.TECAN.plate2) <- "numeric"
PYL1.ABI1.TECAN.plate2 <- as.data.frame(PYL1.ABI1.TECAN.plate2)

## assign mutants
setup2 <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[1:16,1:25]
rownames(setup2) <- setup2[,1]
setup2 <- setup2[,-1]
colnames(setup2) <- 1:24


## 4. Process using growthcurver ##
###################################

## shape format
PYL1.ABI1.TECAN.plate1 <- cbind("time" = colnames(PYL1.ABI1.TECAN.plate1), 
                                t(PYL1.ABI1.TECAN.plate1))
class(PYL1.ABI1.TECAN.plate1) <- "numeric"
PYL1.ABI1.TECAN.plate1 <- as.data.frame(PYL1.ABI1.TECAN.plate1)

PYL1.ABI1.TECAN.plate2 <- cbind("time" = colnames(PYL1.ABI1.TECAN.plate2), 
                                t(PYL1.ABI1.TECAN.plate2))
class(PYL1.ABI1.TECAN.plate2) <- "numeric"
PYL1.ABI1.TECAN.plate2 <- as.data.frame(PYL1.ABI1.TECAN.plate2)

## run
PYL1.ABI1.TECAN.plate1.curves <- SummarizeGrowthByPlate(PYL1.ABI1.TECAN.plate1[1:262,], 
                                                        t_trim = 70,
                                                        bg_correct = "min",
                                                        plot_fit = F)

PYL1.ABI1.TECAN.plate2.curves <- SummarizeGrowthByPlate(PYL1.ABI1.TECAN.plate2[1:262,], 
                                                        t_trim = 70,
                                                        bg_correct = "min",
                                                        plot_fit = F)

## label wells by genotype
for(i in 1:nrow(PYL1.ABI1.TECAN.plate1.curves)){
  PYL1.ABI1.TECAN.plate1.curves[i,1] <- setup1[substr(PYL1.ABI1.TECAN.plate1.curves[i,1], 1, 1),
                                               substr(PYL1.ABI1.TECAN.plate1.curves[i,1], 2, 
                                                      nchar(PYL1.ABI1.TECAN.plate1.curves[i,1]))]
}

for(i in 1:nrow(PYL1.ABI1.TECAN.plate2.curves)){
  PYL1.ABI1.TECAN.plate2.curves[i,1] <- setup2[substr(PYL1.ABI1.TECAN.plate2.curves[i,1], 1, 1),
                                       substr(PYL1.ABI1.TECAN.plate2.curves[i,1], 2, 
                                              nchar(PYL1.ABI1.TECAN.plate2.curves[i,1]))]
}

## categorise
PYL1.ABI1.TECAN.plate1.curves.auc <- vector(mode = "list", length = length(unique(c(setup1[-c(1,17,18),]))))
names(PYL1.ABI1.TECAN.plate1.curves.auc) <- unique(c(setup1[-c(1,17,18),]))
names(PYL1.ABI1.TECAN.plate1.curves.auc) <- c(names(PYL1.ABI1.TECAN.plate1.curves.auc)[1],
                                              names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1][order(as.numeric(substr(names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1], 2, nchar(names(PYL1.ABI1.TECAN.plate1.curves.auc)[-1]) - 1)))])
PYL1.ABI1.TECAN.plate1.curves.auc <- PYL1.ABI1.TECAN.plate1.curves.auc[-which(names(PYL1.ABI1.TECAN.plate1.curves.auc) %in% c("T118N", "S119A", "S119N"))]
for(i in 1:length(PYL1.ABI1.TECAN.plate1.curves.auc)){
  PYL1.ABI1.TECAN.plate1.curves.auc[[i]] <- which(PYL1.ABI1.TECAN.plate1.curves$sample == names(PYL1.ABI1.TECAN.plate1.curves.auc)[i])
  PYL1.ABI1.TECAN.plate1.curves.auc[[i]] <- PYL1.ABI1.TECAN.plate1.curves$auc_e[PYL1.ABI1.TECAN.plate1.curves.auc[[i]]]
}
PYL1.ABI1.TECAN.plate1.curves.auc <- lapply(PYL1.ABI1.TECAN.plate1.curves.auc, rev)

## take mean value of replicates
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

PYL1.ABI1.TECAN.plate2.curves.auc <- vector(mode = "list", length = 31)
names(PYL1.ABI1.TECAN.plate2.curves.auc) <- unique(c(setup2))
names(PYL1.ABI1.TECAN.plate2.curves.auc) <- c(names(PYL1.ABI1.TECAN.plate2.curves.auc)[1],
                                              names(PYL1.ABI1.TECAN.plate2.curves.auc)[-1][order(as.numeric(substr(names(PYL1.ABI1.TECAN.plate2.curves.auc)[-1], 2, nchar(names(PYL1.ABI1.TECAN.plate2.curves.auc)[-1]) - 1)))])
for(i in 1:length(PYL1.ABI1.TECAN.plate2.curves.auc)){
  PYL1.ABI1.TECAN.plate2.curves.auc[[i]] <- which(PYL1.ABI1.TECAN.plate2.curves$sample == names(PYL1.ABI1.TECAN.plate2.curves.auc)[i])
  PYL1.ABI1.TECAN.plate2.curves.auc[[i]] <- PYL1.ABI1.TECAN.plate2.curves$auc_e[PYL1.ABI1.TECAN.plate2.curves.auc[[i]]]
}

## dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
dosages <- rev(dosages)

## fit plate-specific WT curves
WT.PYL1.drc.1 <- cbind(PYL1.ABI1.TECAN.plate1.curves.auc$WT,dosages)
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

WT.PYL1.drc.2 <- cbind(PYL1.ABI1.TECAN.plate2.curves.auc$WT[13:24], dosages)
class(WT.PYL1.drc.2) <- "numeric"
WT.PYL1.drc.2 <- as.data.frame(WT.PYL1.drc.2)
colnames(WT.PYL1.drc.2) <- c("GR", "concentration")
WT.PYL1.drc.2 <- drm(WT.PYL1.drc.2$GR ~ WT.PYL1.drc.2$concentration,
                     fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                     type = 'continuous')
WT.PYL1.drc.2.par <- WT.PYL1.drc.2$fit$par
names(WT.PYL1.drc.2.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.2.par <- WT.PYL1.drc.2.par[c(2:4,1)]
WT.PYL1.drc.2.par[4] <- -WT.PYL1.drc.2.par[4]

## linearly rescale all curves to WT Bmax (100%) and 0%
coefs1 <- lm(c(0, 100) ~ c(mean(PYL1.ABI1.TECAN.plate1.curves.auc$L144A[1:10]), WT.PYL1.drc.1.par["B[inf]"]))
PYL1.ABI1.TECAN.plate1.curves.auc <- lapply(PYL1.ABI1.TECAN.plate1.curves.auc, function(x){y <- x*coefs1$coefficients[[2]] + coefs1$coefficients[[1]]; return(y)})
coefs2 <- lm(c(0, 100) ~ c(mean(PYL1.ABI1.TECAN.plate2.curves.auc$L144A[1:10]), WT.PYL1.drc.2.par["B[inf]"]))
PYL1.ABI1.TECAN.plate2.curves.auc <- lapply(PYL1.ABI1.TECAN.plate2.curves.auc, function(x){y <- x*coefs2$coefficients[[2]] + coefs2$coefficients[[1]]; return(y)})


## 5. Unify TECAN data ##
#########################

## from second data set only take the ones not in the first
PYL1.ABI1.TECAN.curves.auc <- c(PYL1.ABI1.TECAN.plate1.curves.auc,
                                PYL1.ABI1.TECAN.plate2.curves.auc[-which(names(PYL1.ABI1.TECAN.plate2.curves.auc) %in% names(PYL1.ABI1.TECAN.plate1.curves.auc))])
PYL1.ABI1.TECAN.curves.auc <- c(PYL1.ABI1.TECAN.curves.auc[1],
                                PYL1.ABI1.TECAN.curves.auc[order(as.numeric(substr(names(PYL1.ABI1.TECAN.curves.auc)[-1], 2, nchar(names(PYL1.ABI1.TECAN.curves.auc)[-1]) - 1))) + 1])


## 6. Fit all TECAN curves' dose response profiles ##
#####################################################

parameters.Hill.TECAN <- matrix(NA, nrow = length(PYL1.ABI1.TECAN.curves.auc), ncol = 16)
colnames(parameters.Hill.TECAN) <- c("Hill", "B[0]", "B[inf]", "EC50", 
                                     "Hill SE", "B[0] SE", "B[inf] SE", "EC50 SE", 
                                     "Hill P", "B[0] P", "B[inf] P", "EC50 P", 
                                     "Data points", "AIC", "Residual var", "R^2")
rownames(parameters.Hill.TECAN) <- names(PYL1.ABI1.TECAN.curves.auc)

## Calculate and plot all the dose-response curves with the 4-parametric Hill model
for(i in 1:nrow(parameters.Hill.TECAN)){
  
  print(i)
  
  ### DRC modelling, using the R DRM package
  
  ### output parameter translations:
  ### b: Hill coefficient, i.e. steepness of the curve (n)
  ### c: basal binding fitness (B[0])
  ### d: saturated binding fitness (B[inf])
  ### e: curve inflection point (EC50)
  ### model: binding(ABA conc.) = c + ((d-c)/(1 + (e/x)^b))
  
  ### fetch genotype's data
  tmp.PYL1.drc.in <- cbind(PYL1.ABI1.TECAN.curves.auc[[i]][1:12],dosages)
  class(tmp.PYL1.drc.in) <- "numeric"
  tmp.PYL1.drc.in <- as.data.frame(tmp.PYL1.drc.in)
  colnames(tmp.PYL1.drc.in) <- c("GR", "concentration")
  
  ## ignore questionable fits and 'bandstops'
  if(!rownames(parameters.Hill.TECAN)[i] %in% c("T118N", "S119N", "V110H", "A116H", "V193H", "V193W")){

    ### curve fit
    tmp.PYL1.drc <- drm(tmp.PYL1.drc.in$GR ~ tmp.PYL1.drc.in$concentration,
                        fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                        type = 'continuous')
    
    ### parameter take-over
    parameters.Hill.TECAN[i,"Data points"] <- sum(is.na(tmp.PYL1.drc.in$GR) == F)
    parameters.Hill.TECAN[i,1:4] <- summary(tmp.PYL1.drc)$coefficients[,1]
    parameters.Hill.TECAN[i,5:8] <- summary(tmp.PYL1.drc)$coefficients[,2]
    parameters.Hill.TECAN[i,9:12] <- summary(tmp.PYL1.drc)$coefficients[,4]
    if(parameters.Hill.TECAN[i,1] < 0){
      parameters.Hill.TECAN[i,1] <- -parameters.Hill.TECAN[i,1] ## invert Hill parameter output
    }else{
      parameters.Hill.TECAN[i,2:3] <- parameters.Hill.TECAN[i,c(3,2)]
      parameters.Hill.TECAN[i,10:11] <- parameters.Hill.TECAN[i,c(11,10)]
    }
    parameters.Hill.TECAN[i,14:15] <- as.numeric(mselect(tmp.PYL1.drc, icfct = AIC)[c(2,4)])
    
    ### calculation of non-linear R^2
    predicted_values <- predict(tmp.PYL1.drc)
    rss <- sum((tmp.PYL1.drc.in$GR - predicted_values)^2)
    tss <- sum((tmp.PYL1.drc.in$GR - mean(tmp.PYL1.drc.in$GR))^2)
    tmp.r_squared.Hill <- 1 - (rss/tss)
    parameters.Hill.TECAN[i,16] <- tmp.r_squared.Hill
    rm(predicted_values, rss, tss, tmp.r_squared.Hill)
   
    ### predict the full, smoothened curve (using 1000 data points) and confidence interval
    tmp.PYL1.drc.Hill.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
    tmp.PYL1.drc.Hill.predict <- predict(tmp.PYL1.drc, 
                                         newdata = tmp.PYL1.drc.Hill.predict.newdata, 
                                         interval = "confidence")
    
    #### new data with predictions
    tmp.PYL1.drc.Hill.predict.newdata$p <- tmp.PYL1.drc.Hill.predict[,1]
    tmp.PYL1.drc.Hill.predict.newdata$pmin <- tmp.PYL1.drc.Hill.predict[,2]
    tmp.PYL1.drc.Hill.predict.newdata$pmax <- tmp.PYL1.drc.Hill.predict[,3]
    
    ### plot curve
    tmp.PYL1.drc.in$concentration[1] <- 9.062741e-03/3.5/3.5/3.5 ## "0-conc." positioning for log scale
    rm(tmp.PYL1.drc, tmp.PYL1.drc.Hill.predict, tmp.PYL1.drc.Hill.predict.newdata, tmp.PYL1.drc.in)
    
  }

}


## 7. Plot ##
#############

## summarise data
out.all <- matrix(NA, nrow = nrow(parameters.Hill.TECAN), ncol = 26)
colnames(out.all) <- c("bulk_B0", "bulk_B0_SE", "bulk_B0_P", 
                       "bulk_Binf", "bulk_Binf_SE", "bulk_Binf_P", 
                       "bulk_EC50", "bulk_EC50_SE", "bulk_EC50_P", 
                       "bulk_Hill", "bulk_Hill_SE", "bulk_Hill_P",
                       "bulk_R2",
                       "TECAN_B0", "TECAN_B0_SE", "TECAN_B0_P", 
                       "TECAN_Binf", "TECAN_Binf_SE", "TECAN_Binf_P", 
                       "TECAN_EC50", "TECAN_EC50_SE", "TECAN_EC50_P", 
                       "TECAN_Hill", "TECAN_Hill_SE", "TECAN_Hill_P",
                       "TECAN_R2")
rownames(out.all) <- rownames(parameters.Hill.TECAN)

for(i in 1:nrow(out.all)){
  tmp.name <- rownames(out.all)[i]
  tmp.id1 <- match(tmp.name, rownames(parameters.Hill))
  out1 <- parameters.Hill[tmp.id1,c(2,6,10,3,7,11,4,8,12,1,5,9,16)]
  tmp.id2 <- match(tmp.name, rownames(parameters.Hill.TECAN))
  out2 <- parameters.Hill.TECAN[tmp.id2,c(2,6,10,3,7,11,4,8,12,1,5,9,16)]
  out.all[i,] <- c(out1, out2)
}
out.all <- as.data.frame(out.all)

## raw linear regression
out.B0 <- out.all[which(out.all$bulk_B0_P < 0.1 & out.all$bulk_R2 > 0.5 & out.all$TECAN_B0_P < 0.1 & out.all$TECAN_R2 > 0.5),]
p.B0 <- summary(lm(out.B0$bulk_B0 ~ out.B0$TECAN_B0))$coefficients[2,4]
out.Binf <- out.all[which(out.all$bulk_Binf_P < 0.1 & out.all$bulk_R2 > 0.5 & out.all$TECAN_Binf_P < 0.1 & out.all$TECAN_R2 > 0.5),]
p.Binf <- summary(lm(out.Binf$bulk_Binf ~ out.Binf$TECAN_Binf))$coefficients[2,4]
out.EC50 <- out.all[which(out.all$bulk_EC50_P < 0.1 & out.all$bulk_R2 > 0.5 & out.all$TECAN_EC50_P < 0.1 & out.all$TECAN_R2 > 0.5),]
p.EC50 <- summary(lm(log10(out.EC50$bulk_EC50) ~ log10(out.EC50$TECAN_EC50)))$coefficients[2,4]
out.Hill <- out.all[which(out.all$bulk_Hill_P < 0.1 & out.all$bulk_R2 > 0.5 & out.all$TECAN_Hill_P < 0.1 & out.all$TECAN_R2 > 0.5),]
p.Hill <- summary(lm(log10(out.Hill$bulk_Hill) ~ log10(out.Hill$TECAN_Hill)))$coefficients[2,4]

## B0
r.B0 <- cor(x = out.B0$bulk_B0, y = out.B0$TECAN_B0, 
            method = "pearson", use = "complete.obs")

out.S3B_B0 <- ggplot(out.B0, aes(x = `TECAN_B0`, y = `bulk_B0`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 130)) +
  geom_point(data = out.B0,
             mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = `TECAN_B0` - `TECAN_B0_SE`,
                    xmax = `TECAN_B0` + `TECAN_B0_SE`), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = `bulk_B0` - `bulk_B0_SE`,
                    ymax = `bulk_B0` + `bulk_B0_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.B0,
              mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = -20,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.B0, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.B0, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[0] ~ "(individual mutants, AUC)"),
       y = bquote(B[0] ~ "(library sequencing)"))

## Binf
r.Binf <- cor(x = out.Binf$bulk_Binf, y = out.Binf$TECAN_Binf, 
              method = "pearson", use = "complete.obs")

out.S3B_Binf <- ggplot(out.Binf, aes(x = `TECAN_Binf`, y = `bulk_Binf`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-20, 130), ylim = c(-20, 130)) +
  geom_point(data = out.Binf,
             mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = `TECAN_Binf` - `TECAN_Binf_SE`,
                    xmax = `TECAN_Binf` + `TECAN_Binf_SE`), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = `bulk_Binf` - `bulk_Binf_SE`,
                    ymax = `bulk_Binf` + `bulk_Binf_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.Binf,
              mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = -20,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.Binf, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.Binf, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[infinity] ~ "(individual mutants, AUC)"),
       y = bquote(B[infinity] ~ "(library sequencing)"))

## EC50
r.EC50 <- cor(x = log(out.EC50$bulk_EC50), y = log(out.EC50$TECAN_EC50), 
              method = "pearson", use = "complete.obs")

out.S3B_EC50 <- ggplot(out.EC50, aes(x = `TECAN_EC50`, y = `bulk_EC50`)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.001, 10000), ylim = c(0.001, 10000), expand = T) +
  geom_point(data = out.EC50,
             mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = TECAN_EC50 - TECAN_EC50_SE, 
                    xmax = TECAN_EC50 + TECAN_EC50_SE), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = `bulk_EC50` - `bulk_EC50_SE`, 
                    ymax = `bulk_EC50` + `bulk_EC50_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.EC50,
              mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = 0.0018,
           y = 3025,
           label = expr_text(bquote(italic(r) == .(format(r.EC50, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.EC50, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(EC[50] ~ "(individual mutants, AUC)"),
       y = bquote(EC[50] ~ "(library sequencing)"))

## n
r.Hill <- cor(x = log(out.Hill$bulk_Hill), y = log(out.Hill$TECAN_Hill), 
              method = "pearson", use = "complete.obs")

out.S3B_Hill <- ggplot(out.Hill, aes(x = `TECAN_Hill`, y = `bulk_Hill`)) +
  scale_x_log10(breaks = c(0.4, 0.6, 1, 2, 4), 
                labels = c(0.4, 0.6, 1, 2, 4),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.4, 0.6, 1, 2, 4), 
                labels = c(0.4, 0.6, 1, 2, 4),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.4, 5), ylim = c(0.4, 5), expand = T) +
  geom_point(data = out.Hill,
             mapping = aes(x = `TECAN_Hill`, y = `bulk_Hill`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = pmax(TECAN_Hill - TECAN_Hill_SE, 0.1),
                    xmax = TECAN_Hill + TECAN_Hill_SE), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = pmax(`bulk_Hill` - `bulk_Hill_SE`, 0.1),
                    ymax = `bulk_Hill` + `bulk_Hill_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.Hill,
              mapping = aes(x = `TECAN_Hill`, y = `bulk_Hill`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  annotate("text",
           x = 0.4,
           y = 4.22,
           label = expr_text(bquote(italic(r) == .(format(r.Hill, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.Hill, digits = 4, nsmall = 4)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(n ~ "(individual mutants, AUC)"),
       y = bquote(n ~ "(library sequencing)"))

## combined plot
pdf("../../results/FigureS3/FigureS3B_Hill_parameter_validation.pdf", height = 15, width = 55)
plot_grid(out.S3B_B0, out.S3B_Binf, out.S3B_EC50, out.S3B_Hill, align = "hv", axis = "tblr", ncol = 4)
dev.off()


## 8. Version ##
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
# [1] rlang_1.1.6        cowplot_1.2.0      ggtext_0.1.2       ggplot2_4.0.0      reshape_0.8.10     scales_1.4.0       drc_3.0-1          MASS_7.3-65        growthcurver_0.3.1
# [10] readxl_1.4.5       stringr_1.5.2     
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      lattice_0.22-7     magrittr_2.0.4     grid_4.5.1         RColorBrewer_1.1-3
# [10] mvtnorm_1.3-3      cellranger_1.1.0   plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    mgcv_1.9-3         TH.data_1.1-4     
# [19] codetools_0.2-20   abind_1.4-8        cli_3.6.5          crayon_1.5.3       splines_4.5.1      withr_3.0.2        plotrix_3.8-4      tools_4.5.1        minpack.lm_1.2-4  
# [28] dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4    car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6      
# [37] glue_1.8.0         Rcpp_1.1.0         tibble_3.3.0       tidyselect_1.2.1   farver_2.1.2       nlme_3.1-168       carData_3.0-5      compiler_4.5.1     S7_0.2.0          
# [46] gridtext_0.1.5 
