# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 24.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################
## Figure 1E - individual variant validation ##
###############################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "growthrates", "reshape", "drc", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. TECAN data pre-processing ##
##################################

## import yeast colony matrix
setup <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[2:18,1:25]
colnames(setup) <- setup[1,]
setup <- setup[-1,]
rownames(setup) <- setup[,1]
setup <- setup[,-1]

## import results
tecan.30 <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[20:381,1:301]
rownames(tecan.30) <- tecan.30[,1]
tecan.30 <- tecan.30[,-1]

### Summarise all Rapamycine results
results <- vector(mode = 'list', length = 15)
names(results) <- c('WT',
                    'Q34R',
                    'Q34I',
                    'E36R',
                    'H87A',
                    'H87V',
                    'H87P',
                    'A116S',
                    'T118W',
                    'S119A',
                    'S119V',
                    'S119N',
                    'H142C',
                    'L144A')

results$WT <- tecan.30[c('B24', 'B23',
                         'B22', 'B21',
                         'B20', 'B19',
                         'B18', 'B17',
                         'B16', 'B15',
                         'B14', 'B13',
                         'B12', 'B11',
                         'B10', 'B9',
                         'B8', 'B7',
                         'B6', 'B5',
                         'B4', 'B3',
                         'B2', 'B1'),]

results$Q34R <- tecan.30[c('C24', 'C23',
                           'C22', 'C21',
                           'C20', 'C19',
                           'C18', 'C17',
                           'C16', 'C15',
                           'C14', 'C13',
                           'C12', 'C11',
                           'C10', 'C9',
                           'C8', 'C7',
                           'C6', 'C5',
                           'C4', 'C3',
                           'C2', 'C1'),]

results$Q34I <- tecan.30[c('D24', 'D23',
                           'D22', 'D21',
                           'D20', 'D19',
                           'D18', 'D17',
                           'D16', 'D15',
                           'D14', 'D13',
                           'D12', 'D11',
                           'D10', 'D9',
                           'D8', 'D7',
                           'D6', 'D5',
                           'D4', 'D3',
                           'D2', 'D1'),]

results$E36R <- tecan.30[c('E24', 'E23',
                           'E22', 'E21',
                           'E20', 'E19',
                           'E18', 'E17',
                           'E16', 'E15',
                           'E14', 'E13',
                           'E12', 'E11',
                           'E10', 'E9',
                           'E8', 'E7',
                           'E6', 'E5',
                           'E4', 'E3',
                           'E2', 'E1'),]

results$H87A <- tecan.30[c('F24', 'F23',
                           'F22', 'F21',
                           'F20', 'F19',
                           'F18', 'F17',
                           'F16', 'F15',
                           'F14', 'F13',
                           'F12', 'F11',
                           'F10', 'F9',
                           'F8', 'F7',
                           'F6', 'F5',
                           'F4', 'F3',
                           'F2', 'F1'),]

results$H87V <- tecan.30[c('G24', 'G23',
                           'G22', 'G21',
                           'G20', 'G19',
                           'G18', 'G17',
                           'G16', 'G15',
                           'G14', 'G13',
                           'G12', 'G11',
                           'G10', 'G9',
                           'G8', 'G7',
                           'G6', 'G5',
                           'G4', 'G3',
                           'G2', 'G1'),]

results$H87P <- tecan.30[c('H24', 'H23',
                           'H22', 'H21',
                           'H20', 'H19',
                           'H18', 'H17',
                           'H16', 'H15',
                           'H14', 'H13',
                           'H12', 'H11',
                           'H10', 'H9',
                           'H8', 'H7',
                           'H6', 'H5',
                           'H4', 'H3',
                           'H2', 'H1'),]

results$A116S <- tecan.30[c('I24', 'I23',
                            'I22', 'I21',
                            'I20', 'I19',
                            'I18', 'I17',
                            'I16', 'I15',
                            'I14', 'I13',
                            'I12', 'I11',
                            'I10', 'I9',
                            'I8', 'I7',
                            'I6', 'I5',
                            'I4', 'I3',
                            'I2', 'I1'),]

results$T118W <- tecan.30[c('K24', 'K23',
                            'K22', 'K21',
                            'K20', 'K19',
                            'K18', 'K17',
                            'K16', 'K15',
                            'K14', 'K13',
                            'K12', 'K11',
                            'K10', 'K9',
                            'K8', 'K7',
                            'K6', 'K5',
                            'K4', 'K3',
                            'K2', 'K1'),]

results$S119A <- tecan.30[c('L24', 'L23',
                            'L22', 'L21',
                            'L20', 'L19',
                            'L18', 'L17',
                            'L16', 'L15',
                            'L14', 'L13',
                            'L12', 'L11',
                            'L10', 'L9',
                            'L8', 'L7',
                            'L6', 'L5',
                            'L4', 'L3',
                            'L2', 'L1'),]

results$S119V <- tecan.30[c('M24', 'M23',
                            'M22', 'M21',
                            'M20', 'M19',
                            'M18', 'M17',
                            'M16', 'M15',
                            'M14', 'M13',
                            'M12', 'M11',
                            'M10', 'M9',
                            'M8', 'M7',
                            'M6', 'M5',
                            'M4', 'M3',
                            'M2', 'M1'),]

results$S119N <- tecan.30[c('N24', 'N23',
                            'N22', 'N21',
                            'N20', 'N19',
                            'N18', 'N17',
                            'N16', 'N15',
                            'N14', 'N13',
                            'N12', 'N11',
                            'N10', 'N9',
                            'N8', 'N7',
                            'N6', 'N5',
                            'N4', 'N3',
                            'N2', 'N1'),]

results$H142C <- tecan.30[c('O24', 'O23',
                            'O22', 'O21',
                            'O20', 'O19',
                            'O18', 'O17',
                            'O16', 'O15',
                            'O14', 'O13',
                            'O12', 'O11',
                            'O10', 'O9',
                            'O8', 'O7',
                            'O6', 'O5',
                            'O4', 'O3',
                            'O2', 'O1'),]

results$L144A <- tecan.30[c('P24', 'P23',
                            'P22', 'P21',
                            'P20', 'P19',
                            'P18', 'P17',
                            'P16', 'P15',
                            'P14', 'P13',
                            'P12', 'P11',
                            'P10', 'P9',
                            'P8', 'P7',
                            'P6', 'P5',
                            'P4', 'P3',
                            'P2', 'P1'),]

results <- lapply(results, function(x) {colnames(x) <- as.numeric(tecan.30[1,])/3600; class(x) <- 'numeric'; return(x)})

## determine growth rates
PYL1.TECAN.curves.gr <- results
for (i in 1:length(PYL1.TECAN.curves.gr)){
  
  tmp.out <- vector(mode = 'list', length = nrow(PYL1.TECAN.curves.gr[[i]]))
  
  for (k in 1:length(tmp.out)){
    
    tmp.out[[k]] <- rbind(as.numeric(colnames(PYL1.TECAN.curves.gr[[i]])),
                          as.numeric(PYL1.TECAN.curves.gr[[i]][k,]))
    tmp.out[[k]] <- tmp.out[[k]][,!is.na(tmp.out[[k]][2,])]
    tmp.out[[k]] <- fit_easylinear(time = tmp.out[[k]][1,],
                                   y = tmp.out[[k]][2,],
                                   h = 15)
    tmp.out[[k]] <- tmp.out[[k]]@par[['mumax']]
  }
  
  PYL1.TECAN.curves.gr[[i]] <- do.call(c, tmp.out)
  
}

### fit WT curve Hill model
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}

WT.PYL1.drc <- cbind(PYL1.TECAN.curves.gr$WT,
                     rep(sort(dosages), each = 2))
class(WT.PYL1.drc) <- "numeric"
WT.PYL1.drc <- as.data.frame(WT.PYL1.drc)
colnames(WT.PYL1.drc) <- c("GR", "concentration")

WT.PYL1.drc <- drm(WT.PYL1.drc$GR ~ WT.PYL1.drc$concentration,
                   fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                   type = 'continuous')
WT.PYL1.drc.par <- WT.PYL1.drc$fit$par
names(WT.PYL1.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.par <- WT.PYL1.drc.par[c(2:4,1)]
WT.PYL1.drc.par[4] <- -WT.PYL1.drc.par[4]

## standardise TECAN data to the WT curve's B[inf]
PYL1.TECAN.curves.gr <- lapply(PYL1.TECAN.curves.gr, function(x){y <- 100*c(x/WT.PYL1.drc.par["B[inf]"]); return(y)})

## calculate dose-response curves using duplicate measurements, generate the confidence intervals at the key concentrations
PYL1.TECAN.curves.gr.predict <- vector(mode = "list", length = length(PYL1.TECAN.curves.gr))
names(PYL1.TECAN.curves.gr.predict) <- names(PYL1.TECAN.curves.gr)
for(i in 1:length(PYL1.TECAN.curves.gr)){

  print(names(PYL1.TECAN.curves.gr.predict)[i])
    
  ## temporary input data frame
  tmp.drc <- cbind(PYL1.TECAN.curves.gr[[i]],
                   rep(sort(dosages), each = 2))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")
  
  ## temporary Hill fit & parameters 
  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')

  ## temporary prediction
  tmp.drc.predict.newdata <- expand.grid(conc = sort(dosages))
  tmp.drc.predict <- predict(tmp.drc,
                             newdata = tmp.drc.predict.newdata,
                             interval = "confidence",
                             level = 0.95)
  rownames(tmp.drc.predict) <- sort(dosages)

  ## export
  PYL1.TECAN.curves.gr.predict[[i]] <- tmp.drc.predict
  
  ## clean up
  rm(tmp.drc, tmp.drc.predict.newdata, tmp.drc.predict)
  
}

## clean up
rm(results,setup,tecan.30,i,k,tmp.out,WT.PYL1.drc,WT.PYL1.drc.par)


## 2. Pre-processed DiMSum data ##
##################################

## load pre-processed GluePCA DMS data 
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## load dose-response curve metrics
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## which TECAN variants have high confidence DMS library fits?
parameters.Hill.lib <- parameters.Hill[match(names(PYL1.TECAN.curves.gr), rownames(parameters.Hill)),]
parameters.Hill.lib.hq <- parameters.Hill.lib[which(parameters.Hill.lib[,"R^2"] > 0.9),]

## filter the TECAN curvce fits accordingly
PYL1.TECAN.curves.gr <- PYL1.TECAN.curves.gr[match(rownames(parameters.Hill.lib.hq), names(PYL1.TECAN.curves.gr))]
PYL1.TECAN.curves.gr.predict <- PYL1.TECAN.curves.gr.predict[match(rownames(parameters.Hill.lib.hq), names(PYL1.TECAN.curves.gr.predict))]

## summarise the data
PYL1.ABI1.keyvars <- matrix(NA, ncol = 12, nrow = length(PYL1.TECAN.curves.gr))
rownames(PYL1.ABI1.keyvars) <- names(PYL1.TECAN.curves.gr)
colnames(PYL1.ABI1.keyvars) <- names(PYL1.ABI1)
for (i in 1:nrow(PYL1.ABI1.keyvars)){
  if(i == 1){
    PYL1.ABI1.keyvars[i,] <- sapply(PYL1.ABI1, function(x){x <- x[which(x$WT == T)[1],"gr_normalised_WTscaled"]; return(x)})
  }else{
    tmp.wt <- substr(names(PYL1.TECAN.curves.gr)[i], 1, 1)
    tmp.pos <- as.numeric(substr(names(PYL1.TECAN.curves.gr)[i], 2, nchar(names(PYL1.TECAN.curves.gr)[i]) - 1))
    tmp.mut <- substr(names(PYL1.TECAN.curves.gr)[i], nchar(names(PYL1.TECAN.curves.gr)[i]), nchar(names(PYL1.TECAN.curves.gr)[i]))
    PYL1.ABI1.keyvars[i,] <- sapply(PYL1.ABI1, function(x){x <- x[which(x$WT_AA == tmp.wt & x$Pos == tmp.pos & x$Mut == tmp.mut),"gr_normalised_WTscaled"]; return(x)})
  }
}
 
## remake the curves, generate the confidence intervals at the key concentrations
PYL1.ABI1.keyvars.out <- vector(mode = "list", length = nrow(PYL1.ABI1.keyvars))
names(PYL1.ABI1.keyvars.out) <- names(PYL1.TECAN.curves.gr)
for(i in 1:length(PYL1.ABI1.keyvars.out)){

  tmp.drc <- cbind(PYL1.ABI1.keyvars[i,],colnames(PYL1.ABI1.keyvars))
  class(tmp.drc) <- "numeric"
  tmp.drc <- as.data.frame(tmp.drc)
  colnames(tmp.drc) <- c("B", "concentration")

  tmp.drc <- drm(tmp.drc$B ~ tmp.drc$concentration,
                 fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                 type = 'continuous')

  ## predict the full, smoothened curve (using 1000 data points) and confidence interval
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
PYL1.mutant.summary <- vector(mode = "list", length = length(PYL1.TECAN.curves.gr))
names(PYL1.mutant.summary) <- names(PYL1.TECAN.curves.gr)
for(i in 1:length(PYL1.mutant.summary)){
  
  #### matrix setup
  PYL1.mutant.summary[[i]] <- matrix(NA, ncol = 9, nrow = 12)
  colnames(PYL1.mutant.summary[[i]]) <- c("TECAN1_gr",
                                          "TECAN2_gr",
                                          "TECAN_predict",
                                          "TECAN_predict_min",
                                          "TECAN_predict_max",
                                          "Competition_gr",
                                          "Competition_predict",
                                          "Competition_predict_min",
                                          "Competition_predict_max")
  rownames(PYL1.mutant.summary[[i]]) <- sort(dosages)
  
  #### fill in TECAN data
  PYL1.mutant.summary[[i]][,"TECAN1_gr"] <- PYL1.TECAN.curves.gr[[i]][seq(f=1,to=23,by=2)]
  PYL1.mutant.summary[[i]][,"TECAN2_gr"] <- PYL1.TECAN.curves.gr[[i]][seq(f=2,to=24,by=2)]
  PYL1.mutant.summary[[i]][,"TECAN_predict"] <- PYL1.TECAN.curves.gr.predict[[i]][,"Prediction"]
  PYL1.mutant.summary[[i]][,"TECAN_predict_min"] <- PYL1.TECAN.curves.gr.predict[[i]][,"Lower"]
  PYL1.mutant.summary[[i]][,"TECAN_predict_max"] <- PYL1.TECAN.curves.gr.predict[[i]][,"Upper"]

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
TECAN.mean <- do.call(c, lapply(PYL1.mutant.summary, function(x){apply(x[,1:2], 1, function(x){mean(x)})}))
TECAN.lower <- do.call(c, lapply(PYL1.mutant.summary, function(x){apply(x[,1:2], 1, function(x){min(x)})}))
TECAN.upper <- do.call(c, lapply(PYL1.mutant.summary, function(x){apply(x[,1:2], 1, function(x){max(x)})}))
TECAN.pred <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,3]}))
TECAN.pred.min <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,4]}))
TECAN.pred.max <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,5]}))

DMS <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,6]}))
DMS.pred <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,7]}))
DMS.pred.min <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,8]}))
DMS.pred.max <- do.call(c, lapply(PYL1.mutant.summary, function(x){x[,9]}))

out.tecan.df <- cbind("conc" = rep(1:12, 9),
                      "tecan mean" = TECAN.mean, 
                      "tecan min" = TECAN.lower, 
                      "tecan max" = TECAN.upper, 
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
# [1] rlang_1.1.6       ggtext_0.1.2      ggplot2_4.0.0     drc_3.0-1         MASS_7.3-65       reshape_0.8.10    growthrates_0.8.5
# [8] deSolve_1.40      lattice_0.22-7    readxl_1.4.5     
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       FME_1.3.6.4        magrittr_2.0.3    
# [7] grid_4.5.1         RColorBrewer_1.1-3 mvtnorm_1.3-3      cellranger_1.1.0   plyr_1.8.9         Matrix_1.7-4      
# [13] Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    mgcv_1.9-3         scales_1.4.0       TH.data_1.1-4     
# [19] codetools_0.2-20   abind_1.4-8        cli_3.6.5          crayon_1.5.3       splines_4.5.1      withr_3.0.2       
# [25] plotrix_3.8-4      rootSolve_1.8.2.4  tools_4.5.1        parallel_4.5.1     coda_0.19-4.1      minpack.lm_1.2-4  
# [31] minqa_1.2.8        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4   
# [37] car_3.1-3          pkgconfig_2.0.3    pillar_1.11.0      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0        
# [43] tibble_3.3.0       tidyselect_1.2.1   rstudioapi_0.17.1  farver_2.1.2       nlme_3.1-168       labeling_0.4.3    
# [49] carData_3.0-5      compiler_4.5.1     S7_0.2.0           gridtext_0.1.5
