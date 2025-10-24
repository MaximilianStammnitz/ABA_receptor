# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 3B - PYL1-ABI1 TECAN validations ##
###########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthrates", "drc", 
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


## 2. Import TECAN data ##
##########################

## raw data import
PYL1.ABI1.TECAN <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[20:405,1:263]
rownames(PYL1.ABI1.TECAN) <- PYL1.ABI1.TECAN[,1]
colnames(PYL1.ABI1.TECAN) <- PYL1.ABI1.TECAN[1,]
PYL1.ABI1.TECAN <- PYL1.ABI1.TECAN[,-1]
PYL1.ABI1.TECAN <- PYL1.ABI1.TECAN[-c(1:2),]
colnames(PYL1.ABI1.TECAN) <- as.numeric(colnames(PYL1.ABI1.TECAN))/3600 ## convert to hours
class(PYL1.ABI1.TECAN) <- "numeric"
PYL1.ABI1.TECAN <- as.data.frame(PYL1.ABI1.TECAN)

## assign mutants
setup <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[1:16,1:25]
rownames(setup) <- setup[,1]
setup <- setup[,-1]
colnames(setup) <- 1:24
mut.names <- unique(c(setup[,1],setup[,13]))
mut.names <- mut.names[!is.na(mut.names)]
PYL1.ABI1.TECAN.muts <- vector(mode = "list", length = length(mut.names))
names(PYL1.ABI1.TECAN.muts) <- mut.names
for(i in 1:length(PYL1.ABI1.TECAN.muts)){
  PYL1.ABI1.TECAN.muts[[i]] <- which(setup == mut.names[i], arr.ind = TRUE)
  PYL1.ABI1.TECAN.muts[[i]] <- paste0(rownames(PYL1.ABI1.TECAN.muts[[i]]), PYL1.ABI1.TECAN.muts[[i]][,2])
  PYL1.ABI1.TECAN.muts[[i]] <- PYL1.ABI1.TECAN[PYL1.ABI1.TECAN.muts[[i]],]
}

## calculate growth rates
PYL1.ABI1.TECAN.muts.rates <- PYL1.ABI1.TECAN.muts
for (i in 1:length(PYL1.ABI1.TECAN.muts.rates)){
  
  print(i)
  
  tmp.out <- vector(mode = 'list', length = nrow(PYL1.ABI1.TECAN.muts.rates[[i]]))
  
  for (k in 1:length(tmp.out)){
    
    tmp.out[[k]] <- rbind(as.numeric(colnames(PYL1.ABI1.TECAN.muts.rates[[i]])),
                          as.numeric(PYL1.ABI1.TECAN.muts.rates[[i]][k,]))
    tmp.out[[k]] <- tmp.out[[k]][,!is.na(tmp.out[[k]][2,])]
    tmp.out[[k]] <- fit_easylinear(time = tmp.out[[k]][1,],
                                   y = tmp.out[[k]][2,],
                                   h = 15)
    tmp.out[[k]] <- tmp.out[[k]]@par[['mumax']]
    
  }
  
  PYL1.ABI1.TECAN.muts.rates[[i]] <- do.call(c, tmp.out)
  
}

## Dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}

## fit WT curve(s)
WT.PYL1.drc <- cbind(PYL1.ABI1.TECAN.muts.rates$WT[1:12],rev(dosages))
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

## rescale all curves to WT Bmax
PYL1.ABI1.TECAN.muts.rates <- lapply(PYL1.ABI1.TECAN.muts.rates, function(x){y <- 100*x/WT.PYL1.drc.par["B[inf]"]; return(y)})

## clean-up
rm(i,k,tmp.out)


## 3. Fit all TECAN curves' dose response profiles ##
#####################################################

parameters.Hill.TECAN <- matrix(NA, nrow = length(PYL1.ABI1.TECAN.muts.rates), ncol = 16)
colnames(parameters.Hill.TECAN) <- c("Hill", "B[0]", "B[inf]", "EC50", 
                                     "Hill SE", "B[0] SE", "B[inf] SE", "EC50 SE", 
                                     "Hill P", "B[0] P", "B[inf] P", "EC50 P", 
                                     "Data points", "AIC", "Residual var", "R^2")
rownames(parameters.Hill.TECAN) <- names(PYL1.ABI1.TECAN.muts.rates)

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
  tmp.PYL1.drc.in <- cbind(PYL1.ABI1.TECAN.muts.rates[[i]][1:12],rev(dosages))
  class(tmp.PYL1.drc.in) <- "numeric"
  tmp.PYL1.drc.in <- as.data.frame(tmp.PYL1.drc.in)
  colnames(tmp.PYL1.drc.in) <- c("GR", "concentration")
    
  #if(!rownames(parameters.Hill.TECAN)[i] %in% c("V110H", "A116H", "V193H", "V193W", "A190E", "H87P")){
  if(!rownames(parameters.Hill.TECAN)[i] %in% c("V110H", "A116H", "V193H", "V193W", "A116R", "A190E", "L125I")){
  #if(rownames(parameters.Hill.TECAN)[i] %in% c("V110Y", "V110W", "A116W", "H142W")){
    
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
#dev.off()


## 4. Plot ##
#############

## summarise data
out.all <- matrix(NA, nrow = 22, ncol = 26)
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
rownames(out.all) <- c("Q34Y", "Q34I", "E36R", "T38K", "Q39F", 
                       "H48C", "D80M", "P82G", "Y85P", "H87A",
                       "H87P", "A116R", "T118I", "L125I", "V132P", 
                       "I137G", "L144A", "E161A", "E171K", "A190E",
                       "R195M", "A202T")

for(i in 1:nrow(out.all)){
  tmp.name <- rownames(out.all)[i]
  tmp.id1 <- match(tmp.name, rownames(parameters.Hill))
  out1 <- parameters.Hill[tmp.id1,c(2,6,10,3,7,11,4,8,12,1,5,9,16)]
  tmp.id2 <- match(tmp.name, rownames(parameters.Hill.TECAN))
  out2 <- parameters.Hill.TECAN[tmp.id2,c(2,6,10,3,7,11,4,8,12,1,5,9,16)]
  out.all[i,] <- c(out1, out2)
}
out.all <- as.data.frame(out.all)
out.all <- out.all[which(out.all$bulk_R2 > 0.9 & out.all$TECAN_R2 > 0.9),]

## raw linear regression
p.B0 <- summary(lm(out.all$bulk_B0 ~ out.all$TECAN_B0))$coefficients[2,4]
p.Binf <- summary(lm(out.all$bulk_Binf ~ out.all$TECAN_Binf))$coefficients[2,4]
p.EC50 <- summary(lm(log10(out.all$bulk_EC50) ~ log10(out.all$TECAN_EC50)))$coefficients[2,4]
p.Hill <- summary(lm(log10(out.all$bulk_Hill) ~ log10(out.all$TECAN_Hill)))$coefficients[2,4]

## B0
r.B0 <- cor(x = out.all$bulk_B0, y = out.all$TECAN_B0, 
            method = "pearson", use = "complete.obs")

out.S3B_B0 <- ggplot(out.all, aes(x = `TECAN_B0`, y = `bulk_B0`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = pmax(`TECAN_B0` - `TECAN_B0_SE`, 0), 
                    xmax = `TECAN_B0` + `TECAN_B0_SE`), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = pmax(`bulk_B0` - `bulk_B0_SE`, 0),
                    ymax = `bulk_B0` + `bulk_B0_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.B0, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.B0, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[0] ~ "(individual mutants)"),
       y = bquote(B[0] ~ "(library sequencing)"))

## Binf
r.Binf <- cor(x = out.all$bulk_Binf, y = out.all$TECAN_Binf, 
              method = "pearson", use = "complete.obs")

out.S3B_Binf <- ggplot(out.all, aes(x = `TECAN_Binf`, y = `bulk_Binf`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = pmax(`TECAN_Binf` - `TECAN_Binf_SE`, 0),
                    xmax = `TECAN_Binf` + `TECAN_Binf_SE`), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = pmax(`bulk_Binf` - `bulk_Binf_SE`, 0),
                    ymax = `bulk_Binf` + `bulk_Binf_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0,
           y = 120,
           label = expr_text(bquote(italic(r) == .(format(r.Binf, digits = 2, nsmall = 2)) ~ 
                                      ", " ~ italic(P) == .(format(p.Binf, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(B[infinity] ~ "(individual mutants)"),
       y = bquote(B[infinity] ~ "(library sequencing)"))

## EC50
r.EC50 <- cor(x = log(out.all$bulk_EC50), y = log(out.all$TECAN_EC50), 
              method = "pearson", use = "complete.obs")

out.S3B_EC50 <- ggplot(out.all, aes(x = `TECAN_EC50`, y = `bulk_EC50`)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.001, 10000), ylim = c(0.001, 10000), expand = T) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = pmax(TECAN_EC50 - TECAN_EC50_SE, 0.001),
                    xmax = TECAN_EC50 + TECAN_EC50_SE), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = pmax(`bulk_EC50` - `bulk_EC50_SE`, 0.001),
                    ymax = `bulk_EC50` + `bulk_EC50_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
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
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(EC[50] ~ "(individual mutants)"),
       y = bquote(EC[50] ~ "(library sequencing)"))

## n
r.Hill <- cor(x = log(out.all$bulk_Hill), y = log(out.all$TECAN_Hill), 
              method = "pearson", use = "complete.obs")

out.S3B_Hill <- ggplot(out.all, aes(x = `TECAN_Hill`, y = `bulk_Hill`)) +
  scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5), 
                labels = c(0.1, 0.2, 0.5, 1, 2, 5),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5), 
                labels = c(0.1, 0.2, 0.5, 1, 2, 5),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.1, 10), ylim = c(0.1, 10), expand = T) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_Hill`, y = `bulk_Hill`),
             color = "black", size = 5, shape = 16) +
  geom_errorbar(aes(xmin = pmax(TECAN_Hill - TECAN_Hill_SE, 0.1),
                    xmax = TECAN_Hill + TECAN_Hill_SE), 
                color = "black", linewidth = 0.75, height = 0) +
  geom_errorbar(aes(ymin = pmax(`bulk_Hill` - `bulk_Hill_SE`, 0.1),
                    ymax = `bulk_Hill` + `bulk_Hill_SE`), 
                color = "black", linewidth = 0.75, width = 0) +
  annotate("text",
           x = 0.1186,
           y = 7.1,
           label = expr_text(bquote(italic(r) == .(format(r.Hill, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.Hill, digits = 4, nsmall = 4)))),
           parse = T,
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(n ~ "(individual mutants)"),
       y = bquote(n ~ "(library sequencing)"))

## combined plot
pdf("../../results/FigureS3/FigureS3B_Hill_parameter_validation.pdf", height = 15, width = 55)
plot_grid(out.S3B_B0, out.S3B_Binf, out.S3B_EC50, out.S3B_Hill, align = "hv", axis = "tblr", ncol = 4)
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
# [1] rlang_1.1.6       cowplot_1.2.0     ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    scales_1.4.0      drc_3.0-1         MASS_7.3-65       growthrates_0.8.5 deSolve_1.40      lattice_0.22-7    readxl_1.4.5     
# [13] stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      FME_1.3.6.4        magrittr_2.0.4     grid_4.5.1         RColorBrewer_1.1-3 mvtnorm_1.3-3      cellranger_1.1.0   plyr_1.8.9        
# [13] Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    mgcv_1.9-3         TH.data_1.1-4      codetools_0.2-20   abind_1.4-8        cli_3.6.5          crayon_1.5.3       splines_4.5.1      withr_3.0.2       
# [25] plotrix_3.8-4      rootSolve_1.8.2.4  tools_4.5.1        parallel_4.5.1     coda_0.19-4.1      minpack.lm_1.2-4   minqa_1.2.8        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4   
# [37] car_3.1-3          pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         tidyselect_1.2.1   tibble_3.3.0       farver_2.1.2       nlme_3.1-168       carData_3.0-5      compiler_4.5.1    
# [49] S7_0.2.0           gridtext_0.1.5 
