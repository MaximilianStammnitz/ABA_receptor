# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################
## Figure 6B - curves with bandstop response ##
###############################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "drc", "reshape", "ggplot2", "ggtext", 
              "wesanderson", "rlang", "tgp")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Table summary (21 (+)-ABA concentrations) ##
##################################################

load("../../data/DiMSum/PYL1-ABI1_21conc/PYL1-ABI1_preprocessed_21conc.RData")
PYL1.ABI1 <- PYL1_block
rm(PYL1_block)

### remove synonymous mutations, except for the main wildtype
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[-which(x$Nham_aa == 0 & x$Nham_nt > 0),]; return(x)})

## Summarise WT-rescaled binding fitness for each variant vs. ABA combination
## Matrix format
PYL1.ABI1.summary <- matrix(NA, nrow = 32*21, ncol = 21)
colnames(PYL1.ABI1.summary) <- names(PYL1.ABI1)
ids <- c()
for(i in 1:32){
  tmp.pos <- rep(str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"], "", 32)[i], 21)
  ids <- c(ids, paste0(tmp.pos, i + 117))
}
rownames(PYL1.ABI1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
                                                 "I", "F", "Y", "W", "K",
                                                 "R", "H", "D", "E", "S",
                                                 "T", "C", "N", "Q", "P",
                                                 "*"),32))

## Fill
WT <- str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"],"",32)[1,]
for(j in 1:length(PYL1.ABI1)){
  
  tmp.vars <- str_split_fixed(PYL1.ABI1[[j]]$aa_seq,"",32)
  
  for(i in 1:nrow(tmp.vars)){
    
    ### only take the "best" WT
    if(all(tmp.vars[i,] == WT)){
      
      PYL1.ABI1.summary[paste0(WT, 118:149, WT),j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
      next
      
    }
    
    pos.mut <- which(tmp.vars[i,] != WT)
    
    ### skip if higher-order mutant
    if(length(pos.mut) > 1){
      
      next
      
    }
    
    ### otherwise obtain the corresponding fitness
    tmp.mut <- paste0(WT[pos.mut], pos.mut + 117, tmp.vars[i,pos.mut])
    PYL1.ABI1.summary[tmp.mut,j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
    
  }
  
}
rm(i,j,tmp.vars,ids,pos.mut,tmp.mut,tmp.pos)

## Subset genotypes to only occur once (i.e. only one WT)
WTs.pos <- match(paste0(WT, 118:149, WT), rownames(PYL1.ABI1.summary))
PYL1.ABI1.summary.nonWT <- PYL1.ABI1.summary[-WTs.pos[-1],]
PYL1.ABI1.summary.nonWT <- PYL1.ABI1.summary.nonWT[c(WTs.pos[1],1:c(WTs.pos[1]-1),c(WTs.pos[1]+1):nrow(PYL1.ABI1.summary.nonWT)),]
rownames(PYL1.ABI1.summary.nonWT)[1] <- "WT"


## 2. Generate curves ##
########################

## calculate curves
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_21conc.RData")
input <- parameters.Hill[match(c("S119V", "H142S", "R145P"),rownames(parameters.Hill)),]

## dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[21] <- 0
for(i in 2:20){
  dosages[i] <- 2500 * (1/2.5)^c(i-1)
}


## 3. Plot all inverter curves for the 21-conc. DMS series ##
#############################################################

### curves
bandstop.plot.21conc <- function(id){
  
  ## extract curve fit and 95% confidence interval
  # curve.tmp <- curves[grep(id, curves$variable),]
  # curve.tmp$min.value <- curves.min[grep(id, curves.min$variable),"value"]
  # curve.tmp$max.value <- curves.max[grep(id, curves.max$variable),"value"]
  # curve.tmp$variable <- factor(curve.tmp$variable)
  
  ## extract data points
  points.tmp <- PYL1.ABI1.summary.nonWT[grep(id, rownames(PYL1.ABI1.summary.nonWT)),,drop = F]
  points.tmp <- as.data.frame(t(points.tmp))
  dosages[21] <- 6.87194767360001e-05/2.5/2.5/2.5 ## "0-conc." positioning for log scale
  points.tmp$conc <- log(dosages)
  colnames(points.tmp)[1] <- "value"
  
  ## run bayesian gaussian process modelling (4 MCMC chains)
  bgp.model.tmp <- bgp(Z = points.tmp$value,
                       X = points.tmp$conc,
                       BTE = c(5000, 10000, 10),
                       R = 4)
  
  #### predict a full, smoothened curve (using 1000 data points) and some sort of confidence interval
  tmp.newdata <- expand.grid(conc = c(0, seq(from = log(0.000001), to = log(2500), length.out = 999)))
  tmp.pm <- predict(object = bgp.model.tmp, 
                    XX = tmp.newdata, 
                    MAP = T,
                    zcov = T)
  
  #### add 95% confidence interval from model
  tmp.newdata$p <- tmp.pm$ZZ.km
  tmp.newdata$pmin <- tmp.pm$ZZ.km - 1.96*sqrt(tmp.pm$ZZ.ks2)
  tmp.newdata$pmax <- tmp.pm$ZZ.km + 1.96*sqrt(tmp.pm$ZZ.ks2)
  
  #### re-convert
  tmp.newdata$conc <- exp(tmp.newdata$conc)
  points.tmp$conc <- exp(points.tmp$conc)
  
  #### calculate delta binding
  detal.bind <- mean(points.tmp$value[c(1,21)]) - min(points.tmp$value)
  start.bind <- points.tmp$value[21] - min(points.tmp$value)
  end.bind <- points.tmp$value[1] - min(points.tmp$value)
  
  ## plot
  cols <- substr(id, nchar(id), nchar(id))
  if(cols %in% c("F", "Y", "W")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[5]
  }else if(cols %in% c("A", "V", "L", "I", "M")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[4]
  }else if(cols %in% c("Q", "S", "T", "N", "C")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[3]
  }else if(cols %in% c("K", "R", "H")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[2]
  }else if(cols %in% c("D", "E")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[1]
  }else if(cols %in% c("P", "G")){
    cols <- "black"
  }
  
  p <- ggplot(data = points.tmp, aes(x = conc, y = value)) +
    geom_ribbon(data = tmp.newdata,
                aes(x = conc, y = p, ymin = pmin, ymax = pmax),
                alpha = 0.1, fill = "grey50") +
    geom_line(data = tmp.newdata,
              aes(x = conc, y = p),
              linewidth = 4, alpha = 1, colour = cols) +
    geom_point(data = points.tmp, aes(x = conc, y = value),
               size = 13, alpha = 1, colour = cols, shape = 16) +
    annotate("text",
             x = 6.87194767360001e-05/2.5/2.5/2.5, y = 0,
             size = 25, label = expr_text(bquote(
               Delta * "Binding = " *
                 .(format(round(start.bind, 1), nsmall = 1)) *
                 " - " * .(format(round(end.bind, 1), nsmall = 1)) * "%")),
             parse = T,
             color = "black", hjust = 0) +
    scale_x_log10(breaks = c(6.87194767360001e-05/2.5/2.5/2.5, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, "1,000"),
                  limits = c(6.87194767360001e-05/2.5/2.5/2.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-10000, 10000)) +
    coord_cartesian(ylim = c(-5, 100)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(size = 80),
          axis.text = element_text(size = 40),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = -2),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(1, 1, 2, 2),"cm")) +
    labs(title = id,
         x = "(+)-ABA concentration (µM)",
         y = "Relative PYL1/ABI1 Binding")
  print(p)
  
}

pdf("../../results/Figure6/Figure6B_S119V_curve_21conc.pdf", height = 13, width = 28)
bandstop.plot.21conc(id = "S119V")
dev.off()

pdf("../../results/Figure6/Figure6B_H142S_curve_21conc.pdf", height = 13, width = 28)
bandstop.plot.21conc(id = "H142S")
dev.off()

pdf("../../results/Figure6/Figure6B_R145P_curve_21conc.pdf", height = 13, width = 28)
bandstop.plot.21conc(id = "R145P")
dev.off()

## clean up
rm(list=ls())


## 4. Table summary (21 (+)-ABA concentrations) ##
##################################################

## Load
load("../../data/DiMSum/PYL1-ABI1_12conc/PYL1-ABI1_preprocessed_12conc.RData")

## remove synonymous variants
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(x[,"Nham_aa"] == 0 & is.na(x[,"WT"]) == T),]; return(y)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(PYL1.ABI1[[1]][,"Nham_aa"] == 0 & is.na(PYL1.ABI1[[1]][,"WT"]) == F)[-1],]; return(y)})

## summarise in one matrix
PYL1.summary <- matrix(NA, nrow = 177*21, ncol = 12)
colnames(PYL1.summary) <- names(PYL1.ABI1)
ids <- c()
for(i in 1:177){
  tmp.pos <- rep(str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"], "", 177)[i], 21)
  ids <- c(ids, paste0(tmp.pos, i + 32))
}
rownames(PYL1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
                                            "I", "F", "Y", "W", "K",
                                            "R", "H", "D", "E", "S",
                                            "T", "C", "N", "Q", "P",
                                            "*"),177))
WT <- str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"],"",177)[1,]
for(j in 1:length(PYL1.ABI1)){

  tmp.vars <- str_split_fixed(PYL1.ABI1[[j]]$aa_seq,"",177)

  for(i in 1:nrow(tmp.vars)){

    ### only take the "best" WT
    if(all(tmp.vars[i,] == WT)){

      PYL1.summary[paste0(WT, 33:209, WT),j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
      next

    }

    pos.mut <- which(tmp.vars[i,] != WT)

    ### skip if higher-order mutant
    if(length(pos.mut) > 1){

      next

    }

    ### otherwise obtain the corresponding fitness
    tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
    PYL1.summary[tmp.mut,j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]

  }

}

## remove WT repetition
WTs.pos <- match(paste0(WT, 33:209, WT), rownames(PYL1.summary))
PYL1.summary.nonWT <- PYL1.summary[-WTs.pos[-1],]
PYL1.summary.nonWT <- PYL1.summary.nonWT[c(16,1:15,17:nrow(PYL1.summary.nonWT)),]
rownames(PYL1.summary.nonWT)[1] <- "WT"
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)

## clean up
rm(PYL1.ABI1, PYL1.summary, tmp.vars, i, ids, j, pos.mut, tmp.mut, tmp.pos, WT, WTs.pos, packages, install_if_missing)


## 5. Generate curves ##
########################

## calculate curves
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
input <- parameters.Hill[match(c("S119V", "H142S", "R145P"),rownames(parameters.Hill)),]

dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}


## 6. Plot all inverter curves ##
#################################

### curves
bandstop.plot.12conc <- function(id){

  ## extract data points
  points.tmp <- PYL1.summary.nonWT[grep(id, rownames(PYL1.summary.nonWT)),,drop = F]
  points.tmp <- as.data.frame(t(points.tmp))
  dosages[12] <- 6.87194767360001e-05/2.5/2.5/2.5 ## "0-conc." positioning for log scale
  points.tmp$conc <- log(dosages)
  colnames(points.tmp)[1] <- "value"
  
  ## run bayesian gaussian process modelling (4 MCMC chains)
  bgp.model.tmp <- bgp(Z = points.tmp$value,
                       X = points.tmp$conc,
                       BTE = c(5000, 10000, 10),
                       R = 4)
  
  #### predict a full, smoothened curve (using 1000 data points) and some sort of confidence interval
  tmp.newdata <- expand.grid(conc = c(0, seq(from = log(0.000001), to = log(2500), length.out = 999)))
  tmp.pm <- predict(object = bgp.model.tmp, 
                    XX = tmp.newdata, 
                    MAP = T,
                    zcov = T)
  
  #### add 95% confidence interval from model
  tmp.newdata$p <- tmp.pm$ZZ.km
  tmp.newdata$pmin <- tmp.pm$ZZ.km - 1.96*sqrt(tmp.pm$ZZ.ks2)
  tmp.newdata$pmax <- tmp.pm$ZZ.km + 1.96*sqrt(tmp.pm$ZZ.ks2)
  
  #### re-convert
  tmp.newdata$conc <- exp(tmp.newdata$conc)
  points.tmp$conc <- exp(points.tmp$conc)
  
  #### calculate delta binding
  delta.bind <- mean(points.tmp$value[c(1,12)]) - min(points.tmp$value)
  start.bind <- points.tmp$value[12] - min(points.tmp$value)
  end.bind <- points.tmp$value[1] - min(points.tmp$value)
  
  ## plot
  cols <- substr(id, nchar(id), nchar(id))
  if(cols %in% c("F", "Y", "W")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[5]
  }else if(cols %in% c("A", "V", "L", "I", "M")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[4]
  }else if(cols %in% c("Q", "S", "T", "N", "C")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[3]
  }else if(cols %in% c("K", "R", "H")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[2]
  }else if(cols %in% c("D", "E")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[1]
  }else if(cols %in% c("P", "G")){
    cols <- "black"
  }

  p <- ggplot(data = points.tmp, aes(x = conc, y = value)) +
    geom_line(data = tmp.newdata,
              aes(x = conc, y = p),
              linewidth = 4, alpha = 1, colour = cols) +
    geom_point(data = points.tmp, aes(x = conc, y = value),
               size = 8, alpha = 1, colour = cols, shape = 16) +
    annotate("text",
             x = 6.87194767360001e-05/2.5/2.5/2.5, 
             y = 0,
             size = 15, label = expr_text(bquote(
               Delta * "Binding = " *
                 .(format(round(start.bind, 1), nsmall = 1)) *
                 " - " * .(format(round(end.bind, 1), nsmall = 1)) * "%")),
             parse = T,
             color = "black", hjust = 0) +
    scale_x_log10(breaks = c(6.87194767360001e-05/2.5/2.5/2.5, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c("", "", "", "", "", "", "", "", ""),
                  limits = c(6.87194767360001e-05/2.5/2.5/2.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-10000, 10000)) +
    coord_cartesian(ylim = c(-5, 100)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(size = 80),
          axis.text = element_text(size = 40),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(1, 1, 2, 2),"cm"))
  print(p)

}

pdf("../../results/Figure6/Figure6B_S119V_curve_12conc.pdf", height = 6, width = 10)
bandstop.plot.12conc(id = "S119V")
dev.off()

pdf("../../results/Figure6/Figure6B_H142S_curve_12conc.pdf", height = 6, width = 10)
bandstop.plot.12conc(id = "H142S")
dev.off()

pdf("../../results/Figure6/Figure6B_R145P_curve_12conc.pdf", height = 6, width = 10)
bandstop.plot.12conc(id = "R145P")
dev.off()


## 7. Version ##
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
# [1] tgp_2.4-23        rlang_1.1.6       wesanderson_0.3.7 ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    drc_3.0-1         MASS_7.3-65       scales_1.4.0      stringr_1.5.2    
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      lattice_0.22-7     magrittr_2.0.4     grid_4.5.1         RColorBrewer_1.1-3 maptree_1.4-9     
# [11] mvtnorm_1.3-3      plyr_1.8.9         Matrix_1.7-4       Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    TH.data_1.1-4      codetools_0.2-20   abind_1.4-8        cli_3.6.5         
# [21] litedown_0.7       commonmark_2.0.0   splines_4.5.1      withr_3.0.2        plotrix_3.8-4      tools_4.5.1        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           rpart_4.1.24      
# [31] zoo_1.8-14         lifecycle_1.0.4    car_3.1-3          cluster_2.1.8.1    pkgconfig_2.0.3    pillar_1.11.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         xfun_0.53         
# [41] tibble_3.3.0       tidyselect_1.2.1   farver_2.1.2       carData_3.0-5      compiler_4.5.1     S7_0.2.0           markdown_2.0       gridtext_0.1.5   