# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################
## Figure 6A - curves with inverted response ##
###############################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "drc", "reshape", "ggplot2", "ggtext", "wesanderson", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

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


## 2. Generate curves ##
########################

## calculate curves
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")
input <- parameters.Hill[match(c("V110W", "A116W", "S119R", "H142W"),rownames(parameters.Hill)),]

## Define the four-parameter logistic function
four_pl_function <- function(x, Hill, B0, Binf, EC50) {
  B0 + (Binf - B0) / (1 + (EC50 / x)^c(Hill))
}

## Calculate the predicted response over the dose range
curves <- matrix(NA, nrow = nrow(input), ncol = 1000)
rownames(curves) <- rownames(input)
colnames(curves) <- c(0, exp(seq(log(0.001), log(5000), length = 999)))
curves.min <- curves.max <- curves

dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
for(i in 1:nrow(curves)){
  
  print(i)
  
  ## fetch raw data
  tmp.id <- rownames(curves)[i]
  PYL1.summary.nonWT[match(tmp.id, rownames(PYL1.summary.nonWT)),]
  
  ## estimate confidence interval
  tmp.df <- as.data.frame(cbind(conc = dosages,
                                out =  PYL1.summary.nonWT[match(tmp.id, rownames(PYL1.summary.nonWT)),]))
  
  if(i == 1){
    tmp.model <- drm(tmp.df$out ~ tmp.df$conc, 
                     fct = LL.4(fixed = c(NA, NA, NA, NA)),
                     type = 'continuous') 
    par.WT.PYL1.drc.Hill <- tmp.model$fit$par
    names(par.WT.PYL1.drc.Hill) <- c("Hill", "B[0]", "B[inf]", "EC50")
  }else{
    tmp.model <- drm(tmp.df$out ~ tmp.df$conc, 
                     fct = LL.4(fixed = c(NA, NA, NA, NA)),
                     type = 'continuous',
                     start = par.WT.PYL1.drc.Hill)
  }
  tmp.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  tmp.model.predict <- predict(tmp.model,
                               newdata = tmp.predict.newdata,
                               interval = "confidence")
  
  ## fetch curves and confidences interval
  curves[i,] <- tmp.model.predict[,"Prediction"]
  curves.min[i,] <- tmp.model.predict[,"Lower"]
  curves.max[i,] <- tmp.model.predict[,"Upper"]
  
}

## mean
curves <- rbind(curves,
                "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves <- t(curves)
curves <- as.data.frame(curves)
curves <- melt(curves, id.vars = "conc")

## max
curves.max <- rbind(curves.max,
                    "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves.max <- t(curves.max)
curves.max <- as.data.frame(curves.max)
curves.max <- melt(curves.max, id.vars = "conc")

## min
curves.min <- rbind(curves.min,
                    "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves.min <- t(curves.min)
curves.min <- as.data.frame(curves.min)
curves.min <- melt(curves.min, id.vars = "conc")

## clean up
rm(i, four_pl_function, tmp.df, tmp.model, tmp.model.predict, tmp.predict.newdata, tmp.id, par.WT.PYL1.drc.Hill)


## 4. Plot all inverter curves ##
#################################

### curves
inverter.plot <- function(id){
  
  ## extract curve fit and 95% confidence interval
  curve.tmp <- curves[grep(id, curves$variable),]
  curve.tmp$min.value <- curves.min[grep(id, curves.min$variable),"value"]
  curve.tmp$max.value <- curves.max[grep(id, curves.max$variable),"value"]
  curve.tmp$variable <- factor(curve.tmp$variable)
  
  ## extract data points
  points.tmp <- PYL1.summary.nonWT[grep(id, rownames(PYL1.summary.nonWT)),,drop = F]
  points.tmp <- as.data.frame(t(points.tmp))
  dosages[12] <- 9.062741e-03/3.5/3.5/3.5
  points.tmp$conc <- dosages
  colnames(points.tmp)[1] <- "value"
  
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
  
  p <- ggplot(data = curve.tmp, aes(x = conc, y = variable)) +
    geom_ribbon(data = curve.tmp,
                aes(group = variable, x = conc, y = value, ymin = min.value, ymax = max.value),
                alpha = 0.1, fill = "grey50") +
    geom_line(data = curve.tmp,
              aes(group = variable, x = conc, y = value), 
              linewidth = 5, alpha = 1, colour = cols) +
    geom_point(data = points.tmp, aes(x = conc, y = value),
               size = 13, alpha = 1, colour = cols, shape = 16) +
    annotate("text", x = 0.0002113759, y = 10,
             size = 25, label = expr_text(bquote(
               italic(R)^2 == .(format(round(parameters.Hill[grep(id, rownames(parameters.Hill)), "R^2"], 2), digits = 2)) *
                 ", " *
                 Delta * "Binding = " * .(format(round(parameters.Hill[grep(id, rownames(parameters.Hill)), "B[inf]"] -
                                                       parameters.Hill[grep(id, rownames(parameters.Hill)),"B[0]"], 1), nsmall = 1)) * "%"
             ),
             parse = T,
             color = "black", hjust = 0)) +
    scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                  limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-10000, 10000)) +
    coord_cartesian(ylim = c(-5, 120)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(size = 80),
          axis.text = element_text(size = 40),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 70, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(1, 1, 2, 2),"cm")) +
    labs(title = id,
         x = "(+)-ABA conc. (µM)",
         y = "Relative PYL1/ABI1 Binding")
print(p)
  
}

pdf("../../results/Figure6/Figure6A_V110W_curve.pdf", height = 13, width = 22)
inverter.plot(id = "V110W")
dev.off()

pdf("../../results/Figure6/Figure6A_A116W_curve.pdf", height = 13, width = 22)
inverter.plot(id = "A116W")
dev.off()

pdf("../../results/Figure6/Figure6A_S119R_curve.pdf", height = 13, width = 22)
inverter.plot(id = "S119R")
dev.off()

pdf("../../results/Figure6/Figure6A_H142W_curve.pdf", height = 13, width = 22)
inverter.plot(id = "H142W")
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
# [1] rlang_1.1.6       wesanderson_0.3.7 ggtext_0.1.2      ggplot2_4.0.0     reshape_0.8.10    drc_3.0-1        
# [7] MASS_7.3-65       scales_1.4.0      stringr_1.5.1    
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     xml2_1.4.0         gtools_3.9.5       stringi_1.8.7      lattice_0.22-7    
# [7] magrittr_2.0.3     grid_4.5.1         RColorBrewer_1.1-3 mvtnorm_1.3-3      plyr_1.8.9         Matrix_1.7-4      
# [13] Formula_1.2-5      survival_3.8-3     multcomp_1.4-28    TH.data_1.1-4      codetools_0.2-20   abind_1.4-8       
# [19] cli_3.6.5          litedown_0.7       commonmark_2.0.0   splines_4.5.1      withr_3.0.2        plotrix_3.8-4     
# [25] tools_4.5.1        dplyr_1.1.4        vctrs_0.6.5        R6_2.6.1           zoo_1.8-14         lifecycle_1.0.4   
# [31] car_3.1-3          pkgconfig_2.0.3    pillar_1.11.0      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0        
# [37] xfun_0.53          tibble_3.3.0       tidyselect_1.2.1   rstudioapi_0.17.1  farver_2.1.2       carData_3.0-5     
# [43] compiler_4.5.1     S7_0.2.0           markdown_2.0       gridtext_0.1.5   
