# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##############################################################
## Figure 3C - example curves with significantly lower EC50 ##
##############################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "drc", "reshape", "ggplot2", "ggtext")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

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


## 2. Generate curves ##
########################

## calculate curves
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
input <- parameters.Hill[grep("WT|H87|R164",rownames(parameters.Hill)),]

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


## 3. Plot H87 variants ##
##########################

### curves
curves.H87 <- curves[grep("H87P|H87V|H87A|H87G", curves$variable),]
curves.H87$min.value <- curves.min[grep("H87P|H87V|H87A|H87G", curves.min$variable),"value"]
curves.H87$max.value <- curves.max[grep("H87P|H87V|H87A|H87G", curves.max$variable),"value"]
curves.H87$variable <- factor(curves.H87$variable)

### data points
points.H87 <- PYL1.summary.nonWT[grep("H87P|H87V|H87A|H87G", rownames(PYL1.summary.nonWT)),]
points.H87 <- as.data.frame(t(points.H87))
dosages[12] <- 9.062741e-03/3.5/3.5/3.5
points.H87$conc <- dosages

## colours like in panel B
cols.all <- colorRampPalette(c("blue", "gray90", "red"))(1000)
names(cols.all) <- seq(f = -2, t = 2, length.out = 1000)
ec50_logFC <- log10(parameters.Hill[-1,"EC50"] / parameters.Hill["WT","EC50"])
cols.H87 <- ec50_logFC[grep("H87P|H87V|H87A|H87G", names(ec50_logFC))]
cols.H87 <- as.data.frame(cbind("logFC" = cols.H87))
cols.H87$col <- rep(NA, nrow(cols.H87))
cols.H87$col[1] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.H87$logFC[1]))]
cols.H87$col[2] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.H87$logFC[2]))]
cols.H87$col[3] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.H87$logFC[3]))]
cols.H87$col[4] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.H87$logFC[4]))]
curves.H87$col <- cols.H87[match(as.character(curves.H87$variable), rownames(cols.H87)),"col"]

pdf("../../results/Figure3/Figure3C_H87_curves.pdf", height = 13, width = 17)

out.3C.H87 <- ggplot(data = curves.H87, aes(x = conc, y = variable)) +
  geom_ribbon(data = curves.H87,
              aes(group = variable, x = conc, y = value, ymin = min.value, ymax = max.value),
              alpha = 0.1, fill = "grey50") +
  geom_line(data = curves.H87,
            aes(group = variable, x = conc, y = value), 
            linewidth = 3, alpha = 1, colour = curves.H87$col) +
  geom_point(data = points.H87, aes(x = conc, y = H87G),
             size = 6, shape = 16, alpha = 0.6, colour = cols.H87$col[1]) +
  geom_point(data = points.H87, aes(x = conc, y = H87A),
             size = 6, shape = 16, alpha = 0.6, colour = cols.H87$col[2]) +
  geom_point(data = points.H87, aes(x = conc, y = H87V),
             size = 6, shape = 16, alpha = 0.6, colour = cols.H87$col[3]) +
  geom_point(data = points.H87, aes(x = conc, y = H87P),
             size = 6, shape = 16, alpha = 0.6, colour = cols.H87$col[4]) +
  annotate("text", x = 0.0003, y = parameters.Hill["H87P","B[0]"] + 10, 
           size = 15, label = paste0("P (", format(round(parameters.Hill["H87P","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.H87["H87P","col"], hjust = 0) +
  annotate("text", x = 0.0003, y = parameters.Hill["H87V","B[0]"] + 10, 
           size = 15, label = paste0("V (", format(round(parameters.Hill["H87V","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.H87["H87V","col"], hjust = 0) +
  annotate("text", x = 0.0003, y = parameters.Hill["H87A","B[0]"] + 10, 
           size = 15, label = paste0("A (", format(round(parameters.Hill["H87A","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.H87["H87A","col"], hjust = 0) +
  annotate("text", x = 0.0003, y = parameters.Hill["H87G","B[0]"] + 10, 
           size = 15, label = paste0("G (", format(round(parameters.Hill["H87G","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.H87["H87G","col"], hjust = 0) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-150, 150)) +
  coord_cartesian(ylim = c(-5, 120)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(size = 80),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 55, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 55, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 2, 2),"cm")) +
  labs(title = "H87",
       x = "(+)-ABA concentration (µM)",
       y = "Relative PYL1/ABI1 Binding")

print(out.3C.H87)

dev.off()


## 4. Plot R164 variants ##
###########################

### curves
curves.R164 <- curves[grep("R164D|R164N|R164W", curves$variable),]
curves.R164$min.value <- curves.min[grep("R164D|R164N|R164W", curves.min$variable),"value"]
curves.R164$max.value <- curves.max[grep("R164D|R164N|R164W", curves.max$variable),"value"]
curves.R164$variable <- factor(curves.R164$variable)

### data points
points.R164 <- PYL1.summary.nonWT[grep("R164D|R164N|R164W", rownames(PYL1.summary.nonWT)),]
points.R164 <- as.data.frame(t(points.R164))
dosages[12] <- 9.062741e-03/3.5/3.5/3.5
points.R164$conc <- dosages

## colours like in panel B
cols.all <- colorRampPalette(c("blue", "gray90", "red"))(1000)
names(cols.all) <- seq(f = -2, t = 2, length.out = 1000)
ec50_logFC <- log10(parameters.Hill[-1,"EC50"] / parameters.Hill["WT","EC50"])
cols.R164 <- ec50_logFC[grep("R164D|R164N|R164W", names(ec50_logFC))]
cols.R164 <- as.data.frame(cbind("logFC" = cols.R164))
cols.R164$col <- rep(NA, length(cols.R164))
cols.R164$col[1] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.R164$logFC[1]))]
cols.R164$col[2] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.R164$logFC[2]))]
cols.R164$col[3] <- cols.all[which.min(abs(as.numeric(names(cols.all)) -  cols.R164$logFC[3]))]
curves.R164$col <- cols.R164[match(as.character(curves.R164$variable), rownames(cols.R164)),"col"]

pdf("../../results/Figure3/Figure3C_R164_curves.pdf", height = 13, width = 17)

out.3C.R164 <- ggplot(data = curves.R164, aes(x = conc, y = variable)) +
  geom_ribbon(data = curves.R164,
              aes(group = variable, x = conc, y = value, ymin = min.value, ymax = max.value),
              alpha = 0.1, fill = "grey50") +
  geom_line(data = curves.R164,
            aes(group = variable, x = conc, y = value), 
            linewidth = 3, alpha = 1, colour = curves.R164$col) +
  geom_point(data = points.R164, aes(x = conc, y = R164W),
             size = 7, shape = 16, alpha = 0.6, colour = cols.R164$col[1]) +
  geom_point(data = points.R164, aes(x = conc, y = R164D),
             size = 7, shape = 16, alpha = 0.6, colour = cols.R164$col[2]) +
  geom_point(data = points.R164, aes(x = conc, y = R164N),
             size = 7, shape = 16, alpha = 0.6, colour = cols.R164$col[3]) +
  annotate("text", x = 0.0003, y = parameters.Hill["R164W","B[0]"] + 10, 
           size = 15, label = paste0("W (3.80×)"), 
           color = cols.R164["R164W","col"], hjust = 0) +
  annotate("text", x = 0.0003, y = parameters.Hill["R164D","B[0]"] + 14, 
           size = 15, label = paste0("D (", format(round(parameters.Hill["R164D","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.R164["R164D","col"], hjust = 0) +
  annotate("text", x = 0.0003, y = parameters.Hill["R164N","B[0]"] + 16, 
           size = 15, label = paste0("N (", format(round(parameters.Hill["R164N","EC50"]/parameters.Hill["WT","EC50"], digits = 2), nsmall = 2), "×)"), 
           color = cols.R164["R164N","col"], hjust = 0) +
  scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-150, 150)) +
  coord_cartesian(ylim = c(-5, 120)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(size = 80),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(linewidth = 1, color = 'black'),
        axis.line.y = element_line(linewidth = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 55, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 55, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 2, 2),"cm")) +
  labs(title = "R164",
       x = "(+)-ABA concentration (µM)",
       y = "Relative PYL1/ABI1 Binding")

print(out.3C.R164)

dev.off()


## 5. Version ##
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
# [1] ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 drc_3.0-1      MASS_7.3-65    scales_1.4.0   stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       dplyr_1.1.4        compiler_4.5.1     gtools_3.9.5       tidyselect_1.2.1   plotrix_3.8-4      Rcpp_1.1.0         xml2_1.4.0         splines_4.5.1     
# [11] lattice_0.22-7     TH.data_1.1-4      R6_2.6.1           plyr_1.8.9         commonmark_2.0.0   generics_0.1.4     Formula_1.2-5      tibble_3.3.0       car_3.1-3          pillar_1.11.1     
# [21] RColorBrewer_1.1-3 rlang_1.1.6        multcomp_1.4-28    stringi_1.8.7      xfun_0.53          litedown_0.7       S7_0.2.0           cli_3.6.5          withr_3.0.2        magrittr_2.0.4    
# [31] gridtext_0.1.5     grid_4.5.1         mvtnorm_1.3-3      markdown_2.0       sandwich_3.1-1     lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2       codetools_0.2-20  
# [41] zoo_1.8-14         survival_3.8-3     abind_1.4-8        carData_3.0-5      pkgconfig_2.0.3    tools_4.5.1