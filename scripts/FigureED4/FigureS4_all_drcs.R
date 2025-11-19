# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#######################################################
## Supplementary Figure 4 - All dose response curves ##
#######################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "drc", "scales", "reshape", "ggplot2", "ggtext", "cowplot", "wesanderson")

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


## 2. Input and filter dose response curve parameters ##
########################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## unlike the high quality set in Figure 2D, here we only remove those without any fit
parameters.Hill <- parameters.Hill[-which(is.na(parameters.Hill[,"Hill"]) == T),]

## calculate curves

### Define the four-parameter logistic function
four_pl_function <- function(x, Hill, B0, Binf, EC50) {
  B0 + (Binf - B0) / (1 + (EC50 / x)^c(Hill))
}

### Calculate the predicted response over the dose range
curves <- matrix(NA, nrow = nrow(parameters.Hill), ncol = 1000)
rownames(curves) <- rownames(parameters.Hill)
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
rm(i, four_pl_function, tmp.df, tmp.model, tmp.model.predict, tmp.predict.newdata, dosages, tmp.id, par.WT.PYL1.drc.Hill)


## 3. Plot ##
#############

## PYL1 residues = columns

#### WT residues
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))

#### mutant order like in DMS heatmap(s)
mut.order <- c("G", "A", "V", "L", "M",
               "I", "F", "Y", "W", "K",
               "R", "H", "D", "E", "S",
               "T", "C", "N", "Q", "P",
               "*")

#### summarise dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
dosages[12] <- 9.062741e-03/3.5/3.5/3.5

#### mutant colours like in Taraneh's paper
cols <- c("G" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "P" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "*" = wes_palette("Darjeeling2", 6, type = "continuous")[6], 
          "W" = wes_palette("Darjeeling2", 6, type = "continuous")[5], "Y" = wes_palette("Darjeeling2", 6, type = "continuous")[5], 
          "F" = wes_palette("Darjeeling2", 6, type = "continuous")[5],
          "A" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "V" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "L" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
          "I" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "M" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
          "S" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "T" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "C" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
          "Q" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "N" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
          "K" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "R" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "H" = wes_palette("Darjeeling2", 6, type = "continuous")[2],
          "D" = wes_palette("Darjeeling2", 6, type = "continuous")[1], "E" = wes_palette("Darjeeling2", 6, type = "continuous")[1])
cols <- cols[match(mut.order, names(cols))]

#### collect plots (takes a few minutes to loop)
plots.out <- vector(mode = "list", length = 177)
names(plots.out) <- paste0(WT.seq, 33:209)
for(i in 1:length(WT.seq)){
  
  print(i)
  
  ### define summary object per PYL1 residue
  tmp.curves <- vector(mode = "list", length = 21)
  names(tmp.curves) <- paste0(WT.seq[i], i + 32, mut.order)
  
  ### find & rename "WT"
  tmp.WT <- which(names(tmp.curves) == paste0(WT.seq[i], i + 32, WT.seq[i]))
  names(tmp.curves)[tmp.WT] <- "WT"
  tmp.plot <- tmp.points <- tmp.curves
  
  #if(i == 1){
    
  for(j in 1:21){
      
      ### extract curves and raw data points [note: some may not be available (!?])
      tmp.curves[[j]] <- curves[which(curves$variable %in% names(tmp.curves)[j]),]
      tmp.points[[j]] <- as.data.frame(cbind("conc" = rev(dosages),
                                             "value" = rev(PYL1.summary.nonWT[which(rownames(PYL1.summary.nonWT) %in% names(tmp.points)[j]),])))
      
      #### confidence intervals
      tmp.curves[[j]]$min.value <- curves.min[which(curves.min$variable %in% names(tmp.curves)[j]),"value"]
      tmp.curves[[j]]$max.value <- curves.max[which(curves.max$variable %in% names(tmp.curves)[j]),"value"]
      
      ### save the plot(s)
      if(j != tmp.WT){
        
        tmp.plot[[j]] <- ggplot(data = tmp.curves[[j]]) +
          # geom_ribbon(data = tmp.curves[[j]],
          #             aes(x = conc, y = value, ymin = min.value, ymax = max.value),
          #             alpha = 0.1, fill = "grey50") +
          geom_line(data = tmp.curves[[j]], aes(x = conc, y = value), 
                    color = cols[j], linewidth = 2) +
          geom_point(data = tmp.points[[j]], aes(x = conc, y = value),
                     color = cols[j], size = 10) +
          scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                        labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                        limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
          scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-10000, 10000)) +
          coord_cartesian(ylim = c(-5, 115)) +
          theme_classic(base_size = 50) +
          theme(plot.title = element_markdown(),
                plot.subtitle = element_markdown(size = 20),
                title = element_text(size = 40),
                axis.line.x = element_line(size = 1, color = 'black'),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.line.y = element_line(size = 1, color = 'black'),
                axis.title.y = element_text(family = 'Helvetica', colour = cols[j], size = 100, angle = 0, vjust = 0.5, hjust = 0.5),
                axis.text.y = element_text(size = 30),
                legend.position = "none",
                text = element_text(family="Helvetica"),
                plot.margin = unit(c(1, 1, 0.5, 0.5),"cm")) +
          labs(x = NULL,
               y = NULL)
        
      }else{
        
        tmp.plot[[j]] <- ggplot(data = tmp.curves[[j]]) +
          geom_point(data = tmp.points[[j]], aes(x = conc, y = value),
                     color = "white", size = 10) +
          geom_line(data = tmp.curves[[j]], aes(x = conc, y = value), 
                    color = "white", linewidth = 2) +
          scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                        labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                        limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
          scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
          coord_cartesian(ylim = c(-5, 115)) +
          theme_classic(base_size = 50) +
          theme(plot.title = element_markdown(),
                plot.subtitle = element_markdown(size = 20),
                title = element_text(size = 40),
                axis.line.x = element_line(size = 1, color = 'black'),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.line.y = element_line(size = 1, color = 'black'),
                axis.title.y = element_text(family = 'Helvetica', colour = cols[j], size = 100, angle = 0, vjust = 0.5, hjust = 0.5),
                axis.text.y = element_text(size = 30),
                legend.position = "none",
                text = element_text(family="Helvetica"),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm")) +
          labs(x = NULL,
               y = NULL)
        
      }
   
    }   
  
  plots.out[[i]] <- tmp.plot
  
}

#### summarise the lists
out.1_30 <- do.call(c, plots.out[1:30])
out.31_60 <- do.call(c, plots.out[31:60])
out.61_90 <- do.call(c, plots.out[61:90])
out.91_120 <- do.call(c, plots.out[91:120])
out.121_150 <- do.call(c, plots.out[121:150])
out.151_180 <- do.call(c, plots.out[151:180])

pdf("../../results/FigureS4/FigureS4_DRCs_res033-062.pdf", width = 250, height = 100)
plot_grid(plotlist = out.1_30,
          ncol = 30,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()

pdf("../../results/FigureS4/FigureS4_DRCs_res063-093.pdf", width = 250, height = 100)
plot_grid(plotlist = out.31_60,
          ncol = 30,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()

pdf("../../results/FigureS4/FigureS4_DRCs_res93-122.pdf", width = 250, height = 100)
plot_grid(plotlist = out.61_90,
          ncol = 30,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()

pdf("../../results/FigureS4/FigureS4_DRCs_res123-152.pdf", width = 250, height = 100)
plot_grid(plotlist = out.91_120,
          ncol = 30,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()

pdf("../../results/FigureS4/FigureS4_DRCs_res153-182.pdf", width = 250, height = 100)
plot_grid(plotlist = out.121_150,
          ncol = 30,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()

pdf("../../results/FigureS4/FigureS4_DRCs_res183-209.pdf", width = 225, height = 100)
plot_grid(plotlist = out.151_180,
          ncol = 27,
          nrow = 21,
          align = "v", axis = "tb", byrow = F)
dev.off()


## 4. Version ##
################

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
# [1] wesanderson_0.3.7 cowplot_1.1.3     ggtext_0.1.2      ggplot2_3.5.1     reshape_0.8.9     scales_1.3.0      drc_3.0-1        
# [8] MASS_7.3-64       stringr_1.5.1    
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-2      gtable_0.3.6      dplyr_1.1.4       compiler_4.4.1    gtools_3.9.5      tidyselect_1.2.1  plotrix_3.8-4    
# [8] Rcpp_1.0.14       xml2_1.3.6        splines_4.4.1     lattice_0.22-6    TH.data_1.1-3     R6_2.6.1          plyr_1.8.9       
# [15] generics_0.1.3    Formula_1.2-5     tibble_3.2.1      car_3.1-3         munsell_0.5.1     pillar_1.10.1     rlang_1.1.5      
# [22] multcomp_1.4-28   stringi_1.8.4     cli_3.6.4         withr_3.0.2       magrittr_2.0.3    gridtext_0.1.5    grid_4.4.1       
# [29] rstudioapi_0.17.1 mvtnorm_1.3-3     sandwich_3.1-1    lifecycle_1.0.4   vctrs_0.6.5       glue_1.8.0        codetools_0.2-20 
# [36] zoo_1.8-12        survival_3.8-3    abind_1.4-8       carData_3.0-5     colorspace_2.1-1  pkgconfig_2.0.3   tools_4.4.1   