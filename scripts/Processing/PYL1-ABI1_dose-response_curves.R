# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#########################################################
## Generate Hill dose-response models for all variants ##
#########################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "drc", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-process DiMSum data ##
################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

### remove synonymous mutations, except for the main wildtype
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[-which(x$Nham_aa == 0 & x$Nham_nt > 0),]; return(x)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[-which(x$WT == T & x$nt_seq != "actcaagacgaattcacccaactctcccaatcaatcgccgagttccacacgtaccaactcggtaacggccgttgctcatctctcctagctcagcgaatccacgcgccgccggaaacagtatggtccgtggtgagacgt"),]; return(x)})

## Summarise WT-rescaled binding fitness for each variant vs. ABA combination
## Matrix format
PYL1.ABI1.summary <- matrix(NA, nrow = 177*21, ncol = 12)
colnames(PYL1.ABI1.summary) <- names(PYL1.ABI1)
ids <- c()
for(i in 1:177){
  tmp.pos <- rep(str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"], "", 177)[i], 21)
  ids <- c(ids, paste0(tmp.pos, i + 32))
}
rownames(PYL1.ABI1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
                                                 "I", "F", "Y", "W", "K",
                                                 "R", "H", "D", "E", "S",
                                                 "T", "C", "N", "Q", "P",
                                                 "*"),177))

## Fill
WT <- str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"],"",177)[1,]
for(j in 1:length(PYL1.ABI1)){
  
  tmp.vars <- str_split_fixed(PYL1.ABI1[[j]]$aa_seq,"",177)
  
  for(i in 1:nrow(tmp.vars)){
    
    ### only take the "best" WT
    if(all(tmp.vars[i,] == WT)){
      
      PYL1.ABI1.summary[paste0(WT, 33:209, WT),j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
      next
      
    }
    
    pos.mut <- which(tmp.vars[i,] != WT)
    
    ### skip if higher-order mutant
    if(length(pos.mut) > 1){
      
      next
      
    }
    
    ### otherwise obtain the corresponding fitness
    tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
    PYL1.ABI1.summary[tmp.mut,j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
    
  }
  
}
rm(i,j,tmp.vars,ids,pos.mut,tmp.mut,tmp.pos)

## Subset genotypes to only occur once (i.e. only one WT)
WTs.pos <- match(paste0(WT, 33:209, WT), rownames(PYL1.ABI1.summary))
PYL1.ABI1.summary.nonWT <- PYL1.ABI1.summary[-WTs.pos[-1],]
PYL1.ABI1.summary.nonWT <- PYL1.ABI1.summary.nonWT[c(WTs.pos[1],1:c(WTs.pos[1]-1),c(WTs.pos[1]+1):nrow(PYL1.ABI1.summary.nonWT)),]
rownames(PYL1.ABI1.summary.nonWT)[1] <- "WT"


## 2. Calculate and plot (+)-ABA dose-response curves using the Hill equation ##
################################################################################

## Calculate and plot all the dose-response curves with the 4-parametric Hill model
pdf("../../results/All_PYL1-ABI1_dose-response_curves.pdf", height = 15, width = 18)

## setup
failed.drc.Hill <- c()
parameters.Hill <- matrix(NA, nrow = nrow(PYL1.ABI1.summary.nonWT), ncol = 16)
colnames(parameters.Hill) <- c("Hill", "B[0]", "B[inf]", "EC50",
                               "Hill SE", "B[0] SE", "B[inf] SE", "EC50 SE",
                               "Hill P", "B[0] P", "B[inf] P", "EC50 P",
                               "Data points", "AIC", "Residual var", "R^2")
rownames(parameters.Hill) <- rownames(PYL1.ABI1.summary.nonWT)

## loop over all genotypes, this may take a few minutes ...
for(i in 1:nrow(parameters.Hill)){

  print(i)

  ### DRC modelling, using the R DRC package (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021)

  ### output parameter translations:
  ### b: Hill coefficient, i.e. steepness of the curve (n)
  ### c: basal binding fitness (B[0])
  ### d: saturated binding fitness (B[inf])
  ### e: curve inflection point (EC50)
  ### model: binding(ABA conc.) = c + ((d-c)/(1 + (e/x)^b))

  ### fetch genotype's data
  drc.summary.tmp <- cbind(PYL1.ABI1.summary.nonWT[i,],
                           colnames(PYL1.ABI1.summary.nonWT))
  class(drc.summary.tmp) <- "numeric"
  drc.summary.tmp <- as.data.frame(drc.summary.tmp)
  colnames(drc.summary.tmp) <- c("binding", "concentration")

  ### freely estimate the wildtype parameters, used to initialise runs for PYL1 mutants
  if(i == 1){

    WT.PYL1.drc.Hill <- drm(drc.summary.tmp$binding ~ drc.summary.tmp$concentration,
                            fct = LL.4(fixed = c(NA, NA, NA, NA),
                                       names = c("Hill", "B[0]", "B[inf]", "EC50")),
                            type = 'continuous')

    par.WT.PYL1.drc.Hill <- WT.PYL1.drc.Hill$fit$par
    names(par.WT.PYL1.drc.Hill) <- c("Hill", "B[0]", "B[inf]", "EC50")
    rm(WT.PYL1.drc.Hill)
  }

  ### skip variant if there are less than 5 data points
  parameters.Hill[i,"Data points"] <- sum(is.na(drc.summary.tmp$binding) == F)
  if(sum(is.na(drc.summary.tmp$binding) == F) < 5){next}

  ### skip variant if DRC fails to build
  skip_to_next <- FALSE
  tryCatch({

    tmp.PYL1.drc.Hill <- drm(drc.summary.tmp$binding ~ drc.summary.tmp$concentration,
                             fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                             type = 'continuous', start = par.WT.PYL1.drc.Hill)},
    error = function(e){skip_to_next <<- TRUE})

  if(skip_to_next == T) {

    failed.drc.Hill <- c(failed.drc.Hill, i)
    next

  }else{

    ### fetch Hill outputs
    parameters.Hill[i,1:4] <- summary(tmp.PYL1.drc.Hill)$coefficients[,1]
    parameters.Hill[i,5:8] <- summary(tmp.PYL1.drc.Hill)$coefficients[,2]
    parameters.Hill[i,9:12] <- summary(tmp.PYL1.drc.Hill)$coefficients[,4]
    if(parameters.Hill[i,1] < 0){
      parameters.Hill[i,1] <- -parameters.Hill[i,1] ## invert Hill parameter output
    }else{
      parameters.Hill[i,2:3] <- parameters.Hill[i,c(3,2)]
      parameters.Hill[i,10:11] <- parameters.Hill[i,c(11,10)]
    }
    parameters.Hill[i,14:15] <- as.numeric(mselect(tmp.PYL1.drc.Hill, icfct = AIC)[c(2,4)])

    ### calculation of non-linear R^2
    predicted_values <- predict(tmp.PYL1.drc.Hill)
    rss <- sum((drc.summary.tmp$binding - predicted_values)^2)
    tss <- sum((drc.summary.tmp$binding - mean(drc.summary.tmp$binding))^2)
    tmp.r_squared.Hill <- 1 - (rss/tss)
    parameters.Hill[i,16] <- tmp.r_squared.Hill
    rm(predicted_values, rss, tss, tmp.r_squared.Hill)

    ### predict the full, smoothened curve (using 1000 data points) and confidence interval
    tmp.PYL1.drc.Hill.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
    tmp.PYL1.drc.Hill.predict <- predict(tmp.PYL1.drc.Hill,
                                         newdata = tmp.PYL1.drc.Hill.predict.newdata,
                                         interval = "confidence")

    #### new data with predictions
    tmp.PYL1.drc.Hill.predict.newdata$p <- tmp.PYL1.drc.Hill.predict[,1]
    tmp.PYL1.drc.Hill.predict.newdata$pmin <- tmp.PYL1.drc.Hill.predict[,2]
    tmp.PYL1.drc.Hill.predict.newdata$pmax <- tmp.PYL1.drc.Hill.predict[,3]

    ### plot curve
    drc.summary.tmp$concentration[12] <- 9.062741e-03/3.5/3.5/3.5 ## "0-conc." positioning for log scale
    print(ggplot(drc.summary.tmp, aes(x = concentration, y = binding)) +
            geom_ribbon(data = tmp.PYL1.drc.Hill.predict.newdata,
                        aes(x = conc, y = p, ymin = pmin, ymax = pmax),
                        alpha = 0.2, fill = "grey50") +
            geom_point(data = drc.summary.tmp, aes(x = concentration, y = binding),
                       color = "black", size = 8, shape = 16) +
            geom_line(data = tmp.PYL1.drc.Hill.predict.newdata, aes(x = conc, y = p), linewidth = 1.5) +
            scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                          labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                          limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
            scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limits = c(-5, 115)) +
            coord_cartesian(ylim = c(-15, 115)) +
            annotate("text",
                     x = 9.062741e-03/3.5/3.5/3.5,
                     y = 105,
                     label = expr_text(bquote(italic(R)^2 == .(format(parameters.Hill[i,"R^2"], digits = 2)))),
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
                 y = "Relative PYL1/ABI1 Binding (library sequencing)",
                 title = rownames(PYL1.ABI1.summary.nonWT)[i],
                 subtitle = paste0("Hill parameters: B[0] = ", round(parameters.Hill[i,"B[0]"],3), if(!is.na(parameters.Hill[i,"B[0] P"])){if(parameters.Hill[i,"B[0] P"] < 0.05){"*"}},
                                   ", B[inf] = ", round(parameters.Hill[i,"B[inf]"],3), if(!is.na(parameters.Hill[i,"B[inf] P"])){if(parameters.Hill[i,"B[inf] P"] < 0.05){"*"}},
                                   ", EC50 = ", round(parameters.Hill[i,"EC50"],3), if(!is.na(parameters.Hill[i,"EC50 P"])){if(parameters.Hill[i,"EC50 P"] < 0.05){"*"}},
                                   ", Hill = ", round(parameters.Hill[i,"Hill"],3), if(!is.na(parameters.Hill[i,"Hill P"])){if(parameters.Hill[i,"Hill P"] < 0.05){"*"}})))

    rm(tmp.PYL1.drc.Hill, tmp.PYL1.drc.Hill.predict, tmp.PYL1.drc.Hill.predict.newdata, drc.summary.tmp)

  }

}

dev.off()

## Save parameters.Hill dataframe as an .Rdata file
save(parameters.Hill, file = "../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")


## 3. Version ##
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
# [1] rlang_1.1.6   ggtext_0.1.2  ggplot2_4.0.0 drc_3.0-1     MASS_7.3-65   scales_1.4.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       crayon_1.5.3       dplyr_1.1.4        compiler_4.5.1     gtools_3.9.5       Rcpp_1.1.0         tidyselect_1.2.1  
# [9] plotrix_3.8-4      xml2_1.4.0         splines_4.5.1      lattice_0.22-7     TH.data_1.1-4      R6_2.6.1           commonmark_2.0.0   generics_0.1.4    
# [17] Formula_1.2-5      tibble_3.3.0       car_3.1-3          pillar_1.11.0      RColorBrewer_1.1-3 multcomp_1.4-28    stringi_1.8.7      xfun_0.53         
# [25] litedown_0.7       S7_0.2.0           cli_3.6.5          withr_3.0.2        magrittr_2.0.3     gridtext_0.1.5     grid_4.5.1         rstudioapi_0.17.1 
# [33] mvtnorm_1.3-3      markdown_2.0       sandwich_3.1-1     lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2       codetools_0.2-20  
# [41] zoo_1.8-14         survival_3.8-3     abind_1.4-8        carData_3.0-5      tools_4.5.1        pkgconfig_2.0.3
