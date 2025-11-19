# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 29.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 5A - PYL1-PYL1 TECAN validations ##
###########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthcurver", "scales", 
              "reshape", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input PYL1-PYL1 data ##
#############################

## import

load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]


## 2. Import TECAN data ##
##########################

## raw data import
PYL1.PYL1.TECAN <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-PYL1-PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[22:407,]
rownames(PYL1.PYL1.TECAN) <- PYL1.PYL1.TECAN[,1]
colnames(PYL1.PYL1.TECAN) <- PYL1.PYL1.TECAN[1,]
PYL1.PYL1.TECAN <- PYL1.PYL1.TECAN[,-1]
PYL1.PYL1.TECAN <- PYL1.PYL1.TECAN[-c(1:2),]
colnames(PYL1.PYL1.TECAN) <- as.numeric(colnames(PYL1.PYL1.TECAN))/3600 ## convert to hours
class(PYL1.PYL1.TECAN) <- "numeric"
PYL1.PYL1.TECAN <- as.data.frame(PYL1.PYL1.TECAN)

## assign mutants
setup <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-PYL1-PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[2:18,1:25]
colnames(setup) <- setup[1,]
setup <- setup[-1,]
rownames(setup) <- setup[,1]
setup <- setup[,-1]


## 3. Process using growthcurver ##
###################################

PYL1.PYL1.TECAN <- cbind("time" = colnames(PYL1.PYL1.TECAN), 
                         t(PYL1.PYL1.TECAN))
class(PYL1.PYL1.TECAN) <- "numeric"
PYL1.PYL1.TECAN <- as.data.frame(PYL1.PYL1.TECAN)
colnames(PYL1.PYL1.TECAN)[26] <- "blank"

## run
PYL1.PYL1.TECAN.curves <- SummarizeGrowthByPlate(PYL1.PYL1.TECAN[1:275,], 
                                                 t_trim = 75,
                                                 bg_correct = "blank",
                                                 plot_fit = F)
PYL1.PYL1.TECAN.curves <- PYL1.PYL1.TECAN.curves[-384,]

## label wells by genotype
for(i in 1:nrow(PYL1.PYL1.TECAN.curves)){
  PYL1.PYL1.TECAN.curves[i,1] <- setup[substr(PYL1.PYL1.TECAN.curves[i,1], 1, 1),
                                       substr(PYL1.PYL1.TECAN.curves[i,1], 2, 
                                              nchar(PYL1.PYL1.TECAN.curves[i,1]))]
}

## remove unannotated ones (empty wells)
PYL1.PYL1.TECAN.curves <- PYL1.PYL1.TECAN.curves[-which(is.na(PYL1.PYL1.TECAN.curves$sample)),]

## categorise
PYL1.PYL1.TECAN.curves.auc <- vector(mode = "list", length = 40)
names(PYL1.PYL1.TECAN.curves.auc) <- unique(c(setup))[-2]
names(PYL1.PYL1.TECAN.curves.auc) <- c(names(PYL1.PYL1.TECAN.curves.auc)[1],
                                       names(PYL1.PYL1.TECAN.curves.auc)[-1][order(as.numeric(substr(names(PYL1.PYL1.TECAN.curves.auc)[-1], 2, nchar(names(PYL1.PYL1.TECAN.curves.auc)[-1]) - 1)))])
for(i in 1:length(PYL1.PYL1.TECAN.curves.auc)){
  PYL1.PYL1.TECAN.curves.auc[[i]] <- which(PYL1.PYL1.TECAN.curves$sample == names(PYL1.PYL1.TECAN.curves.auc)[i])
  PYL1.PYL1.TECAN.curves.auc[[i]] <- PYL1.PYL1.TECAN.curves$auc_e[PYL1.PYL1.TECAN.curves.auc[[i]]]
}

## clean-up
rm(i,packages, setup)


## 4. Plot ##
#############

## summarise data
out.all <- matrix(NA, nrow = 40, ncol = 4)
colnames(out.all) <- c("bulk_mean", "bulk_sd", "tecan_mean", "tecan_sd")
rownames(out.all) <- names(PYL1.PYL1.TECAN.curves.auc)

### add bulk library mean/sd
for(i in 1:nrow(out.all)){
  
  if(i == 1){
    
    tmp.id <- which(PYL1.PYL1.0uM.ABA$WT == T)
    
  }else{
    
    tmp.name <- rownames(out.all)[i]
    tmp.id <- match(tmp.name, paste0(PYL1.PYL1.0uM.ABA[,"WT_AA"],
                                     PYL1.PYL1.0uM.ABA[,"Pos"],
                                     PYL1.PYL1.0uM.ABA[,"Mut"]))
    
  }
  
  out.all[i,"bulk_mean"] <- PYL1.PYL1.0uM.ABA[tmp.id,"gr_normalised_WTscaled"]
  out.all[i,"bulk_sd"] <- PYL1.PYL1.0uM.ABA[tmp.id,"gr_sigma_normalised_WTscaled"]

}

### add TECAN mean/sd
out.all[,"tecan_mean"] <- sapply(PYL1.PYL1.TECAN.curves.auc, mean)
out.all[,"tecan_sd"] <- sapply(PYL1.PYL1.TECAN.curves.auc, sd)

### rescale
out.all[,"tecan_sd"] <- out.all[,"tecan_sd"] * c(out.all["WT","bulk_mean"] / out.all["WT","tecan_mean"])
out.all[,"tecan_mean"] <- out.all[,"tecan_mean"] * c(out.all["WT","bulk_mean"] / out.all["WT","tecan_mean"])
out.all <- as.data.frame(out.all)

## to calculate Pearson's coefficients
r.out <- cor(x = out.all[-1,"bulk_mean"],
         y = out.all[-1,"tecan_mean"],
         method = "pearson", use = "complete.obs")
p.out <- summary(lm(out.all[-1,"bulk_mean"] ~ out.all[-1,"tecan_mean"]))$coefficients[2,4]

pdf("../../results/FigureS5/FigureS5A_abundancePCA_validation.pdf", height = 15, width = 18)

out.S5A <- ggplot(out.all, aes(x = `tecan_mean`, y = `bulk_mean`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-1000, 1000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-1000, 1000)) +
  coord_cartesian(xlim = c(-20, 150), ylim = c(-20, 150)) +
  geom_errorbar(aes(xmin = `tecan_mean` - `tecan_sd`, xmax = `tecan_mean` + `tecan_sd`),
                color = "black", linewidth = 0.5, height = 0) +
  geom_errorbar(aes(ymin = `bulk_mean` - `bulk_sd`, ymax = `bulk_mean` + `bulk_sd`),
                color = "black", linewidth = 0.5, width = 0) +
  geom_point(data = out.all,
             mapping = aes(x = `tecan_mean`, y = `bulk_mean`),
             color = "black", size = 5, shape = 16) +
  geom_smooth(data = out.all[-1,],
              mapping = aes(x = `tecan_mean`, y = `bulk_mean`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = "grey50") +
  annotate("text",
           x = -20,
           y = 140,
           label = expr_text(bquote(italic(r) == .(format(r.out, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.out, digits = 2, nsmall = 2)))),
           parse = T,
           hjust = 0, size = 20, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),        
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Relative PYL1-PYL1 Abundance (individual mutants, AUC)",
       y = "Relative PYL1-PYL1 (library sequencing)")

print(out.S5A)

dev.off()


## 5. Version ##
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
# [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
# [1] rlang_1.1.6        ggtext_0.1.2       ggplot2_4.0.0     
# [4] reshape_0.8.10     scales_1.4.0       growthcurver_0.3.1
# [7] readxl_1.4.5       stringr_1.5.1     
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-4       gtable_0.3.6       dplyr_1.1.4       
# [4] compiler_4.5.1     crayon_1.5.3       tidyselect_1.2.1  
# [7] Rcpp_1.1.0         minpack.lm_1.2-4   xml2_1.4.0        
# [10] splines_4.5.1      lattice_0.22-7     R6_2.6.1          
# [13] plyr_1.8.9         labeling_0.4.3     generics_0.1.4    
# [16] tibble_3.3.0       pillar_1.11.0      RColorBrewer_1.1-3
# [19] stringi_1.8.7      S7_0.2.0           cli_3.6.5         
# [22] withr_3.0.2        magrittr_2.0.3     mgcv_1.9-3        
# [25] grid_4.5.1         gridtext_0.1.5     rstudioapi_0.17.1 
# [28] nlme_3.1-168       lifecycle_1.0.4    vctrs_0.6.5       
# [31] glue_1.8.0         farver_2.1.2       cellranger_1.1.0  
# [34] tools_4.5.1        pkgconfig_2.0.3  
