# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

################################################
## Script Figure 4F - aPCA vs ThermoMPNN ∆∆Gf ##
################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "ggplot2", "ggtext", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input and filter dose response curve parameters ##
########################################################

## Input 12-dose measurements
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

## remove stops and any mutants with incomplete curves
PYL1.summary.nonWT <- PYL1.summary.nonWT[-grep("[*]", rownames(PYL1.summary.nonWT)),]
PYL1.summary.nonWT <- PYL1.summary.nonWT[-sort(unique(do.call(c,apply(PYL1.summary.nonWT[,1:12], 2, function(x){which(is.na(x) == T)})))),1:12]
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)

## add aPCA data to the table
load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]
PYL1.summary.nonWT <- cbind(PYL1.summary.nonWT, "aPCA" = rep(NA, nrow(PYL1.summary.nonWT)))
for(i in 1:nrow(PYL1.summary.nonWT)){
  
  print(i)
  tmp.res <- rownames(PYL1.summary.nonWT)[i]
  if(tmp.res == "WT"){
    
    PYL1.summary.nonWT[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT == T)[1],"gr_normalised_WTscaled"]

  }else{
    
    tmp.pos <- as.numeric(substr(tmp.res, 2, nchar(tmp.res)-1))
    tmp.WT <- substr(tmp.res, 1, 1)
    tmp.MUT <- substr(tmp.res, nchar(tmp.res), nchar(tmp.res))
    rep <- which(PYL1.PYL1.0uM.ABA$Pos == tmp.pos & PYL1.PYL1.0uM.ABA$Mut == tmp.MUT)
    if(length(rep) != 0){
      
      PYL1.summary.nonWT[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[rep,"gr_normalised_WTscaled"]
      
    }

  }
  
}


## 2. PCA ##
############

PCA <- prcomp(PYL1.summary.nonWT[,1:12], scale = T, center = T)
sum.PCA <- summary(PCA)
y.PCA.components.vars <- sum.PCA$importance[2,]*100
y.PCA.components.vars <- round(y.PCA.components.vars, 1)


## 3. Input and in silico ddGf predictions ## 
#############################################

## ThermoMPNN ddG
ThermoMPNN_raw <- read.csv("../../data/Predictions/ThermoMPNN/output/ThermoMPNN_inference_3kay_homodimer_chainApredictions.csv")
ThermoMPNN <- cbind(paste0(ThermoMPNN_raw$wildtype, ThermoMPNN_raw$position + 35, ThermoMPNN_raw$mutation), ThermoMPNN_raw$ddG_pred)
ThermoMPNN <- ThermoMPNN[-grep("^S35", ThermoMPNN[,1]),]
ThermoMPNN <- ThermoMPNN[-grep("^R210", ThermoMPNN[,1]),]
ThermoMPNN <- ThermoMPNN[-which(ThermoMPNN[,1] %in% paste0(unique(substr(ThermoMPNN[,1], 1, nchar(ThermoMPNN[,1]) - 1)), substr(unique(substr(ThermoMPNN[,1], 1, nchar(ThermoMPNN[,1]) - 1)), 1, 1))),]
colnames(ThermoMPNN) <- c("Mutation", "ddGf")
rm(ThermoMPNN_raw)


## 4. Correlate aPCA vs. ThermoMPNN ∆∆Gf ##
###########################################

ddGf.out <- ThermoMPNN
ddGf.out <- cbind(ddGf.out, 
                  "aPCA" = PYL1.summary.nonWT[match(ddGf.out[,"Mutation"], rownames(PYL1.summary.nonWT)),"aPCA"])
ddGf.out <- ddGf.out[,-1]
class(ddGf.out) <- "numeric"
ddGf.out <- as.data.frame(ddGf.out)
rho_ddGf <- cor(x = ddGf.out$ddGf, y = ddGf.out$aPCA, 
                method = "spearman", use = "complete.obs")

pdf('../../results/Figure4/Figure4F_aPCA_vs_ddGf.pdf', height = 10, width = 14)

out.4F <- ggplot(ddGf.out, aes(x = ddGf, y = aPCA)) +
  geom_density2d_filled(aes(x = ddGf, y = aPCA), bins = 6, inherit.aes = F) +
  scale_fill_grey(start = 1, end = 0.1) +
  annotate("text",
           x = -0.5,
           y = -12.5,
           label = expr_text(bquote(rho == .(format(rho_ddGf, digits = 2)))),
           parse = T,
           hjust = 0, size = 20, color = "black") +
  scale_x_continuous(breaks = seq(from = 0, to = 3, length.out = 4), limit = c(-0.5, 3.5)) +
  scale_y_continuous(breaks = seq(from = 0, to = 100, length.out = 6), limit = c(-20, 125)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 30),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 2),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(3, 2, 0.7, 2),"cm")) +
  labs(x = expression(atop("Folding "*Delta*Delta*"G by ThermoMPNN", "(kcal mol"^{-1}*")")),
       y = "Relative PYL1 Abundance")

print(out.4F)

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
# [1] rlang_1.1.6   ggtext_0.1.2  ggplot2_4.0.0 scales_1.4.0  stringr_1.5.2
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          stringi_1.8.7      generics_0.1.4     S7_0.2.0           glue_1.8.0         isoband_0.2.7      gridtext_0.1.5     grid_4.5.1        
# [11] tibble_3.3.0       MASS_7.3-65        lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0         pkgconfig_2.0.3    farver_2.1.2       R6_2.6.1          
# [21] tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     tools_4.5.1        withr_3.0.2        gtable_0.3.6       xml2_1.4.0     