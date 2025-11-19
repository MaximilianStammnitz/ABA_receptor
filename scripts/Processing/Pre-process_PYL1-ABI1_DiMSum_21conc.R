# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#############################################
## Pre-processing of PYL1-ABI1 DiMSum data ##
#############################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "drc", "reshape")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input DiMSum tables ##
############################

## Input function
load_blocks <- function(conc){
  load(paste0("../../data/DiMSum/PYL1-ABI1_21conc/conc_", conc, "_fitness_replicates.RData"))
  data <- rbind(all_variants, synonymous)
  data <- data[-which(data$WT == T)[-1],] ## remove duplicated WT
  return(data)}

PYL1_block <- list("1" = load_blocks("01"), "2" = load_blocks("02"),
                   "3" = load_blocks("03"), "4" = load_blocks("04"),
                   "5" = load_blocks("05"), "6" = load_blocks("06"),
                   "7" = load_blocks("07"), "8" = load_blocks("08"),
                   "9" = load_blocks("09"), "10" = load_blocks("10"),
                   "11" = load_blocks("11"), "12" = load_blocks("12"),
                   "13" = load_blocks("13"), "14" = load_blocks("14"),
                   "15" = load_blocks("15"), "16" = load_blocks("16"),
                   "17" = load_blocks("17"), "18" = load_blocks("18"),
                   "19" = load_blocks("19"), "20" = load_blocks("20"),
                   "21" = load_blocks("21"))

## Clean up environment
rm(packages, install_if_missing, load_blocks)


## 2. Normalise growth rates, using stops & synonymous WT variants ##
#####################################################################

## Error-weighted mean of growth rate estimates
PYL1_block <- lapply(PYL1_block, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
PYL1_block <- lapply(PYL1_block, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})

## Error-weighted mean of stop mutations
stops.block <- lapply(PYL1_block, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## WT/synonymous muts. growth rate for centering (obtain scaled growth rate)
syn.block <- lapply(PYL1_block, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## Calculate coefficients for linear transformation
scaling <- cbind(data.frame(rbind(do.call(c, stops.block), do.call(c, syn.block))))
rownames(scaling) <- c("stops", "syn")
colnames(scaling) <- 1:21

## Linear transformation (stop mode consistently to 0)
coefficients <- matrix(NA, ncol = 21, nrow = 2)
colnames(coefficients) <-  colnames(scaling)
for (i in 1:21){

    tmp.lm <- lm(formula = c(0, scaling[2,i]) ~ scaling[,i])
    tmp.lm <- summary(tmp.lm)
    coefficients[1,i] <- tmp.lm$coefficients[[2]]
    coefficients[2,i] <- tmp.lm$coefficients[[1]]  
  
}
for (i in 1:length(PYL1_block)){
  
  PYL1_block[[i]]$gr_normalised <- PYL1_block[[i]]$growthrate*coefficients[1,i] + coefficients[2,i]
  PYL1_block[[i]]$gr_sigma_normalised <- PYL1_block[[i]]$growthrate_sigma*coefficients[1,i]
  
}

## Clean up environment
rm(i, tmp.lm)


## 3. Combine PYL1 blocks ##
############################

## Only keep mutations with AA Hamming distance of 1
PYL1_block <- lapply(PYL1_block, function(x){x <- x[which(x$Nham_aa == 1 | x$Nham_aa == 0),]; return(x)})

## Annotate the individual mutation types
PYL1_block <- lapply(PYL1_block, function(x){x <- cbind("Pos" = rep(NA, nrow(x)),
                                                        "WT_AA" = rep(NA, nrow(x)),
                                                        "Mut" = rep(NA, nrow(x)),
                                                        x); return(x)})
for(j in 1:length(PYL1_block)){
  
  for(i in 1:nrow(PYL1_block[[j]])){
    
    ### Do not act on WT seqs
    if(PYL1_block[[j]][i,"Nham_aa"] == 0){
      
      next
      
    }
    
    ### Find the position and mutation
    WT.seq <- PYL1_block[[j]][which(PYL1_block[[j]][,"WT"] == T)[1],"aa_seq"]
    WT.seq <- as.character(str_split_fixed(WT.seq, "", 32))
    mut.seq <- PYL1_block[[j]][i,"aa_seq"]
    mut.seq <- as.character(str_split_fixed(mut.seq, "", 32))
    PYL1_block[[j]][i,"Pos"] <- which(mut.seq != WT.seq)
    PYL1_block[[j]][i,"WT_AA"] <- WT.seq[PYL1_block[[j]][i,"Pos"]]
    PYL1_block[[j]][i,"Mut"] <- mut.seq[PYL1_block[[j]][i,"Pos"]]
    
  } 
  
}

## Positional adjustment with respect to the full-length PYL1 protein
PYL1_block <- lapply(PYL1_block, function(x){x$Pos <- x$Pos + 117; return(x)})

## WT specification
PYL1_block <- lapply(PYL1_block, function(x){x[which(x[,"Nham_aa"] == 0),"Pos"] <- "WT"; x[which(x[,"Nham_aa"] == 0),"WT_AA"] <- "WT"; x[which(x[,"Nham_aa"] == 0),"Mut"] <- "WT"; return(x)})

## Dosages
dosages <- rep(NA, 21)
dosages[1] <- 2500
dosages[21] <- 0
for(i in 2:20){
  dosages[i] <- 2500 * (1/2.5)^c(i-1)
}
names(PYL1_block) <- dosages


## 4. Scale growth rates by wildtype B[inf] ##
##############################################

## Extract WT aa variants (incl. synonymous ones)
PYL1.ABI1.WT <- lapply(PYL1_block, function(x){y <- x[which(x[,"Nham_aa"] == 0),]; return(y)})

## Build a dose-response matrix for each nucleotide sequence
PYL1.ABI1.WT.mat <- matrix(NA, ncol = 21, nrow = length(unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))))
colnames(PYL1.ABI1.WT.mat) <- names(PYL1.ABI1.WT)
rownames(PYL1.ABI1.WT.mat) <- unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))
for(i in 1:21){
  PYL1.ABI1.WT.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_normalised"]
}

## Remove synonymous variants not fully covered across (+)-ABA concentrations
PYL1.ABI1.WT.mat <- na.omit(PYL1.ABI1.WT.mat)

## Generate the dose response curve of the (nucleotide-level) WT
WT.PYL1.drc <- cbind(PYL1.ABI1.WT.mat[1,],colnames(PYL1.ABI1.WT.mat))
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

## Scale all growth rates to the wildtype B[inf]
PYL1_block <- lapply(PYL1_block, function(x){x$gr_normalised_WTscaled <- 100*x$gr_normalised/WT.PYL1.drc.par["B[inf]"]; return(x)})
PYL1_block <- lapply(PYL1_block, function(x){x$gr_sigma_normalised_WTscaled <- 100*x$gr_sigma_normalised/WT.PYL1.drc.par["B[inf]"]; return(x)})

## Save PYL1.ABI1 list as an .Rdata file
save(PYL1_block, file = "../../data/DiMSum/PYL1-ABI1_21conc/PYL1-ABI1_preprocessed_21conc.RData")


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
# [1] reshape_0.8.10 drc_3.0-1      MASS_7.3-65    stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        TH.data_1.1-4      Formula_1.2-5      stringi_1.8.7      car_3.1-3          gtools_3.9.5       zoo_1.8-14         glue_1.8.0        
# [11] plyr_1.8.9         scales_1.4.0       grid_4.5.1         carData_3.0-5      abind_1.4-8        mvtnorm_1.3-3      lifecycle_1.0.4    compiler_4.5.1     multcomp_1.4-28    codetools_0.2-20  
# [21] RColorBrewer_1.1-3 sandwich_3.1-1     Rcpp_1.1.0         farver_2.1.2       lattice_0.22-7     R6_2.6.1           splines_4.5.1      magrittr_2.0.4     Matrix_1.7-4       tools_4.5.1       
# [31] plotrix_3.8-4      survival_3.8-3  