# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 21.10.2025
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
load_blocks <- function(block, conc){
  load(paste0("../../data/DiMSum/PYL1-ABI1/block", block, "_", conc, "_fitness_replicates.RData"))
  data <- rbind(all_variants, synonymous)
  data <- data[-which(data$WT == T)[-1],] ## remove duplicated WT
  return(data)}

PYL1_block1 <- list("1" = load_blocks(1, "01"), "2" = load_blocks(1, "02"),
                    "3" = load_blocks(1, "03"), "4" = load_blocks(1, "04"),
                    "5" = load_blocks(1, "05"), "6" = load_blocks(1, "06"),
                    "7" = load_blocks(1, "07"), "8" = load_blocks(1, "08"),
                    "9" = load_blocks(1, "09"), "10" = load_blocks(1, "10"),
                    "11" = load_blocks(1, "11"), "12" = load_blocks(1, "12"))

PYL1_block2 <- list("1" = load_blocks(2, "01"), "2" = load_blocks(2, "02"),
                    "3" = load_blocks(2, "03"), "4" = load_blocks(2, "04"),
                    "5" = load_blocks(2, "05"), "6" = load_blocks(2, "06"),
                    "7" = load_blocks(2, "07"), "8" = load_blocks(2, "08"),
                    "9" = load_blocks(2, "09"), "10" = load_blocks(2, "10"),
                    "11" = load_blocks(2, "11"), "12" = load_blocks(2, "12"))

### Manually removing Y8*, likely this is a DiMSum count/filtering artefact unique to only two concentrations:
### Variant appears fully viable/neutral and skews the normalisation process
PYL1_block2[[9]] <- PYL1_block2[[9]][-which(PYL1_block2[[9]]$aa_seq == "FDRPQI*KHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD"),]
PYL1_block2[[12]] <- PYL1_block2[[12]][-which(PYL1_block2[[12]]$aa_seq == "FDRPQI*KHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD"),]

PYL1_block3 <- list("1" = load_blocks(3, "01"), "2" = load_blocks(3, "02"),
                    "3" = load_blocks(3, "03"), "4" = load_blocks(3, "04"),
                    "5" = load_blocks(3, "05"), "6" = load_blocks(3, "06"),
                    "7" = load_blocks(3, "07"), "8" = load_blocks(3, "08"),
                    "9" = load_blocks(3, "09"), "10" = load_blocks(3, "10"),
                    "11" = load_blocks(3, "11"), "12" = load_blocks(3, "12"))

PYL1_block4 <- list("1" = load_blocks(4, "01"), "2" = load_blocks(4, "02"),
                    "3" = load_blocks(4, "03"), "4" = load_blocks(4, "04"),
                    "5" = load_blocks(4, "05"), "6" = load_blocks(4, "06"),
                    "7" = load_blocks(4, "07"), "8" = load_blocks(4, "08"),
                    "9" = load_blocks(4, "09"), "10" = load_blocks(4, "10"),
                    "11" = load_blocks(4, "11"), "12" = load_blocks(4, "12"))

## Clean up environment
rm(packages, install_if_missing, load_blocks)


## 2. Normalise growth rates between blocks, using stops & synonymous variants ##
#################################################################################

## Error-weighted mean of growth rate estimates
PYL1_block1 <- lapply(PYL1_block1, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
PYL1_block1 <- lapply(PYL1_block1, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})
PYL1_block2 <- lapply(PYL1_block2, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
PYL1_block2 <- lapply(PYL1_block2, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})
PYL1_block3 <- lapply(PYL1_block3, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
PYL1_block3 <- lapply(PYL1_block3, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})
PYL1_block4 <- lapply(PYL1_block4, function(x){x$gr_over_sigmasquared <- x$growthrate/(x$growthrate_sigma)**2; return(x)})
PYL1_block4 <- lapply(PYL1_block4, function(x){x$one_over_sigmasquared <- 1/(x$growthrate_sigma)**2; return(x)})

## Error-weighted mean of stop mutations
stops.block1 <- lapply(PYL1_block1, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
stops.block2 <- lapply(PYL1_block2, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
stops.block3 <- lapply(PYL1_block3, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
stops.block4 <- lapply(PYL1_block4, function(x){x <- x[which(x[,"STOP"] == T),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE) / sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## WT/synonymous muts. growth rate for centering (obtain scaled growth rate)
syn.block1 <- lapply(PYL1_block1, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
syn.block2 <- lapply(PYL1_block2, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
syn.block3 <- lapply(PYL1_block3, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})
syn.block4 <- lapply(PYL1_block4, function(x){x <- x[which(x$Nham_aa == 0 & (x$Nham_nt != 1)),,drop = F]; x <- sum(x$gr_over_sigmasquared, na.rm = TRUE)/sum(x$one_over_sigmasquared, na.rm = TRUE); return(x)})

## Calculate coefficients for linear transformation
scaling <- cbind(data.frame(rbind(do.call(c, stops.block1), do.call(c, syn.block1))),
                 data.frame(rbind(do.call(c, stops.block2), do.call(c, syn.block2))),
                 data.frame(rbind(do.call(c, stops.block3), do.call(c, syn.block3))),
                 data.frame(rbind(do.call(c, stops.block4), do.call(c, syn.block4))))
rownames(scaling) <- c("stops", "syn")
colnames(scaling) <- c(paste0(rep("block1_", 12), 1:12),
                       paste0(rep("block2_", 12), 1:12),
                       paste0(rep("block3_", 12), 1:12),
                       paste0(rep("block4_", 12), 1:12))

## Linear transformation
coefficients <- matrix(NA, ncol = 48, nrow = 2)
colnames(coefficients) <-  colnames(scaling)
for (i in 1:12){
  
  ### Block
  for(j in 1:4){
    
    if(any(colnames(scaling) %in% paste0("block", j, "_", i))){
      
      tmp.lm <- lm(formula = c(0, scaling[2,paste0("block1_", i)]) ~ scaling[,paste0("block", j, "_", i)])
      tmp.lm <- summary(tmp.lm)
      coefficients[1,paste0("block", j, "_", i)] <- tmp.lm$coefficients[[2]]
      coefficients[2,paste0("block", j, "_", i)] <- tmp.lm$coefficients[[1]]  
      
    }
    
  }
  
}
for (i in 1:length(PYL1_block1)){
  
  PYL1_block1[[i]]$gr_normalised <- PYL1_block1[[i]]$growthrate*coefficients[1,i] + coefficients[2,i]
  PYL1_block1[[i]]$gr_sigma_normalised <- PYL1_block1[[i]]$growthrate_sigma*coefficients[1,i]
  PYL1_block2[[i]]$gr_normalised <- PYL1_block2[[i]]$growthrate*coefficients[1,c(i+12)] + coefficients[2,c(i+12)]
  PYL1_block2[[i]]$gr_sigma_normalised <- PYL1_block2[[i]]$growthrate_sigma*coefficients[1,c(i+12)]
  PYL1_block3[[i]]$gr_normalised <- PYL1_block3[[i]]$growthrate*coefficients[1,c(i+24)] + coefficients[2,c(i+24)]
  PYL1_block3[[i]]$gr_sigma_normalised <- PYL1_block3[[i]]$growthrate_sigma*coefficients[1,c(i+24)]
  PYL1_block4[[i]]$gr_normalised <- PYL1_block4[[i]]$growthrate*coefficients[1,c(i+36)] + coefficients[2,c(i+36)]
  PYL1_block4[[i]]$gr_sigma_normalised <- PYL1_block4[[i]]$growthrate_sigma*coefficients[1,c(i+36)]
  
}

## Clean up environment
rm(i, j, tmp.lm)


## 3. Combine PYL1 blocks ##
############################

## Extend the WT aa sequence of each PYL1 block by the other blocks
PYL1_block1 <- lapply(PYL1_block1, function(x) {x$aa_seq <- paste0(x$aa_seq,
                                                                   "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                                                   "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                                                   "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN");
return(x)})
PYL1_block2 <- lapply(PYL1_block2, function(x) {x$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                                                   x$aa_seq,
                                                                   "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                                                   "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN");
return(x)})
PYL1_block3 <- lapply(PYL1_block3, function(x) {x$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                                                   "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                                                   x$aa_seq,
                                                                   "ESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN");
return(x)})
PYL1_block4 <- lapply(PYL1_block4, function(x) {x$aa_seq <- paste0("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRR",
                                                                   "FDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLD",
                                                                   "LLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVL",
                                                                   x$aa_seq);
return(x)})

## Merge PYL1 blocks
PYL1.ABI1 <- vector(mode = "list", length = 12)
names(PYL1.ABI1) <- names(PYL1_block1)
for(i in 1:12){
  
  PYL1.ABI1[[i]] <- rbind(PYL1_block1[[i]],
                          PYL1_block2[[i]],
                          PYL1_block3[[i]],
                          PYL1_block4[[i]])
  
}

## Only keep mutations with AA Hamming distance of 1
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- x[which(x$Nham_aa == 1 | x$Nham_aa == 0),]; return(x)})

## Annotate the individual mutation types
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x <- cbind("Pos" = rep(NA, nrow(x)),
                                                      "WT_AA" = rep(NA, nrow(x)),
                                                      "Mut" = rep(NA, nrow(x)),
                                                      x); return(x)})
for(j in 1:length(PYL1.ABI1)){
  
  for(i in 1:nrow(PYL1.ABI1[[j]])){
    
    ### Do not act on WT seqs
    if(PYL1.ABI1[[j]][i,"Nham_aa"] == 0){
      
      next
      
    }
    
    ### Find the position and mutation
    WT.seq <- PYL1.ABI1[[j]][which(PYL1.ABI1[[j]][,"WT"] == T)[1],"aa_seq"]
    WT.seq <- as.character(str_split_fixed(WT.seq, "", 177))
    mut.seq <- PYL1.ABI1[[j]][i,"aa_seq"]
    mut.seq <- as.character(str_split_fixed(mut.seq, "", 177))
    PYL1.ABI1[[j]][i,"Pos"] <- which(mut.seq != WT.seq)
    PYL1.ABI1[[j]][i,"WT_AA"] <- WT.seq[PYL1.ABI1[[j]][i,"Pos"]]
    PYL1.ABI1[[j]][i,"Mut"] <- mut.seq[PYL1.ABI1[[j]][i,"Pos"]]
    
  } 
  
}

## Positional adjustment with respect to the full-length PYL1 protein
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x$Pos <- x$Pos + 32; return(x)})

## WT specification
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x[which(x[,"Nham_aa"] == 0),"Pos"] <- "WT"; x[which(x[,"Nham_aa"] == 0),"WT_AA"] <- "WT"; x[which(x[,"Nham_aa"] == 0),"Mut"] <- "WT"; return(x)})

## Dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}
names(PYL1.ABI1) <- dosages


## 4. Scale growth rates by wildtype B[inf] ##
##############################################

## Extract WT aa variants (incl. synonymous ones)
PYL1.ABI1.WT <- lapply(PYL1.ABI1, function(x){y <- x[which(x[,"Nham_aa"] == 0),]; return(y)})

## Build a dose-response (and error) matrix for each nucleotide sequence
PYL1.ABI1.WT.mat <- PYL1.ABI1.WT.sigma.mat <- matrix(NA, ncol = 12, nrow = length(unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))))
colnames(PYL1.ABI1.WT.mat) <- colnames(PYL1.ABI1.WT.sigma.mat) <- names(PYL1.ABI1.WT)
rownames(PYL1.ABI1.WT.mat) <- rownames(PYL1.ABI1.WT.sigma.mat) <- unique(do.call(c, sapply(PYL1.ABI1.WT, function(x){x$nt_seq})))
for(i in 1:12){
  PYL1.ABI1.WT.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_normalised"]
  PYL1.ABI1.WT.sigma.mat[,i] <- PYL1.ABI1.WT[[i]][match(rownames(PYL1.ABI1.WT.mat), PYL1.ABI1.WT[[i]]$nt_seq),"gr_sigma_normalised"]
}

## Remove synonymous variants not fully covered across (+)-ABA concentrations
PYL1.ABI1.WT.mat <- na.omit(PYL1.ABI1.WT.mat)
PYL1.ABI1.WT.sigma.mat <- na.omit(PYL1.ABI1.WT.sigma.mat)

## Generate the dose response curve of the (nucleotide-level) WT
WT.PYL1.drc <- cbind(PYL1.ABI1.WT.mat[1,],PYL1.ABI1.WT.sigma.mat[1,],colnames(PYL1.ABI1.WT.mat))
class(WT.PYL1.drc) <- "numeric"
WT.PYL1.drc <- as.data.frame(WT.PYL1.drc)
colnames(WT.PYL1.drc) <- c("GR", "GR sigma", "concentration")

WT.PYL1.drc <- drm(WT.PYL1.drc$GR ~ WT.PYL1.drc$concentration,
                   weights = 1/WT.PYL1.drc$`GR sigma`,
                   fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                   type = 'continuous')
WT.PYL1.drc.par <- WT.PYL1.drc$fit$par
names(WT.PYL1.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.par <- WT.PYL1.drc.par[c(2:4,1)]
WT.PYL1.drc.par[4] <- -WT.PYL1.drc.par[4]

## Scale all growth rates to the wildtype B[inf]
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x$gr_normalised_WTscaled <- 100*x$gr_normalised/WT.PYL1.drc.par["B[inf]"]; return(x)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){x$gr_sigma_normalised_WTscaled <- 100*x$gr_sigma_normalised/WT.PYL1.drc.par["B[inf]"]; return(x)})

## Save PYL1.ABI1 list as an .Rdata file
save(PYL1.ABI1, file = "../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")


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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] reshape_0.8.10 drc_3.0-1      MASS_7.3-65    stringr_1.5.1 
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5        cli_3.6.5          rlang_1.1.6        TH.data_1.1-4      Formula_1.2-5      stringi_1.8.7      car_3.1-3         
# [8] gtools_3.9.5       glue_1.8.0         zoo_1.8-14         plyr_1.8.9         scales_1.4.0       grid_4.5.1         carData_3.0-5     
# [15] abind_1.4-8        mvtnorm_1.3-3      lifecycle_1.0.4    compiler_4.5.1     multcomp_1.4-28    codetools_0.2-20   RColorBrewer_1.1-3
# [22] sandwich_3.1-1     Rcpp_1.1.0         rstudioapi_0.17.1  farver_2.1.2       lattice_0.22-7     R6_2.6.1           splines_4.5.1     
# [29] magrittr_2.0.3     Matrix_1.7-4       tools_4.5.1        plotrix_3.8-4      survival_3.8-3
