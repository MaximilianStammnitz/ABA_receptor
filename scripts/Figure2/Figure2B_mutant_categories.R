# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################
## Figure 2B - Mutant category splits ##
########################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "ggplot2")

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

## export into Supplementary Table 2
## PYL1.summary.nonWT <- PYL1.summary.nonWT[,12:1]
## write.table(PYL1.summary.nonWT, file = "../../data/Supplementary Tables/TableS2A.txt", sep = "\t", quote = F, col.names = T, row.names = T)

## remove any mutants with incomplete curves
PYL1.summary.nonWT <- PYL1.summary.nonWT[-sort(unique(do.call(c,apply(PYL1.summary.nonWT[,1:12], 2, function(x){which(is.na(x) == T)})))),1:12]
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)


## 2. Doughnut distribution of mutants ##
#########################################

## input Hill parameters
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")

## remove mutants with incomplete curve data
parameters.Hill <- parameters.Hill[which(rownames(parameters.Hill) %in% rownames(PYL1.summary.nonWT)),]

## minimal binders
min.binders <- parameters.Hill[which(parameters.Hill[,"B[inf]"] < 25 | parameters.Hill[,"EC50"] > 1000 | is.na(parameters.Hill[,"B[0]"])),]
min.binders <- min.binders[order(min.binders[,"B[0]"]),]

## hypersensitive binders
hyper.binders <- parameters.Hill[which(parameters.Hill[,"B[inf]"] >= 50 & parameters.Hill[,"B[0]"] >= 50),]
hyper.binders <- hyper.binders[which(hyper.binders[,"R^2"] > 0.85),]
hyper.binders <- hyper.binders[order(hyper.binders[,"B[0]"]),]

## inverters
inv.binders <- parameters.Hill[which(c(parameters.Hill[,"B[0]"] - parameters.Hill[,"B[inf]"]) >= 20 & parameters.Hill[,"B[0]"] > 30),]
inv.binders <- inv.binders[which(inv.binders[,"R^2"] > 0.85),]
inv.binders <- inv.binders[order(inv.binders[,"B[0]"]),]

## band-stops
detect_bandstop <- function(y, min_points = 2, min_step = 5) {
  
  # y: numeric vector of responses (ordered by dose)
  # min_points: minimum consecutive points for each down/up segment
  # min_step: minimum absolute change between consecutive points
  # gap_tolerance: max number of intermediate steps allowed between down and up
  
  if(y[1] < 30 & y[2] < 30){
    
    return(FALSE)
    
  }
  
  ## look at the point-to-point changes
  dy <- diff(y)
  
  ## only keep steps > min_step
  
  ## are there three ups and downs of length >= min_points
  dy.abs <- dy
  dy.abs[which(dy < 0)] <- "down"
  dy.abs[which(dy >= 0)] <- "up"
  dy.abs[which(abs(dy) < min_step)] <- "noise"
  r <- rle(dy.abs)
  
  ## go over each sequence of increases/decreases
  if(length(r$lengths) > 1){
    
    for(j in 1:c(length(r$lengths) - 1)){
      
      if(r$values[j] == "down" & r$lengths[j] >= min_points){
        
        if(r$values[j + 1] == "up" & r$lengths[j + 1] >= min_points){
          
          return(TRUE)
          
        }else if(!is.na(r$values[j + 2])){
          
          if(r$values[j + 2] == "up" & r$lengths[j + 2] >= min_points){
            
            return(TRUE)
            
          }else{
            
            next
            
          }
          
        }else{
          
          next
          
        }
        
      }else{
        
        next
        
      }
      
    } 
    
    return(FALSE)
    
  }else{
    
    return(FALSE)
    
  }
  
  
}
PYL1.summary.nonWT.nostop <- PYL1.summary.nonWT[-grep("[*]", rownames(PYL1.summary.nonWT)),]
bandstop.binders.id <- c()
for(i in 1:nrow(PYL1.summary.nonWT.nostop)){
  
  bandstop.tmp <- detect_bandstop(y = as.numeric(rev(PYL1.summary.nonWT.nostop[i,])),
                                  min_points = 2, 
                                  min_step = 5)
  if(bandstop.tmp == T){
    bandstop.binders.id <- c(bandstop.binders.id, rownames(PYL1.summary.nonWT.nostop)[i])
  }else{
    next
  }
  
}

### manual checks of the same positions (remove, add)
bandstop.binders.id <- bandstop.binders.id[-grep("G113P|A116P|S119C|R131H|H142G|H142E|
                                                  |P178V|N181A|N181I|N181H|
                                                  |N181T|N181Q|N181E|T186I|
                                                  |A190H|V193L",bandstop.binders.id)]
bandstop.binders.id <- c(bandstop.binders.id,
                         "S112M", "L114G", "S119V", "S119I", "H142A", "H142C", "H142Q", "R145P", "V193W", "V193H")

## summarise
inv.binders.id <- rownames(inv.binders)
inv.binders.id <- c(inv.binders.id, "V110F")
hyper.binders.id <- rownames(hyper.binders)
hyper.binders.id <- hyper.binders.id[-which(hyper.binders.id %in% inv.binders.id | hyper.binders.id %in% bandstop.binders.id)]
min.binders.id <- rownames(min.binders)
min.binders.id <-  min.binders.id[-which(min.binders.id %in% hyper.binders.id | min.binders.id %in% inv.binders.id | min.binders.id %in% bandstop.binders.id)]

## add all the "normal" DRCs
parameters.Hill <- parameters.Hill[order(parameters.Hill[,"B[0]"]),]
others.id <- rownames(parameters.Hill)
others.id <-  others.id[-which(others.id %in% min.binders.id | others.id %in% hyper.binders.id | others.id %in% inv.binders.id | others.id %in% bandstop.binders.id)]

## summary
mut.clust <- c(length(others.id),
               length(min.binders.id),
               length(hyper.binders.id),
               length(bandstop.binders.id),
               length(inv.binders.id))

# Prepare the data frame
doughnut.in <- data.frame(category = paste0(c("Hill", "Minimal", "Hypersensitive", "Bandstop", "Inverter"), " (", 
                                            round(100*mut.clust/sum(mut.clust), 2), "%)"),
                          value = as.numeric(sort(mut.clust, decreasing = T)))
doughnut.in$category <- factor(doughnut.in$category, 
                               levels = rev(doughnut.in$category))

# Create a basic doughnut plot
custom_colors <- c("grey90", "grey70", "grey50", "grey30", "grey10")

pdf('../../results/Figure2/Figure2B_mutant_distribution_doughnut.pdf', width = 10, height = 10)

ggplot(doughnut.in, aes(x = 2, y = value, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_polar(theta = "y") +
  xlim(0.1, 2.9) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")

dev.off()


## 3. Version ##
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
# [1] ggplot2_4.0.0 scales_1.4.0  stringr_1.5.2
# 
# loaded via a namespace (and not attached):
# [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.4     gtable_0.3.6       glue_1.8.0         tibble_3.3.0       pkgconfig_2.0.3   
# [11] generics_0.1.4     dplyr_1.1.4        lifecycle_1.0.4    cli_3.6.5          S7_0.2.0           grid_4.5.1         vctrs_0.6.5        withr_3.0.2        compiler_4.5.1     tools_4.5.1       
# [21] pillar_1.11.1      rlang_1.1.6        stringi_1.8.7 