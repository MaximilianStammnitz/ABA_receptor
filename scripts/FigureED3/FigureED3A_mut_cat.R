# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##################################################################
## Extended Data Figure 3A - manual clustering of variant types ##
##################################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "pheatmap")

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

## remove any mutants with incomplete curves
PYL1.summary.nonWT <- PYL1.summary.nonWT[-sort(unique(do.call(c,apply(PYL1.summary.nonWT[,1:12], 2, function(x){which(is.na(x) == T)})))),1:12]
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)


## 2. Heatmap and (manual) clustering stats ##
##############################################

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
hyper.binders <- hyper.binders[order(hyper.binders[,"EC50"], decreasing = T),]

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
parameters.Hill <- parameters.Hill[order(parameters.Hill[,"EC50"], decreasing = T),]
others.id <- rownames(parameters.Hill)
others.id <-  others.id[-which(others.id %in% min.binders.id | others.id %in% hyper.binders.id | others.id %in% inv.binders.id | others.id %in% bandstop.binders.id)]

### heatmap
ids.all <- c(min.binders.id, others.id, hyper.binders.id, inv.binders.id, bandstop.binders.id)
pheatmap(t(PYL1.summary.nonWT[match(ids.all, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40, 
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_mutant_clust.pdf")

pheatmap(t(PYL1.summary.nonWT[match(min.binders.id, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40*length(min.binders.id)/length(ids.all), 
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_mutant_minimal.pdf")

pheatmap(t(PYL1.summary.nonWT[match(others.id, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40*length(others.id)/length(ids.all), 
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_mutant_others.pdf")

pheatmap(t(PYL1.summary.nonWT[match(hyper.binders.id, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40*length(hyper.binders.id)/length(ids.all),
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_hyper.pdf")

pheatmap(t(PYL1.summary.nonWT[match(inv.binders.id, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40*length(inv.binders.id)/length(ids.all),
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_inverters.pdf")

pheatmap(t(PYL1.summary.nonWT[match(bandstop.binders.id, rownames(PYL1.summary.nonWT)),]),
         color = colorRampPalette(c("white", "darkgreen"))(n = 1000),
         breaks = seq(f = 0, to = 120, length.out = 1000),
         show_rownames = F,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         legend = F,
         na_col = "grey",
         display_numbers = F,
         border_color = NA, 
         fontsize_row = 25, 
         fontsize_col = 3,
         treeheight_row = 0, 
         treeheight_col = 100,
         fontsize = 20, 
         height = 11, 
         width = 40*length(bandstop.binders.id)/length(ids.all), 
         cellwidth = 0.78,
         main = "", 
         lwd = 0.5,
         filename = "../../results/FigureED3/FigureED3A_bandstops.pdf")


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
# [1] pheatmap_1.0.13 scales_1.4.0    stringr_1.5.2  
# 
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1     R6_2.6.1           magrittr_2.0.4     cli_3.6.5          tools_4.5.1        gtable_0.3.6       RColorBrewer_1.1-3 glue_1.8.0         farver_2.1.2       vctrs_0.6.5       
# [11] stringi_1.8.7      grid_4.5.1         lifecycle_1.0.4    rlang_1.1.6   