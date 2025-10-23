# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 23.10.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################################
## Supplementary Figure 2D & 2E - PYL1 fitness distributions ##
###############################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "ggplot2", "ggtext", "cowplot")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed_v2.RData")


## 2. Stratify mutants ##
#########################

## Build matrices
PYL1.ABI1.gr <- vector(mode = "list", length = 12)
names(PYL1.ABI1.gr) <- names(PYL1.ABI1)
PYL1.ABI1.gr <- lapply(PYL1.ABI1.gr, function(x){x <- matrix(NA, nrow = 177, ncol = 21); 
rownames(x) <- 1:177; colnames(x) <- c("G", "A", "V", "L", "M", "I", "F", 
                                       "Y", "W", "K", "R", "H", "D", "E", 
                                       "S", "T", "C", "N", "Q", "P","*"); 
return(x)})

## fill-in fitness matrix & rename position IDs by the actual positions in PYL1
for(j in 1:length(PYL1.ABI1)){
  
  for(i in 1:nrow(PYL1.ABI1[[j]])){
    
    if(PYL1.ABI1[[j]][i,"Pos"] == "WT"){
      
      next
      
    }else{
      
      PYL1.ABI1.gr[[j]][as.character(c(as.numeric(PYL1.ABI1[[j]][i,"Pos"]) - 32)),PYL1.ABI1[[j]][i,"Mut"]] <- as.numeric(PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"])
      
    }
    
  }
  
}
WT.seq <- c(str_split_fixed("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "", 177))
PYL1.ABI1.gr <- lapply(PYL1.ABI1.gr, function(x){rownames(x) <- paste0(WT.seq,33:209); class(x) <- "numeric"; return(x)})

## key annotations
### 1. ABA contacts vs. others
### 2. ABI1 contacts vs. others
### 3. alpha vs. beta vs. loops
alpha <- beta <- loop <- core <- surface <- ABA_cont <- ABI1_cont <- rep(0, 177)
names(alpha) <- names(beta) <- names(loop) <- names(core) <- names(surface) <- names(ABA_cont) <- names(ABI1_cont) <-  33:209
ABI1_cont[c(c(87,88,90,111:117,142:144,178,180,181,185,186,188,189,192,193,196,199)-32)] <- 1
ABA_cont[c(c(86,88,110,116,121,135,137,143,144,147,149,171,189,193,197)-32)] <- 1
loop[c(c(48:54,65:68,77:81,85:90,95:104,111:118,128:131,138:147,158:164,177:183)-32)] <- 1
beta[c(c(55:64,91:94,105:110,119:127,132:137,148:157,165:176)-32)] <- 1
alpha[c(c(33:47,69:76,82:84,184:208)-32)] <- 1

### 4. core vs. surface
SASA <- read.table("../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/SASA_PyMol.txt", header = F)[,1:2]
SASA[,2] <- as.numeric(gsub("%", "", SASA[,2]))/100
SASA <- SASA[-c(1:2),]
SASA <- SASA[grep("/A/", SASA[,1]),]
core[which(SASA[,"V2"] <= 0.25)] <- 1
surface[which(SASA[,"V2"] > 0.25)] <- 1

### 5. Automised thresholding of data into "stop-like" vs. "WT-like" (unfortunately, Gaussian mixture models fail at low conc.)
bin.split <- vector(mode = "list", length = 12)
for(i in 1:12){
  
  # Kernel densities, equivalent to Figure 1G
  tmp.dist <- density(PYL1.ABI1[[i]][,"gr_normalised_WTscaled"], 
                      n = 1000, kernel = 'gaussian', 
                      from = -20, to = 130)
  tmp.dy <- diff(tmp.dist$y)
  
  # Find indices where derivative changes sign from positive to negative (local max)
  tmp.peak_indices <- which(diff(sign(tmp.dy)) == -2) + 1  # +1 to correct index shift
  tmp.peak_x <- tmp.dist$x[tmp.peak_indices]  # X-coordinates of peaks
  tmp.peak_y <- tmp.dist$y[tmp.peak_indices]  # Y-coordinates of peaks
  
  # Indices between the two main peaks
  tmp.peak_indices <- sort(tmp.peak_indices[order(tmp.peak_y, decreasing = T)[1:2]])
  tmp.min_x <- tmp.dist$x[tmp.peak_indices[1] + which.min(tmp.dist$y[tmp.peak_indices[1]:tmp.peak_indices[2]]) - 1]
  
  # Plot, sanity check of separations
  plot(tmp.dist)
  abline(v = sort(tmp.peak_x[order(tmp.peak_y, decreasing = T)[1:2]], decreasing = T), col = "red")
  abline(v = tmp.min_x, col = "blue")
  
  # With this threshold, how many missense variants are above/below
  tmp.total.missense <- length(which(PYL1.ABI1[[i]][,"Nham_aa"] == 1 & PYL1.ABI1[[i]][,"STOP"] != T))
  tmp.total.missense.WTlike <- length(which(PYL1.ABI1[[i]][which(PYL1.ABI1[[i]][,"Nham_aa"] == 1 & PYL1.ABI1[[i]][,"STOP"] != T),"gr_normalised_WTscaled"] >= tmp.min_x))
  
  # Output
  bin.split[[i]] <- tmp.min_x
  
  # Clean-up
  rm(tmp.min_x, tmp.total.missense, tmp.total.missense.WTlike, tmp.peak_indices, 
     tmp.peak_y, tmp.peak_x, tmp.dy, tmp.dist)
  
  ## check if distribution is not unimodal
  #tmp.dip <- dip.test(PYL1.ABI1[[i]][,"gr_normalised_WTscaled"])
  # if(i < 8){
  #   
  #   ### run GMM
  #   bin.split[[i]] <- normalmixEM(x = PYL1.ABI1[[i]][,"gr_normalised_WTscaled"],
  #                                 mu = c(mean(PYL1.ABI1[[i]][which(PYL1.ABI1[[i]][,"STOP"] == T),"gr_normalised_WTscaled"]),
  #                                        mean(PYL1.ABI1[[i]][which(PYL1.ABI1[[i]][,"Nham_aa"] == 0),"gr_normalised_WTscaled"])),
  #                                 sigma = c(sd(PYL1.ABI1[[i]][which(PYL1.ABI1[[i]][,"STOP"] == T),"gr_normalised_WTscaled"]),
  #                                           sd(PYL1.ABI1[[i]][which(PYL1.ABI1[[i]][,"Nham_aa"] == 0),"gr_normalised_WTscaled"])),
  #                                 k = 2,
  #                                 maxrestarts = 1000)
  #   
  #   ### depict threshold GR
  #   lower <- max(bin.split[[i]]$x[which(bin.split[[i]]$posterior[,1] > 0.9)])
  #   upper <- min(bin.split[[i]]$x[which(bin.split[[i]]$posterior[,2] > 0.9)])
  #   bin.split[[i]] <- c(lower, upper)
  #   
  # }
  
}


## 3. Plot raw binding distributions at 0, 16.7 and 2,500 µM (+)-ABA ##
#######################################################################

## Data separation at 0 µM (+)-ABA
tmp.dist <- density(PYL1.ABI1[[12]]$gr_normalised_WTscaled, 
                    n = 10000, kernel = 'gaussian', 
                    from = -1000, to = 1000)
plot_data <- data.frame(x = tmp.dist$x, y = tmp.dist$y,
                        dosage = as.numeric(names(PYL1.ABI1))[12])
plot_data$category <- ifelse(plot_data$x >= bin.split[[12]], "Above", "Below")
plot_data$y[1] <- 0
plot_data$y[c(which(plot_data$category == "Above")[1] - 1)] <- 0
plot_data$y[which(plot_data$category == "Above")[1]] <- 0
plot_data$y[1000] <- 0

out.S2D_low <- ggplot(plot_data, aes(x = x, y = y, fill = category)) +
  geom_polygon(aes(group = category),
               color = "black",
               linewidth = 0.5, na.rm=T) +
  geom_vline(xintercept = bin.split[[12]], color = "black", size = 2) + 
  annotate("text", x = bin.split[[12]] - 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`0` < bin.split[[12]])) / length(PYL1.ABI1.gr$`0`[which(is.na(PYL1.ABI1.gr$`0`) == F)]), 1), "%"),
           size = 7, color = "black") + 
  annotate("text", x = bin.split[[12]] + 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`0` >= bin.split[[12]])) / length(PYL1.ABI1.gr$`0`[which(is.na(PYL1.ABI1.gr$`0`) == F)]), 1), "%"),
           size = 7, color = "darkgreen") +
  scale_fill_manual(values = c("Above" = "darkgreen", "Below" = "white"),
                    name = "Density Threshold") +
  scale_x_continuous(limits = c(-1000, 1000),
                     breaks = seq(0, 100, length.out = 6),
                     labels = seq(0, 100, length.out = 6)) +
  coord_cartesian(xlim = c(-35, 135)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(size = 40, hjust = 0.5),
        plot.subtitle = element_markdown(size = 20),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(title = "0 µM (+)-ABA",
       x = "Binding", y = "")

pdf("../../results/FigureS2/FigureS2D_variant_distr_low_v2.pdf", width = 6, height = 10)
print(out.S2D_low)
dev.off()

## Data separation at 16.7 µM (+)-ABA
tmp.dist <- density(PYL1.ABI1[[5]]$gr_normalised_WTscaled, 
                    n = 10000, kernel = 'gaussian', 
                    from = -1000, to = 1000)
plot_data <- data.frame(x = tmp.dist$x, y = tmp.dist$y,
                        dosage = as.numeric(names(PYL1.ABI1))[1])
plot_data$category <- ifelse(plot_data$x >= bin.split[[5]], "Above", "Below")
plot_data$y[1] <- 0
plot_data$y[c(which(plot_data$category == "Above")[1] - 1)] <- 0
plot_data$y[which(plot_data$category == "Above")[1]] <- 0
plot_data$y[1000] <- 0

out.S2D_mid <- ggplot(plot_data, aes(x = x, y = y, fill = category)) +
  geom_polygon(aes(group = category),
               color = "black",
               linewidth = 0.5, na.rm=T) +
  geom_vline(xintercept = bin.split[[5]], color = "black", size = 2) +  # Add vertical line
  annotate("text", x = bin.split[[5]] - 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`16.6597251145356` < bin.split[[5]])) / length(PYL1.ABI1.gr$`16.6597251145356`[which(is.na(PYL1.ABI1.gr$`16.6597251145356`) == F)]), 1), "%"),
           size = 7, color = "black") + 
  annotate("text", x = bin.split[[5]] + 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`16.6597251145356` >= bin.split[[5]])) / length(PYL1.ABI1.gr$`16.6597251145356`[which(is.na(PYL1.ABI1.gr$`16.6597251145356`) == F)]), 1), "%"),
           size = 7, color = "darkgreen") +
  scale_fill_manual(values = c("Above" = "darkgreen", "Below" = "white"),
                    name = "Density Threshold") +
  scale_x_continuous(limits = c(-1000, 1000),
                     breaks = seq(0, 100, length.out = 6),
                     labels = seq(0, 100, length.out = 6)) +
  coord_cartesian(xlim = c(-35, 135)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(size = 40, hjust = 0.5),
        plot.subtitle = element_markdown(size = 20),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(title = "16.7 µM (+)-ABA",
       x = "Binding", y = "")

pdf("../../results/FigureS2/FigureS2D_variant_distr_mid_v2.pdf", width = 6, height = 10)
print(out.S2D_mid)
dev.off()

## Data separation at 2,500 µM (+)-ABA
tmp.dist <- density(PYL1.ABI1[[1]]$gr_normalised_WTscaled, 
                    n = 10000, kernel = 'gaussian', 
                    from = -1000, to = 1000)
plot_data <- data.frame(x = tmp.dist$x, y = tmp.dist$y,
                        dosage = as.numeric(names(PYL1.ABI1))[1])
plot_data$category <- ifelse(plot_data$x >= bin.split[[1]], "Above", "Below")
plot_data$y[1] <- 0
plot_data$y[c(which(plot_data$category == "Above")[1] - 1)] <- 0
plot_data$y[which(plot_data$category == "Above")[1]] <- 0
plot_data$y[1000] <- 0

out.S2D_high <- ggplot(plot_data, aes(x = x, y = y, fill = category)) +
  geom_polygon(aes(group = category),
               color = "black",
               linewidth = 0.5, na.rm=T) +
  geom_vline(xintercept = bin.split[[1]], color = "black", size = 2) +  # Add vertical line
  annotate("text", x = bin.split[[1]] - 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`2500` < bin.split[[1]])) / length(PYL1.ABI1.gr$`2500`[which(is.na(PYL1.ABI1.gr$`2500`) == F)]), 1), "%"),
           size = 7, color = "black") + 
  annotate("text", x = bin.split[[1]] + 30, 
           y = max(plot_data$y) * 1.1, 
           label = paste0(round(100 * length(which(PYL1.ABI1.gr$`2500` >= bin.split[[1]])) / length(PYL1.ABI1.gr$`2500`[which(is.na(PYL1.ABI1.gr$`2500`) == F)]), 1), "%"),
           size = 7, color = "darkgreen") +
  scale_fill_manual(values = c("Above" = "darkgreen", "Below" = "white"),
                    name = "Density Threshold") +
  scale_x_continuous(limits = c(-1000, 1000),
                     breaks = seq(0, 100, length.out = 6),
                     labels = seq(0, 100, length.out = 6)) +
  coord_cartesian(xlim = c(-35, 135)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(size = 40, hjust = 0.5),
        plot.subtitle = element_markdown(size = 20),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(title = "2,500 µM (+)-ABA",
       x = "Binding", y = "")

pdf("../../results/FigureS2/FigureS2D_variant_distr_high_v2.pdf", width = 6, height = 10)
print(out.S2D_high)
dev.off()


## 4. For 16.7 µM (+)-ABA: how do core and binding residues perform ##
######################################################################

i <- 5

## combine data and annotations
df <- PYL1.ABI1[[i]]

### 1. stop-like vs. WT-like split
df$impact <- rep(NA,nrow(df))
df$impact[which(df$gr_normalised_WTscaled < bin.split[[i]])] <- "stop-like"
df$impact[which(df$gr_normalised_WTscaled >= bin.split[[i]])] <- "WT-like"
df$impact <- factor(df$impact, levels = c("stop-like", "WT-like"))

### 2. mutation types
df$type <- rep(NA,nrow(df))
df$type[which(df$Nham_aa == 1 & df$STOP == F)] <- "missense"
df$type[which(df$Nham_aa == 1 & df$STOP == T)] <- "stops"
df$type[which(df$Nham_aa == 0 & is.na(df$WT) == T)] <- "synonymous"
df$type <- factor(df$type, levels = c("missense", "stops", "synonymous"))

### 3. ABA contacts vs. others
df$ABAcontact <- rep(NA,nrow(df))
df$ABAcontact[which(df$Pos %in% as.numeric(names(which(ABA_cont == 1))))] <- "ABA contact"
df$ABAcontact[which(!df$Pos %in% as.numeric(names(which(ABA_cont == 1))))] <- "other"
df$ABAcontact <- factor(df$ABAcontact, levels = c("ABA contact", "other"))

### 4. ABI1 contacts vs. others
df$ABI1contact <- rep(NA,nrow(df))
df$ABI1contact[which(df$Pos %in% as.numeric(names(which(ABI1_cont == 1))))] <- "ABI1 contact"
df$ABI1contact[which(!df$Pos %in% as.numeric(names(which(ABI1_cont == 1))))] <- "other"
df$ABI1contact <- factor(df$ABI1contact, levels = c("ABI1 contact", "other"))

### 5. alpha vs. beta vs. loops
df$secondary <- rep(NA,nrow(df))
df$secondary[which(df$Pos %in% as.numeric(names(which(alpha == 1))))] <- "alpha"
df$secondary[which(df$Pos %in% as.numeric(names(which(beta == 1))))] <- "beta"
df$secondary[which(df$Pos %in% as.numeric(names(which(loop == 1))))] <- "loop"
df$secondary <- factor(df$secondary, levels = c("alpha", "beta", "loop"))

### 6. core vs. surface
df$region <- rep(NA,nrow(df))
df$region[which(df$Pos %in% as.numeric(names(which(core == 1))))] <- "core"
df$region[which(df$Pos %in% as.numeric(names(which(surface == 1))))] <- "surface"
df$region <- factor(df$region, levels = c("core", "surface"))

### 7. clean up data frame
df <- df[-which(is.na(df$type) == T),]

## now let's make contigency tables

### 1. core vs. surface (missense only)
cont.1 <- rbind(c(length(which(df$type == "missense" & df$region == "core" & df$impact == "stop-like")), length(which(df$type == "missense" & df$region == "surface" & df$impact == "stop-like"))),
                c(length(which(df$type == "missense" & df$region == "core" & df$impact == "WT-like")), length(which(df$type == "missense" & df$region == "surface" & df$impact == "WT-like"))))
fisher.test(cont.1)

### 2. ABA vs. no ABA contacts (missense only)
cont.2 <- rbind(c(length(which(df$type == "missense" & df$ABAcontact == "ABA contact" & df$impact == "stop-like")), length(which(df$type == "missense" & df$ABAcontact == "other" & df$impact == "stop-like"))),
                c(length(which(df$type == "missense" & df$ABAcontact == "ABA contact" & df$impact == "WT-like")), length(which(df$type == "missense" & df$ABAcontact == "other" & df$impact == "WT-like"))))
fisher.test(cont.2)

### 3. ABI1 vs. no ABI1 contacts (missense only)
cont.3 <- rbind(c(length(which(df$type == "missense" & df$ABI1contact == "ABI1 contact" & df$impact == "stop-like")), length(which(df$type == "missense" & df$ABI1contact == "other" & df$impact == "stop-like"))),
                c(length(which(df$type == "missense" & df$ABI1contact == "ABI1 contact" & df$impact == "WT-like")), length(which(df$type == "missense" & df$ABI1contact == "other" & df$impact == "WT-like"))))
fisher.test(cont.3)


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
# [1] cowplot_1.2.0 ggtext_0.1.2  ggplot2_4.0.0 scales_1.4.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3       vctrs_0.6.5        cli_3.6.5          xfun_0.53          rlang_1.1.6       
# [6] stringi_1.8.7      generics_0.1.4     S7_0.2.0           litedown_0.7       labeling_0.4.3    
# [11] glue_1.8.0         markdown_2.0       gridtext_0.1.5     grid_4.5.1         tibble_3.3.0      
# [16] lifecycle_1.0.4    compiler_4.5.1     dplyr_1.1.4        RColorBrewer_1.1-3 Rcpp_1.1.0        
# [21] pkgconfig_2.0.3    rstudioapi_0.17.1  farver_2.1.2       R6_2.6.1           tidyselect_1.2.1  
# [26] pillar_1.11.0      commonmark_2.0.0   magrittr_2.0.3     tools_4.5.1        withr_3.0.2       
# [31] gtable_0.3.6       xml2_1.4.0  
