# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 19.11.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#######################################################
## Figure 5B - aPCA residual distance to ABA and PPI ##
#######################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel", "cowplot", "rlang")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input aPCA data and bPCA Hill parameters ##
#################################################

load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Mut"] == "*"),]

## input Hill parameter distributions from dose-response curve fits, filtering like in Figure 2
load("../../data/DRCs/PYL1-ABI1_parameters_Hill_12conc.RData")
parameters.Hill <- parameters.Hill[which(parameters.Hill[,"R^2"] > 0.9),]
parameters.Hill <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]


## 2. Combine data ##
#####################

parameters.Hill <- cbind(parameters.Hill, "aPCA" = rep(NA, nrow(parameters.Hill)))
for (i in 1:nrow(parameters.Hill)){
  
  print(i)
  tmp.res <- rownames(parameters.Hill)[i]
  if(tmp.res == "WT"){
    
    parameters.Hill[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[which(PYL1.PYL1.0uM.ABA$WT == T)[1],"gr_normalised_WTscaled"]
    
  }else{
    
    tmp.pos <- as.numeric(substr(tmp.res, 2, nchar(tmp.res)-1))
    tmp.WT <- substr(tmp.res, 1, 1)
    tmp.MUT <- substr(tmp.res, nchar(tmp.res), nchar(tmp.res))
    parameters.Hill[i,"aPCA"] <- PYL1.PYL1.0uM.ABA[which(as.numeric(PYL1.PYL1.0uM.ABA$Pos) == tmp.pos & PYL1.PYL1.0uM.ABA$Mut == tmp.MUT),"gr_normalised_WTscaled"]
    
  }
  
}


## 3. Calculate residuals & significance ##
###########################################

## LOESS of all four parameters vs. pseudo-abundance
loess.B0 <- loess(`B[0]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Binf <- loess(`B[inf]` ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.EC50 <- loess(log(EC50) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.Hill <- loess(log(Hill) ~ aPCA, as.data.frame(parameters.Hill), span = 0.5, family = "symmetric")
loess.res <- as.data.frame(cbind("B0" = loess.B0$residuals,
                                 "Binf" = loess.Binf$residuals,
                                 "EC50" = loess.EC50$residuals,
                                 "Hill" = loess.Hill$residuals))

## summarise by position
loess.res.pos <- vector(mode = "list", length = 4)
names(loess.res.pos) <- colnames(loess.res)
loess.res.pos <- lapply(loess.res.pos, 
                        function(x){x <- matrix(NA, nrow = 176, ncol = 20); 
                        rownames(x) <- unique(substr(rownames(loess.res), 1, nchar(rownames(loess.res))-1))[-1]; 
                        colnames(x) <- c("G", "A", "V", "L", "M", "I", "F", 
                                         "Y", "W", "K", "R", "H", "D", "E", 
                                         "S", "T", "C", "N", "Q", "P"); return(x)})
for(i in 2:nrow(loess.res)){
  print(i)
  tmp.pos <- substr(rownames(loess.res)[i], 1, nchar(rownames(loess.res)[i])-1)
  tmp.mut <- substr(rownames(loess.res)[i], nchar(rownames(loess.res)[i]), nchar(rownames(loess.res)[i]))
  loess.res.pos$B0[tmp.pos,tmp.mut] <- loess.res[i,"B0"]
  loess.res.pos$Binf[tmp.pos,tmp.mut] <- loess.res[i,"Binf"]
  loess.res.pos$EC50[tmp.pos,tmp.mut] <- loess.res[i,"EC50"]
  loess.res.pos$Hill[tmp.pos,tmp.mut] <- loess.res[i,"Hill"]
}

## position-wise significance test
loess.res.pos.sig <- matrix(NA, nrow = 176, ncol = 4)
rownames(loess.res.pos.sig) <- rownames(loess.res.pos$B0)
colnames(loess.res.pos.sig) <- names(loess.res.pos)
for(i in 1:nrow(loess.res.pos$B0)){
  loess.res.pos.sig[i,"B0"] <- wilcox.test(x = loess.res.pos$B0[-i,], 
                                           y = loess.res.pos$B0[i,], 
                                           alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"Binf"] <- wilcox.test(x = loess.res.pos$Binf[-i,], 
                                             y = loess.res.pos$Binf[i,], 
                                             alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"EC50"] <- wilcox.test(x = loess.res.pos$EC50[-i,], 
                                             y = loess.res.pos$EC50[i,], 
                                             alternative = "two.sided")$p.value
  loess.res.pos.sig[i,"Hill"] <- wilcox.test(x = loess.res.pos$Hill[-i,], 
                                             y = loess.res.pos$Hill[i,], 
                                             alternative = "two.sided")$p.value
}
loess.res.pos.sig.fdr <- loess.res.pos.sig
loess.res.pos.sig.fdr[,"B0"] <- p.adjust(loess.res.pos.sig[,"B0"], method = "BH")
loess.res.pos.sig.fdr[,"Binf"] <- p.adjust(loess.res.pos.sig[,"Binf"], method = "BH")
loess.res.pos.sig.fdr[,"EC50"] <- p.adjust(loess.res.pos.sig[,"EC50"], method = "BH")
loess.res.pos.sig.fdr[,"Hill"] <- p.adjust(loess.res.pos.sig[,"Hill"], method = "BH")


## 4. Plot ##
#############

## input residue minimal distances from ABA and ABI1
load("../../data/PYL1_distances_to_interfaces.Rdata")

## summarise data
out <- cbind("(+)-ABA dist [Å]" = PYL1_ABA_ABI1_dist$ABA$ABA_vs_PYL1_HAmin_ligand[-c(1:2)][-76],
             "ABI1 dist [Å]" = PYL1_ABA_ABI1_dist$ABI1$ABI1_vs_PYL1_HAmin_ligand[-c(1:2)][-76],
             "Mean B[0] res." = apply(loess.res.pos$B0, 1, mean, na.rm = T),
             "Mean B[inf] res." = apply(loess.res.pos$Binf, 1, mean, na.rm = T),
             "Mean EC50 res." = apply(loess.res.pos$EC50, 1, mean, na.rm = T),
             "Mean Hill res." = apply(loess.res.pos$Hill, 1, mean, na.rm = T),
             "p.adj B[0]" = loess.res.pos.sig.fdr[,"B0"],
             "p.adj B[inf]" = loess.res.pos.sig.fdr[,"Binf"],
             "p.adj EC50" = loess.res.pos.sig.fdr[,"EC50"],
             "p.adj Hill" = loess.res.pos.sig.fdr[,"Hill"])
out <- as.data.frame(out)
out$names <- rownames(out)

## plot
out.5B.B0 <- ggplot(out, aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`)) +
  geom_point(aes(size = -log10(`p.adj B[0]`), 
                 color = squish(`Mean B[0] res.`, range = c(-30,30))), 
             alpha = 0.8, shape = 16) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-30,30),
                        oob = squish) +
  geom_point(data = out[grep("^T33|^E36|^F37|^Q39|^S43", rownames(out)),],
             aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`,
                 size = -log10(`p.adj B[0]`)), 
             shape = 21, fill = NA , color = "black") +
  scale_size_continuous(name = bquote(italic(P) ~ "(FDR adjusted)"),
                        breaks = c(1, 4, 7, 10),
                        labels = expression(10^-1, 10^-4, 10^-7, 10^-10),
                        range = c(0, 11)) +
  scale_x_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  scale_y_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.title = element_text(size = 30, family = "Helvetica", face = "bold", margin = margin(b = 10, t = 10)),
        legend.text = element_text(size = 25, family = "Helvetica", margin = margin(l = 7)),
        legend.key.size = unit(2, "lines"),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(color = "", 
       title = "",
       x = "(+)-ABA distance [Å]",
       y = "PPI distance [Å]") +
  guides(color = guide_colorbar(title = bquote(B[0] ~ "residual"), order = 1), 
         size = guide_legend(title = bquote(italic(P) ~ "(FDR adjusted)"), order = 2))

out.5B.Binf <- ggplot(out, aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`)) +
  geom_point(aes(size = -log10(`p.adj B[inf]`), 
                 color = squish(`Mean B[inf] res.`, 
                                range = c(-30,30))), 
             alpha = 0.8, stroke = 0, shape = 16) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-30,30),
                        oob = squish) +
  geom_point(data = out[grep("^D129|^R130|^R131|^E160|^E161|^E162|^E163|^R164", rownames(out)),],
             aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`,
                 size = -log10(`p.adj B[inf]`)), 
             shape = 21, fill = NA , color = "black") +
  scale_size_continuous(name = bquote(italic(P) ~ "(FDR adjusted)"),
                        breaks = c(1, 4, 7, 10),
                        labels = expression(10^-1, 10^-4, 10^-7, 10^-10),
                        range = c(0, 11)) +
  scale_x_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  scale_y_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.title = element_text(size = 30, family = "Helvetica", face = "bold", margin = margin(b = 10, t = 10)),
        legend.text = element_text(size = 25, family = "Helvetica", margin = margin(l = 7)),
        legend.key.size = unit(2, "lines"),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(color = "", 
       title = "",
       x = "(+)-ABA distance [Å]",
       y = "PPI distance [Å]") +
  guides(color = guide_colorbar(title = bquote(B[infinity] ~ "residual"), order = 1),
         size = "none")

out.5B.EC50 <- ggplot(out, aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`)) +
  geom_point(aes(size = -log10(`p.adj EC50`), 
                 color = squish(`Mean EC50 res.`, 
                                range = c(-5.5, 5.5))), 
             alpha = 0.8, stroke = 0, shape = 16) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-5.5, 5.5),
                        oob = squish,
                        breaks = log10(c(10^-4, 10^-2, 10^0, 10^2, 10^4)),
                        labels = c(expression(10^-4), 
                                   expression(10^-2), 
                                   expression(10^-0),
                                   expression(10^2), 
                                   expression(10^4))) +
  geom_point(data = out[grep("^Y50|^H87|^E171|^L196|^K200", rownames(out)),],
             aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`,
                 size = -log10(`p.adj EC50`)), 
             shape = 21, fill = NA , color = "black") +
  scale_size_continuous(name = bquote(italic(P) ~ "(FDR adjusted)"),
                        breaks = c(1, 4, 7, 10),
                        labels = expression(10^-1, 10^-4, 10^-7, 10^-10),
                        range = c(0, 11)) +
  scale_x_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  scale_y_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.title = element_text(size = 30, family = "Helvetica", face = "bold", margin = margin(b = 10, t = 10)),
        legend.text = element_text(size = 25, family = "Helvetica", margin = margin(l = 7)),
        legend.key.size = unit(2, "lines"),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(color = "", 
       title = "",
       x = "(+)-ABA distance [Å]",
       y = "PPI distance [Å]") +
  guides(color = guide_colorbar(title = bquote(EC[50] ~ "residual"), order = 1),
         size = "none")

out.5B.n <- ggplot(out, aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`)) +
  geom_point(aes(size = -log10(`p.adj Hill`), 
                 color = squish(`Mean Hill res.`, 
                                range = c(-0.7, 0.7))), 
             alpha = 0.8, stroke = 0, shape = 16) +
  scale_color_gradient2(low = "red", mid = "gray90", high = "blue", 
                        midpoint = 0, 
                        limits = c(-0.7, 0.7),
                        oob = squish,
                        breaks = log10(c(10^-0.6, 10^-0.3, 10^0, 10^0.3, 10^0.6)),
                        labels = c(expression(10^-0.6), 
                                   expression(10^-0.3), 
                                   expression(10^0), 
                                   expression(10^0.3), 
                                   expression(10^0.6))) +
  geom_point(data = out[grep("^S41|^Q42|^S43|^P68|^H154", rownames(out)),],
             aes(x = `(+)-ABA dist [Å]`, y = `ABI1 dist [Å]`,
                 size = -log10(`p.adj Hill`)), 
             shape = 21, fill = NA , color = "black") +
  scale_size_continuous(name = bquote(italic(P) ~ "(FDR adjusted)"),
                        breaks = c(1, 4, 7, 10),
                        labels = expression(10^-1, 10^-4, 10^-7, 10^-10),
                        range = c(0, 11)) +
  scale_x_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  scale_y_continuous(breaks = seq(f = 0, t = 35, by = 5), 
                     limits = c(0, 36)) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 35, vjust = 3),
        legend.title = element_text(size = 30, family = "Helvetica", face = "bold", margin = margin(b = 10, t = 10)),
        legend.text = element_text(size = 25, family = "Helvetica", margin = margin(l = 7)),
        legend.key.size = unit(2, "lines"),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(color = "", 
       title = "",
       x = "(+)-ABA distance [Å]",
       y = "PPI distance [Å]") +
  guides(color = guide_colorbar(title = bquote(Hill ~ "residual"), order = 1), 
         size = "none")

pdf('../../results/Figure5/Figure5B_distances_vs_aPCA_residual.pdf', width = 13, height = 35)
plot_grid(out.5B.B0, out.5B.Binf, out.5B.EC50, out.5B.n, align = "hv", axis = "tblr", ncol = 1)
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
# [1] rlang_1.1.6    cowplot_1.2.0  ggrepel_0.9.6  ggtext_0.1.2   ggplot2_4.0.0  reshape_0.8.10 scales_1.4.0   stringr_1.5.2 
# 
# loaded via a namespace (and not attached):
# [1] gtable_0.3.6       dplyr_1.1.4        compiler_4.5.1     crayon_1.5.3       tidyselect_1.2.1   Rcpp_1.1.0         xml2_1.4.0         R6_2.6.1           plyr_1.8.9         commonmark_2.0.0  
# [11] labeling_0.4.3     generics_0.1.4     tibble_3.3.0       pillar_1.11.1      RColorBrewer_1.1-3 stringi_1.8.7      litedown_0.7       xfun_0.53          S7_0.2.0           cli_3.6.5         
# [21] withr_3.0.2        magrittr_2.0.4     grid_4.5.1         gridtext_0.1.5     markdown_2.0       lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2       tools_4.5.1       
# [31] pkgconfig_2.0.3 