# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#######################################################
## Supplementary Figure 1 - GluePCA TECAN benchmarks ##
#######################################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "growthrates", "drc", "scales")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Import TECAN data ##
##########################

## PYL1-ABI1
pyl1.abi1 <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/PYL1-ABI1/GluePCA_benchmark_PYL1_ABI1_TECAN.xlsx', sheet = 1))
rownames(pyl1.abi1) <- pyl1.abi1[,1]
pyl1.abi1 <- pyl1.abi1[,-c(1:2)]
pyl1.abi1 <- pyl1.abi1[-1,]
colnames(pyl1.abi1) <- as.numeric(colnames(pyl1.abi1))/3600 ## convert to hours
class(pyl1.abi1) <- "numeric"
pyl1.abi1 <- as.data.frame(pyl1.abi1)

## GID1A-GAI ##
gid1a.gai <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/GID1A-GAI/GluePCA_benchmark_GID1A-GAI_TECAN.xlsx', sheet = 1))
rownames(gid1a.gai) <- gid1a.gai[,1]
gid1a.gai <- gid1a.gai[,-c(1:2)]
gid1a.gai <- gid1a.gai[-1,]
colnames(gid1a.gai) <- as.numeric(colnames(gid1a.gai))/3600 ## convert to hours
class(gid1a.gai) <- "numeric"
gid1a.gai <- as.data.frame(gid1a.gai)

## FKBP12-FRB ##
fkbp12.frb <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/FKBP12-FRB/GluePCA_benchmark_FKBP12_FRB_TECAN.xlsx', sheet = 1))
rownames(fkbp12.frb) <- fkbp12.frb[,1]
fkbp12.frb <- fkbp12.frb[,-c(1:2)]
fkbp12.frb <- fkbp12.frb[-1,]
colnames(fkbp12.frb) <- as.numeric(colnames(fkbp12.frb))/3600 ## convert to hours
class(fkbp12.frb) <- "numeric"
fkbp12.frb <- as.data.frame(fkbp12.frb)


## 2. Plot ##
#############

## PYL1-ABI1
pdf("../../results/FigureS1/FigureS1A_GluePCA_benchmark_PYL1-ABI1.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without ABA
plot(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.6), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep1",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep2",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep3",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep4",])

#### with ABA
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.6, by = 0.1), labels = seq(f = 0.1, t = 0.6, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("380 µM (+)-ABA", "no ABA"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### PYL1-ABI1 (ABI1 W300A)

#### without ABA
plot(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.6), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", main = "ABI1 W300A", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep4",], col = "black")

#### with ABA
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.6, by = 0.1), labels = seq(f = 0.1, t = 0.6, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## GID1A-GAI ##
pdf("../../results/FigureS1/FigureS1B_GluePCA_benchmark_GID1A-GAI.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without GA3
plot(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.15,1.2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep1",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep2",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep3",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep4",])

#### with GA3
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep3",],pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.2, t = 1.2, by = 0.2), labels = seq(f = 0.2, t = 1.2, by = 0.2), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("1 mM GA3", "no GA3"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### GID1A-GAI (GAI E54R)

#### without GA3
plot(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.15,1.2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", 
     main = "GAI E54R", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep4",], col = "black")

#### with GA3
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.2, t = 1.2, by = 0.2), labels = seq(f = 0.2, t = 1.2, by = 0.2), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## FKBP12-FRB ##
pdf("../../results/FigureS1/FigureS1C_GluePCA_benchmark_FRB-FKBP12.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without Rapamycin
plot(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.5), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep1",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep2",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep3",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep4",])

#### with Rapamycin
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uµM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.5, by = 0.1), labels = seq(f = 0.1, t = 0.5, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("10 µM Rapamycin", "no Rapamycin"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### FKBP12-FRB (FRB PLF + Y2088A)

#### without Rapamycin
plot(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.5), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", main = "FRB PLF Y2088A", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep4",], col = "black")

#### with Rapamycin
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.5, by = 0.1), labels = seq(f = 0.1, t = 0.5, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## 3. Version ##
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
# [1] scales_1.3.0      drc_3.0-1         MASS_7.3-64       growthrates_0.8.4 deSolve_1.40      lattice_0.22-6    readxl_1.4.3     
# 
# loaded via a namespace (and not attached):
# [1] Matrix_1.7-2      compiler_4.4.1    gtools_3.9.5      plotrix_3.8-4     Rcpp_1.0.14       minpack.lm_1.2-4  parallel_4.4.1   
# [8] splines_4.4.1     coda_0.19-4.1     TH.data_1.1-3     R6_2.6.1          Formula_1.2-5     tibble_3.2.1      car_3.1-3        
# [15] munsell_0.5.1     minqa_1.2.8       pillar_1.10.1     rlang_1.1.5       multcomp_1.4-28   cli_3.6.4         FME_1.3.6.3      
# [22] magrittr_2.0.3    grid_4.4.1        rstudioapi_0.17.1 mvtnorm_1.3-3     rootSolve_1.8.2.4 sandwich_3.1-1    lifecycle_1.0.4  
# [29] vctrs_0.6.5       glue_1.8.0        cellranger_1.1.0  codetools_0.2-20  zoo_1.8-12        survival_3.8-3    abind_1.4-8      
# [36] carData_3.0-5     colorspace_2.1-1  pkgconfig_2.0.3   tools_4.4.1 