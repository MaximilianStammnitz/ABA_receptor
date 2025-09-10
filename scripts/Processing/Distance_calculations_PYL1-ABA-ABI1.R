# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 10.09.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

################################################
## Pre-processing of PYL1-(+)ABA-ABI1 complex ##
################################################


## 0. Environment ##
####################

## Libraries
packages <- c("bio3d", "data.table")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Atomic Distance Calculation ##
####################################

## Calculation of minimum PYL1 residue distance to both (+)-ABA and ABI1
PYL1_minimum_ABA_ABI1_distances_from_PDB <- function(input_file, PYL1_chain_query, ABI1_chain_query, ABA_query){

  # load PDB structure
  PYL1.ABI1.pdb <- read.pdb(input_file, rm.alt = TRUE)

  ## Atom selections

  # PYL1 atoms
  sele_PYL1 <- atom.select(PYL1.ABI1.pdb,
                           string = "protein",
                           chain = PYL1_chain_query,
                           verbose = FALSE)

  # PYL1 side-chain atoms
  sele_PYL1_sc <- atom.select(PYL1.ABI1.pdb,
                              string = "sidechain",
                              chain = PYL1_chain_query,
                              verbose=FALSE)

  # PYL1 C-alpha atoms
  sele_PYL1_ca <- atom.select(PYL1.ABI1.pdb,
                              string = "calpha",
                              chain = PYL1_chain_query,
                              verbose=FALSE)

  # PYL1 glycine c-alpha atoms
  sele_PYL1_glyca <- atom.select(PYL1.ABI1.pdb,
                                 resid = "GLY",
                                 string = "calpha",
                                 chain = PYL1_chain_query,
                                 verbose = FALSE)

  # Water atoms
  sele_water <- atom.select(PYL1.ABI1.pdb,
                            string = "water",
                            verbose = FALSE)

  # ABA atoms
  sele_ABA <- atom.select(PYL1.ABI1.pdb,
                          resid = ABA_query,
                          verbose = FALSE)

  # ABI1 atoms
  sele_ABI1 <- atom.select(PYL1.ABI1.pdb,
                           string = "protein",
                           chain = ABI1_chain_query,
                           verbose = FALSE)


  ## Combine atom selections ##
  #############################

  # Heavy atoms
  sele_PYL1_HA <- combine.select(sele_PYL1, sele_water, operator = "-", verbose=FALSE)

  # Side chain heavy atoms + c-alpha for glycine
  sele_PYL1_prot_sc <- combine.select(sele_PYL1, sele_PYL1_sc, operator = "AND", verbose=FALSE)
  sele_PYL1_prot_sc_glyca <- combine.select(sele_PYL1_prot_sc, sele_PYL1_glyca, operator = "OR", verbose=FALSE)
  sele_PYL1_scHA <- combine.select(sele_PYL1_prot_sc_glyca, sele_water, operator = "-", verbose=FALSE)

  # heavy PYL1 atoms + ABA/ABI1
  sele_ABA_PYL1_HA <- combine.select(sele_PYL1_HA, sele_ABA, operator = "OR",verbose = FALSE)
  sele_ABI1_PYL1_HA <- combine.select(sele_PYL1_HA, sele_ABI1, operator = "OR",verbose = FALSE)

  # side chain PYL1 heavy atoms + c-alpha for glycine + ABA/ABI1
  sele_ABA_PYL1_scHA <- combine.select(sele_PYL1_scHA, sele_ABA, operator = "OR",verbose = FALSE)
  sele_ABI1_PYL1_scHA <- combine.select(sele_PYL1_scHA, sele_ABI1, operator = "OR",verbose = FALSE)

  # Lists
  sele_list_ABA <- list("ABA_vs_PYL1_HA" = sele_ABA_PYL1_HA, "ABA_vs_PYL1_scHA" = sele_ABA_PYL1_scHA)
  sele_list_ABI1 <- list("ABI1_vs_PYL1_HA" = sele_ABI1_PYL1_HA, "ABI1_vs_PYL1_scHA" = sele_ABI1_PYL1_scHA)


  ## Calculate minimum target chain distances to ABA ##
  #####################################################

  result_ABA_dt <- data.table()
  for(metric in names(sele_list_ABA)){

    # Distance matrix
    PYL1.ABI1.pdb_sub <- trim.pdb(PYL1.ABI1.pdb, sele_list_ABA[[metric]])
    dist_mat <- dm.xyz(PYL1.ABI1.pdb_sub$xyz, grpby = apply(PYL1.ABI1.pdb_sub$atom[,c("resno", "chain")], 1, paste, collapse = "_"), scut = 0, mask.lower = FALSE)
    resno_sub <- unique(PYL1.ABI1.pdb_sub$atom[,c("resno","chain","resid")])

    # Ligand distance matrix
    ABA_dist <- dist_mat[resno_sub[,"chain"] == PYL1_chain_query &! resno_sub[,"resid"] %in% ABA_query, resno_sub[,"resid"] %in% ABA_query]

    # Absolute residue number
    ABA_dist_dt <- data.table(Pos = resno_sub[resno_sub[,"chain"] == PYL1_chain_query &! resno_sub[,"resid"] %in% c(ABA_query),"resno"])

    # Minimum ABA distance
    if(!is.null(dim(ABA_dist))){
      ABA_dist_dt[, min_dist := apply(ABA_dist, 1, min)]
    }else{
      ABA_dist_dt[, min_dist := ABA_dist]
    }

    names(ABA_dist_dt)[2] <- paste0(metric, "min_ligand")
    if(nrow(result_ABA_dt)==0){
      result_ABA_dt <- ABA_dist_dt
    }else{
      result_ABA_dt <- merge(result_ABA_dt, ABA_dist_dt, by = "Pos", all = T)
    }
  }


  ## Calculate minimum target chain distances to ABI1 ##
  ######################################################

  result_ABI1_dt <- data.table()
  for(metric in names(sele_list_ABI1)){

    # Distance matrix
    PYL1.ABI1.pdb_sub <- trim.pdb(PYL1.ABI1.pdb, sele_list_ABI1[[metric]])
    dist_mat <- dm.xyz(PYL1.ABI1.pdb_sub$xyz, grpby = apply(PYL1.ABI1.pdb_sub$atom[,c("resno", "chain")], 1, paste, collapse = "_"), scut = 0, mask.lower = FALSE)
    resno_sub <- unique(PYL1.ABI1.pdb_sub$atom[,c("resno","chain","resid")])

    # Ligand distance matrix
    ABI1_dist <- dist_mat[resno_sub[,"chain"] == PYL1_chain_query &! resno_sub[,"chain"] %in% ABI1_chain_query, resno_sub[,"chain"] %in% ABI1_chain_query]

    # Absolute residue number
    ABI1_dist_dt <- data.table(Pos = resno_sub[resno_sub[,"chain"] == PYL1_chain_query &! resno_sub[,"chain"] %in% c(ABI1_chain_query), "resno"])

    # Minimum ABA distance
    if(!is.null(dim(ABI1_dist))){
      ABI1_dist_dt[, min_dist := apply(ABI1_dist, 1, min)]
    }else{
      ABI1_dist_dt[, min_dist := ABI1_dist]
    }

    names(ABI1_dist_dt)[2] <- paste0(metric, "min_ligand")
    if(nrow(result_ABI1_dt)==0){
      result_ABI1_dt <- ABI1_dist_dt
    }else{
      result_ABI1_dt <- merge(result_ABI1_dt, ABI1_dist_dt, by = "Pos", all = T)
    }
  }

  # Return
  out <- list("ABA" = result_ABA_dt,
              "ABI1" = result_ABI1_dt)
  return(out)

}
PYL1_ABA_ABI1_dist <- PYL1_minimum_ABA_ABI1_distances_from_PDB(input_file = "pdb3kdj.ent",
                                                               PYL1_chain_query = "A",
                                                               ABI1_chain_query = "B",
                                                               ABA_query = "A8S")

## Save list as an .Rdata file
save(PYL1_ABA_ABI1_dist, file = "../../data/PYL1_distances_to_interfaces.Rdata")


## 2. Version ##
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
# loaded via a namespace (and not attached):
# [1] compiler_4.5.1    tools_4.5.1       rstudioapi_0.17.1
