#!/usr/bin/env Rscript

### ### ### ### ### ### ### ### libraries ### ### ### ### ### ### ### ### ### 

suppressMessages(library('ape'))
suppressMessages(library('readr'))
suppressMessages(library('ggplot2'))
suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))
suppressMessages(library('gridGraphics'))
suppressMessages(library('grid'))
suppressMessages(library(phylosignal))
suppressMessages(library(phylobase))
suppressMessages(library(caper))
suppressMessages(library(phytools))

### ### ### ### ### ### ### ### aesthetics ### ### ### ### ### ### ### ### ### 

pop_cols <- c(
    "pop0" = "#d61ee3",
    "pop1" = "#cf8cd4",
    "pop2" = "#aa7fba",
    "pop3" = "#6e4080",
    "pop4" = "#cb52fa",
    "pop5" = "#1b966f",
    "pop6" = "#3477c9",
    "pop7" = "#94e366",
    "pop8" = "#944521",
    "pop9" = "#6dd6de",
    "pop10" = "#c9be44",
    "pop11" = "#b30e60",
    "pop12" = "#26f0b7")

### ### ### ### ### ### ### ### ### SIMULATION RUNS ### ### ### ### ### ### ### ### 

Path_to_AMY1_CN_info <- ".../project_amylase/simulations/01_simulated_VCFs_and_CN/CN_files/"
Path_to_PHYLO_results <- ".../project_amylase/simulations/03_results_from_simulation/phylosignal/"

SIMULATION_RUNS <- list.files(Path_to_AMY1_CN_info, pattern = "\\.txt$", full.names = FALSE)
SIMULATION_RUNS <- SIMULATION_RUNS[file.info(file.path(Path_to_AMY1_CN_info, SIMULATION_RUNS))$size > 0]
SIMULATION_RUNS <- sub("\\.txt$", "", SIMULATION_RUNS)


### ### ### ### ### ### ### ### loop over simulation runs ### ### ### ### ### ### ### 

for (SIMULATION_RUN in SIMULATION_RUNS) {

  ### VARIABLES

  Path_to_AMY1_CN_info=paste0(".../project_amylase/simulations/01_simulated_VCFs_and_CN/CN_files/", SIMULATION_RUN, ".txt", sep="")

  Path_to_MDIST_info=paste0(".../project_amylase/simulations/02_simulated_MDIST_PCA/MDIST/", SIMULATION_RUN, sep="")

  Path_to_PHYLO_results=paste0(".../project_amylase/simulations/03_results_from_simulation/phylosignal/phylo_", SIMULATION_RUN, sep="")

  ## CONSTANT
  Path_to_pop_metadata=".../project_amylase/simulations/00_simulated_pop_names_and_metadata/POPULATION_METADATA.txt"


  ### ### ###  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  ### POP METADATA
  BRMS_pop_metadata <- read.table(Path_to_pop_metadata, header = TRUE, sep = ";", strip.white = TRUE, stringsAsFactors = FALSE)

  ### AMY1 INFO
  BRMS_CN_and_pop_info <- read.table(Path_to_AMY1_CN_info, header=FALSE, fill=TRUE)
  colnames(BRMS_CN_and_pop_info) <- c('sample', 'SLIM_ID', 'AMY1_CN', 'pop')

  ### (MERGED)
  brms_data <- merge(BRMS_CN_and_pop_info, BRMS_pop_metadata, by = "pop", all.x = TRUE)

  ### MDIST INFO
  IBS_mdist <- as.matrix(read_table(paste0(Path_to_MDIST_info, ".mdist"), col_names = FALSE, show_col_types = FALSE))
  IBS_mdist.id <- read_table(paste0(Path_to_MDIST_info, ".mdist.id"), col_names = FALSE, show_col_types = FALSE)

  rownames(IBS_mdist) <- IBS_mdist.id$X2
  colnames(IBS_mdist) <- IBS_mdist.id$X2
  keep <- rownames(IBS_mdist) %in% brms_data$sample
  IBS_mdist <- IBS_mdist[keep, keep]
  modeling_data <- brms_data[match(rownames(IBS_mdist), brms_data$sample),]
  modeling_data <- modeling_data[!is.na(modeling_data$sample), ]


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  
  # Make tree
  tree <- nj(as.dist(IBS_mdist))

  # Root tree by outgroup
  tip_pop <- brms_data[match(tree$tip.label, brms_data$sample), "pop"]
  outgroup <- tree$tip.label[tip_pop %in% c("pop0","pop1")]
  outgroup_mrca <- getMRCA(tree, outgroup)
  edge_index <- which(tree$edge[, 2] == outgroup_mrca) ## find egde that leads to MRCA pop0-pop1

  tree <- reroot(tree, node = outgroup_mrca, position = tree$edge.length[edge_index] / 2) 

  # Make (slightly) negative branch lengths zero
  if (any(tree$edge.length < 0, na.rm = TRUE)) {
    neg_lengths <- tree$edge.length[tree$edge.length < 0]
    cat("******* Negative branch lengths found.\n******* Summary before correction:\n")
    print(summary(neg_lengths))
    cat("******* Setting negative edge lengths to very small positive number.\n")
    tree$edge.length[tree$edge.length < 0] <- 1e-8}

  # Ultrametricize tree
  tree_chron <- chronos(tree, model="correlated", control = chronos.control(iter.max=10000))
  tree_chron <- structure(unclass(tree_chron), class = "phylo")

  # Plot tree quickly: 'phylogram', 'cladogram', 'fan', 'unrooted', 'radial', 'tidy'
  tip_pop <- brms_data$pop[match(tree_chron$tip.label, brms_data$sample)]
  tip_colors <- pop_cols[as.character(tip_pop)]

  pdf(paste0(Path_to_PHYLO_results, '_phyloTree_REROOT.pdf'), width=7, height=12)
  plot(tree_chron, tip.color=tip_colors, cex=0.3, 'phylogram', show.node.label=TRUE, use.edge.length=TRUE, node.depth=2, align.tip.label=TRUE)
  dev.off()

  # Make phy object
  phy <- comparative.data(phy=tree_chron, data=modeling_data, names.col="sample", vcv=TRUE, na.omit=FALSE)
  phy$data$sample <- phy$phy$tip.label   # Create column with sample names
  phy$data <- droplevels(phy$data)       # Drop levels that are not used

  phy$phy$edge.length[phy$phy$edge.length < 0 & abs(phy$phy$edge.length) < 1e-12] <- 1e-8 ### fix neg branch lengths again (sometimes they re-appear after chrono)

    
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  ### PHYLOGENETIC SIGNAL
  signal_results <- phyloSignal(phylo4d(phy$phy, phy$data["AMY1_CN"]), reps=1000, W=NULL, methods="all")

  signal_summary <- data.frame(SIMULATION_RUN = SIMULATION_RUN, 
                                Statistic = names(signal_results$stat), 
                                Statistic_value = as.numeric(signal_results$stat), 
                                p_value = as.numeric(signal_results$pvalue))

  write.table(signal_summary, paste0(Path_to_PHYLO_results, "_phyloSignal.txt"), sep = "\t", row.names = FALSE, quote = FALSE )


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  ### CORRELOGRAMS
  AMY1.correlogram <- phyloCorrelogram(phylo4d(phy$phy, phy$data["AMY1_CN"]), trait = "AMY1_CN")

  correlogram_summary <- data.frame(
      SIMULATION_RUN = SIMULATION_RUN,
      distance = AMY1.correlogram$res[, 1],
      CI_lower = AMY1.correlogram$res[, 2],
      CI_upper = AMY1.correlogram$res[, 3],
      estimate = AMY1.correlogram$res[, 4])

  write.table(correlogram_summary, paste0(Path_to_PHYLO_results, "_phyloCorr.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

}

