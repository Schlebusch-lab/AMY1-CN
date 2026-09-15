#!/usr/bin/env Rscript

### ### ### ### ### ### ### ### libraries ### ### ### ### ### ### ### ### ### 

suppressMessages(library('brms'))
suppressMessages(library('ape'))
suppressMessages(library('zCompositions'))
suppressMessages(library('compositions'))
suppressMessages(library('readr'))
suppressMessages(library('ggplot2'))
suppressMessages(library('tidybayes'))
suppressMessages(library('ggridges'))
suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))
suppressMessages(library('gridGraphics'))
suppressMessages(library('grid'))
suppressMessages(library('phytools'))

### ### ### ### ### ### ### ### aesthetics ### ### ### ### ### ### ### ### ### 

effect_labels <- c(
  "b_AgrVSNonAgrNon-Agriculture"= "Agr. vs. non-agriculture",
  "AgrVSNonAgrNon-Agriculture"= "Agr. vs. non-agriculture",
  "b_AgrVSNonAgrNonMAgriculture"= "Agr. vs. non-agriculture",
  "sd_pop__Intercept" = "Population identity (sd)", 
  "sd_sample__Intercept" = "Individual-level VCV (sd)")

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
Path_to_BRMS_results <- ".../project_amylase/simulations/03_results_from_simulation/brms/"

SIMULATION_RUNS <- list.files(Path_to_AMY1_CN_info, pattern = "\\.txt$", full.names = FALSE)
SIMULATION_RUNS <- SIMULATION_RUNS[file.info(file.path(Path_to_AMY1_CN_info, SIMULATION_RUNS))$size > 0]
SIMULATION_RUNS <- sub("\\.txt$", "", SIMULATION_RUNS)

### ### ### ### ### ### ### ### loop over simulation runs ### ### ### ### ### ### ### 

for (SIMULATION_RUN in SIMULATION_RUNS) {

  ### VARIABLES

  Path_to_AMY1_CN_info=paste0(".../project_amylase/simulations/01_simulated_VCFs_and_CN/CN_files/", SIMULATION_RUN, ".txt", sep="")

  Path_to_MDIST_info=paste0(".../project_amylase/simulations/02_simulated_MDIST_PCA/MDIST/", SIMULATION_RUN, sep="")

  Path_to_BRMS_results=paste0(".../project_amylase/simulations/03_results_from_simulation/brms/brms_", SIMULATION_RUN, sep="")

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


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  #### MODELING SPECIFICATIONS

  chains=4
  cores=12
  iter=12000
  warmup=4000
  adapt_delta=0.99 
  max_treedepth=10 
  family=poisson()
  seed=1234

  formula <- paste("AMY1_CN ~","(1|gr(sample,cov=A)) + (1|pop)")
  base_formula <- as.formula(formula)

  priors_all_default_except_sd_sample <- set_prior("normal(0, 0.5)", class = "sd", group = "sample")  

  A <- vcv.phylo(tree_chron)

  #### MODELS

  brms_null <- brm(
    base_formula,
    data = brms_data, data2 = list(A=A),
    prior=priors_all_default_except_sd_sample,
    save_pars = save_pars(all = TRUE),
    family = family, sample_prior=TRUE, seed=seed, chains=chains, cores=cores, iter=iter, warmup=warmup,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    silent=2, refresh=0)

  brms_Agr.vs.NonAgr <- brm(
    update(base_formula, . ~ . + AgrVSNonAgr),
    data = brms_data, data2 = list(A=A),
    prior=priors_all_default_except_sd_sample,
    save_pars = save_pars(all = TRUE),
    family = family, sample_prior=TRUE, seed=seed, chains=chains, cores=cores, iter=iter, warmup=warmup,
    control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth),
    silent=2, refresh=0)

  models <- list(Null=brms_null, Agr.vs.NonAgr=brms_Agr.vs.NonAgr)
  model_names <- names(models)

  cat('******* Models run for simulation: ', SIMULATION_RUN, '\n')


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


  ### SUMMARY OF MODEL RESULTS

  model_summary <- list()
  posterior_plots <- list()

  for (name in model_names) {
    
    model <- models[[name]]
    draws <- as.data.frame(as_draws_df(model))
    
    vars <- c(grep("^b_", names(draws), value = TRUE), "sd_pop__Intercept", "sd_sample__Intercept")
    vars <- vars[vars != "b_Intercept"]
    
    draws <- draws[, vars, drop = FALSE]
    
    ### SUMMARY
    x <- data.frame(
      Variable = names(draws),
      Estimate = sapply(draws, mean),
      Est.Error = sapply(draws, sd),
      CI_lower = sapply(draws, quantile, probs = 0.025),
      CI_upper = sapply(draws, quantile, probs = 0.975))
    
    x$model <- name
    x$simulation_run <- SIMULATION_RUN
    x$CIs_include_zero <- ifelse(x$CI_lower <= 0 & x$CI_upper >= 0, "YES", "NO")

    model_summary[[name]] <- x

  }


  ### SAVE SUMMARY IN BETTER FORMAT

  model_summary <- do.call(rbind, model_summary)

  model_summary <- model_summary[, c("model", "Variable", "simulation_run", setdiff(names(model_summary), c("model", "Variable", "simulation_run")))]
  model_summary$Label <- effect_labels[model_summary$Variable]
  model_summary$Label[is.na(model_summary$Label)] <- model_summary$Variable[is.na(model_summary$Label)]

  write.table(model_summary, paste0(Path_to_BRMS_results, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)

}

