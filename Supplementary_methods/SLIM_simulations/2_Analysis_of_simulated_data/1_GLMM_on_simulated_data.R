#!/usr/bin/env Rscript

### ### ### ### ### ### ### ### libraries ### ### ### ### ### ### ### ### ### 

suppressMessages(library("glmmTMB"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))


### ### ### ### ### ### ### ### aesthetics ### ### ### ### ### ### ### ### ### 

effect_labels <- c(
  "b_AgrVSNonAgrNon-Agriculture"= "Agr. vs. non-agriculture",
  "AgrVSNonAgrNon-Agriculture"= "Agr. vs. non-agriculture",
  "b_AgrVSNonAgrNonMAgriculture"= "Agr. vs. non-agriculture")

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
Path_to_GLMM_results <- ".../project_amylase/simulations/03_results_from_simulation/GLMM/"

SIMULATION_RUNS <- list.files(Path_to_AMY1_CN_info, pattern = "\\.txt$", full.names = FALSE)
SIMULATION_RUNS <- SIMULATION_RUNS[file.info(file.path(Path_to_AMY1_CN_info, SIMULATION_RUNS))$size > 0]
SIMULATION_RUNS <- sub("\\.txt$", "", SIMULATION_RUNS)


### ### ### ### ### ### ### ### loop over simulation runs ### ### ### ### ### ### ### 

for (SIMULATION_RUN in SIMULATION_RUNS) {

  ### VARIABLES
  Path_to_AMY1_CN_info=paste0(".../project_amylase/simulations/01_simulated_VCFs_and_CN/CN_files/", SIMULATION_RUN, ".txt", sep="")

  Path_to_PCA_info=paste0(".../project_amylase/simulations/02_simulated_MDIST_PCA/PCA/", SIMULATION_RUN, sep="")

  Path_to_GLMM_results=paste0(".../project_amylase/simulations/03_results_from_simulation/GLMM/glmm_", SIMULATION_RUN, sep="")

  ## CONSTANT
  Path_to_pop_metadata=".../project_amylase/simulations/00_simulated_pop_names_and_metadata/POPULATION_METADATA.txt"


  ### ### ###  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  ### POP METADATA
  GLMM_pop_metadata <- read.table(Path_to_pop_metadata, header = TRUE, sep = ";", strip.white = TRUE, stringsAsFactors = FALSE)

  ### AMY1 INFO
  GLMM_CN_and_pop_info <- read.table(Path_to_AMY1_CN_info, header=FALSE, fill=TRUE)
  colnames(GLMM_CN_and_pop_info) <- c('sample', 'SLIM_ID', 'AMY1_CN', 'pop')

  ### (MERGED)
  GLMM_data <- merge(GLMM_CN_and_pop_info, GLMM_pop_metadata, by = "pop", all.x = TRUE)

  ### PCA INFO
  GLMM_PCAeigenvec <- read.table(paste0(Path_to_PCA_info, ".eigenvec", sep=""), header = FALSE, fill = TRUE)
  colnames(GLMM_PCAeigenvec) <- c("sample", "sample2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

  GLMM_PCAeigenval <- read.table(paste0(Path_to_PCA_info, ".eigenval", sep=""), header = FALSE)
  GLMM_PCAeigenval$PC <- as.numeric(rownames(GLMM_PCAeigenval))
  colnames(GLMM_PCAeigenval) <- c("Eigenvalue", "PC")

  ### (FINAL DATAFRAME)
  GLMM <- merge(GLMM_data, GLMM_PCAeigenvec, by = "sample", all = FALSE)


  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  #### MODELING SPECIFICATIONS

  base_vars <- paste0("PC", 1:5)
  formula <- paste("AMY1_CN ~", paste(base_vars, collapse = " + "), "+ (1 | pop)")
  base_formula <- as.formula(formula)
  family='compois'

  #### RUN MODELS
  
  glmm_null <- glmmTMB(base_formula, data=GLMM, family=family)
  glmm_Agr.vs.NonAgr <- glmmTMB(update(base_formula, . ~ . + AgrVSNonAgr), data=GLMM, family=family)

  models <- list(Null=glmm_null, Agr.vs.NonAgr=glmm_Agr.vs.NonAgr)
  model_names <- names(models)

  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

  #### SAVE GLMM MODELING RESULTS

  model_summary <- list()

  for (name in model_names) {
    model <- models[[name]]
    
    x <- as.data.frame(summary(model)$coefficients$cond)
    names(x)[names(x) == "Pr(>|z|)"] <- "p_value"
    ci <- as.data.frame(confint(model))
    x$Variable <- rownames(x)
    x <- x[!x$Variable %in% c("(Intercept)", "Std.Dev.(Intercept)|pop"), ]
    
    x$CI_lower <- ci[x$Variable, "2.5 %"]
    x$CI_upper <- ci[x$Variable, "97.5 %"]
    x$model <- name
    x$simulation_run <- SIMULATION_RUN
    
    model_summary[[name]] <- x
  }

  model_summary <- do.call(rbind, model_summary)
  model_summary <- model_summary[, c("model", "Variable", "simulation_run", setdiff(names(model_summary), c("model", "Variable", "simulation_run")))]

  model_summary$Label <- effect_labels[model_summary$Variable]
  model_summary$Label[is.na(model_summary$Label)] <- 
    model_summary$Variable[is.na(model_summary$Label)]
  }


  write.table(model_summary,paste0(Path_to_GLMM_results, ".txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)

}

