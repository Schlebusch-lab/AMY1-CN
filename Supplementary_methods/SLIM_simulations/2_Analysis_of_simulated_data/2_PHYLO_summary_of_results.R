#!/usr/bin/env Rscript

### ### ### ### ### ### ### ### libraries ### ### ### ### ### ### ### ### ### 

suppressMessages(library("ggplot2"))
suppressMessages(library("patchwork"))


#### READ ALL MODEL RESULTS

Path_to_PHYLO_results <- ".../project_amylase/simulations/03_results_from_simulation/phylosignal/"

files_phyloCorr <- list.files(
  Path_to_PHYLO_results,
  pattern = "^phylo_.*\\_phyloCorr\\.txt$",
  full.names = TRUE)

files_phyloSignal <- list.files(
  Path_to_PHYLO_results,
  pattern = "^phylo_.*\\_phyloSignal\\.txt$",
  full.names = TRUE)

summary_phyloCorr <- dplyr::bind_rows(lapply(files_phyloCorr, function(file) {
  x <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  x$file <- basename(file)
  x$s_value <- sub(".*_(S[0-9]{3})\\_phyloCorr\\.txt$", "\\1", basename(file))
  x}))

summary_phyloSignal <- dplyr::bind_rows(lapply(files_phyloSignal, function(file) {
  x <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  x$file <- basename(file)
  x$s_value <- sub(".*_(S[0-9]{3})\\_phyloSignal\\.txt$", "\\1", basename(file))
  x}))

rownames(summary_phyloCorr) <- NULL
rownames(summary_phyloSignal) <- NULL


##### PHYLOGENETIC SIGNAL

summary_phyloSignal$significant <- summary_phyloSignal$p_value < 0.05

ggplot(summary_phyloSignal, aes(x = significant, y = Statistic_value, color = significant)) +
  geom_jitter(width = 0.3, alpha = 0.6, size = 1.5) +
  facet_grid(Statistic ~ s_value, scales = "free_y", 
    labeller = labeller(s_value = function(x)
      paste0("s=", sprintf("%.3f", as.numeric(sub("^S", "", x)) / 1000)))) +
  scale_x_discrete(labels = c("FALSE" = "p ≥ 0.05", "TRUE" = "p < 0.05")) +
  scale_color_manual(values = c("TRUE" = "cornflowerblue", "FALSE" = "black"), labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"), name='Significance') +
  labs(x=NULL, y = "Statistic value") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


##### CORRELOGRAM

run_colors <- setNames(rainbow(length(unique(summary_phyloCorr$SIMULATION_RUN))), sort(unique(summary_phyloCorr$SIMULATION_RUN)))

ggplot(summary_phyloCorr, aes(x = distance, y = estimate, group = SIMULATION_RUN, color = SIMULATION_RUN, fill = SIMULATION_RUN)) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.15, linewidth = 0) +
  geom_line(alpha = 0.8, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ s_value, nrow=3,
    labeller = labeller(s_value = function(x)
      paste0("s=", sprintf("%.3f", as.numeric(sub("^S", "", x)) / 1000)))) +
  scale_color_manual(values = run_colors) +
  scale_fill_manual(values = run_colors) +
  labs(x = "Phylogenetic distance",  y = "Correlation", color = "Simulation run", fill = "Simulation run") +
  theme_bw() + theme(legend.position = "none")

