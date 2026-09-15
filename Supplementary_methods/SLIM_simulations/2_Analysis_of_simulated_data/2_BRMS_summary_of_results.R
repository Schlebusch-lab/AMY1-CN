#!/usr/bin/env Rscript

### ### ### ### ### ### ### ### libraries ### ### ### ### ### ### ### ### ### 

suppressMessages(library("ggplot2"))
suppressMessages(library("patchwork"))
suppressMessages(library("dplyr"))
suppressMessages(library(tidyr))

### ### ### ### ### ### ### ### aesthetics ### ### ### ### ### ### ### ### ### 

model_palette <- c(Null = "orange2", Agr.vs.NonAgr = "red3")

design <- "
    AAAAACCCCC
    AAAAACCCCC
    AAAAACCCCC
    AAAAACCCCC
    AAAAACCCCC
    AAAAACCCCC
    AAAAACCCCC
    BBBBBCCCCC
    BBBBBCCCCC
    BBBBBCCCCC
  "


#### READ ALL MODEL SUMMARIES

BRMS_results_path <- ".../simulations/03_results_from_simulation/brms/"

files <- list.files(
  BRMS_results_path,
  pattern = "^brms_.*\\.txt$",
  full.names = TRUE)

model_summary_all <- dplyr::bind_rows(lapply(files, function(file) {

  x <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  x$file <- basename(file)
  x$s_value <- sub(".*_(S[0-9]{3})\\.txt$", "\\1", basename(file))

  x
}))

rownames(model_summary_all) <- NULL


#### SIMULATION SUMMARY

significance_summary <- model_summary_all %>%

  group_by(s_value, model, Label) %>% summarise(
    n = sum(!is.na(CIs_include_zero)),

    n_CI_excludes_zero = sum(CIs_include_zero == "NO", na.rm = TRUE),
    percent_CI_excludes_zero = round(100 * sum(CIs_include_zero == "NO", na.rm = TRUE) / n, 2),

    n_negative_significant = sum(CIs_include_zero == "NO" & Estimate < 0, na.rm = TRUE),
    percent_CI_excludes_zero_and_negative = round(100 * n_negative_significant / n, 2),
    
    percent_CI_excludes_zero_and_positive = round(percent_CI_excludes_zero - percent_CI_excludes_zero_and_negative, 2),

    .groups = "drop")


#### PLOT

plot_dat <- bind_rows(

  significance_summary %>% filter(Label == "Agr. vs. non-agriculture") %>%
    transmute(s_value, model, Label,
    direction = "Significant (pos)",
    percent = percent_CI_excludes_zero_and_positive),

  significance_summary %>% filter(Label == "Agr. vs. non-agriculture") %>%
    transmute(s_value, model, Label,
     direction = "Significant (neg)",
     percent = percent_CI_excludes_zero_and_negative),

  significance_summary %>% filter(Label != "Agr. vs. non-agriculture") %>%
    transmute(s_value, model, Label,
    direction = "Significant",
    percent = percent_CI_excludes_zero))

ggplot(plot_dat, aes(x = Label, y = percent, fill = direction)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_grid(model ~ s_value, scales = "free_x",
    labeller = labeller(s_value = function(x)
      paste0("s=", sprintf("%.3f", as.numeric(sub("^S", "", x)) / 1000)))) +
  scale_fill_manual(values = c("Significant (pos)" = "#ddd79f", "Significant (neg)" = "#bbaa16", "Significant" = "#6a620b"), name='Direction') +
  labs(x = "", y = "Proportion of replicates (%) with significant results") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))









