#' ---
#' title: Figure supplements
#' date: ""
#' author: ""
#' output:
#'    rmarkdown::pdf_document:
#'        number_sections: FALSE
#'        toc: FALSE
#' header-includes:
#' - \usepackage{caption}
#' - \captionsetup[figure]{labelformat=empty}
#' bibliography: main/references.bib
#' ---

#+ include=FALSE
knitr::opts_chunk$set(echo = FALSE, message = FALSE, comment = FALSE, warning = FALSE)
if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open="wt")
  sink(log, type = "message")
}


#+ libs
library(stats) # masks filter
library(readr)
library(dplyr)
library(ggplot2)
library(extrafont)
library(stringr)
library(cowplot)
library(patchwork)
library(brms)
library(rstan)
library(tidybayes)
library(here)
library(glue)
old <- theme_set(theme_cowplot(font_size = 12, font_family = "Helvetica"))

#'
#+ params
is_ci <- function() {
  "CI" %in% Sys.getenv()
}
chains <- ifelse(is_ci(), 1, 4)
cores <- chains
refresh = 0
rstan_options(auto_write = TRUE, javascript = FALSE)
if (!dir.exists("results/models")) {
    message("Creating results/models dir..")
    dir.create("results/models", recursive = TRUE)
}

#+ data
pvalues <- read_csv(here("results/pvalues.csv")) %>% 
  rename(de_tool = analysis_platform)
pvalues_sample <- read_csv(here("results/pvalues_sample.csv")) %>% 
  rename(de_tool = analysis_platform)
sequencing_metadata <- read_csv(here("results/sequencing_metadata_unique_platform.csv"))
de_simulation_results <- read_csv(here("results/data/de_simulation_results.csv"))


#+ fig2supp1
qc_threshold <- function(n_pvals, bins, fdr) {
  qbinom(1 - 1 / bins * fdr, n_pvals, 1 / bins)
}

qc_thr <- qc_threshold(n_pvals = 20000, bins = 40, fdr = 0.05)

p <- de_simulation_results %>% 
  mutate(
    reps = factor(reps, levels = c("N=3", "N=6", "N=10")),
    n_eff = factor(n_eff, levels = c("Effects=100", "Effects=200", "Effects=400", "Effects=800"))
  ) %>% 
  ggplot() +
  geom_histogram(aes(pvalue), binwidth = 1/40, center = 1/80) +
  geom_hline(yintercept = qc_thr, color = "red") +
  facet_grid(reps~n_eff) +
  scale_x_continuous(breaks = c(0, 0.5, 1))
ggsave(here("figures/figure_2_figure_supplement_1.tiff"), plot = p, width = 12, height = 8, units = "cm", dpi = 300)

fig_cap <- glue('Figure 2--figure supplement 1. Simulated RNA-seq data shows that histograms from p value sets with around one hundred 
                 true effects out of 20,000 features can be classified as "uniform".
                RNA-seq data was simulated with polyester R package [@frazee2015] on 20,000 transcripts from human transcriptome 
                 using grid of 3, 6, and 10 replicates and 100, 200, 400, and 800 effects for two groups. 
                 Fold changes were set to 0.5 and 2.
                Differential expression was assessed using DESeq2 R package [@love2014] using default settings 
                and group 1 versus group 2 contrast. 
                Effects denotes in facet labels the number of true effects and N denotes number of replicates.
                Red line denotes QC threshold used for dividing p histograms into discrete classes.
                Code and workflow used to run these simulations is available on Github: https://github.com/rstats-tartu/simulate-rnaseq.
                Raw data of the figure is available on Zenodo https://zenodo.org with doi: 10.5281/zenodo.4463803.')

#+ fig.cap=fig_cap
p

#+ FigS1
f <- anticons ~ year
family <- bernoulli()
data <- pvalues_sample
conditions <- NULL
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_year.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions),
          plot = FALSE)

fig_cap <- glue("Figure 3--figure supplement 1. The increasing proportion of anti-conservative histograms. 
                Binomial logistic model: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. 
                Lines denote best fit of linear model. Shaded area denotes 95% credible region.")

#+ fig.cap=fig_cap 
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Year")
ggsave(here("figures/figure_3_figure_supplement_1.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ FigS2
f <- anticons ~ year + (year | de_tool)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_year__year_detool.rds"))
conditions <- make_conditions(data, vars = "de_tool")
row.names(conditions) <- conditions$de_tool
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions, 
                              re_formula = NULL),
          plot = FALSE)

fig_cap <- glue("Figure 3--figure supplement 2. A 2-level binomial logistic model *{deparse(f)}* 
                reveals that all differential expression analysis tools are associated with 
                temporally increasing anti-conservative p value histograms, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. 
                Lines denote best fit of linear model. Shaded area denotes 95% credible region.")

#+fig.cap=fig_cap
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative\np value histograms",
       x = "Year") +
  facet_wrap(~ de_tool) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(here("figures/figure_3_figure_supplement_2.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

#' 
#+ FigS3
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata)
f <- anticons ~ year + (year | model)
mod <- brm(formula = f, 
           data = data, 
           family = family,
           iter = ifelse(is_ci(), 400, 4000),
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_year__year_model.rds"))
conditions <- make_conditions(data, vars = "model")
row.names(conditions) <- str_wrap(conditions$model, width = 20)
p <- plot(conditional_effects(mod, 
                              effects = "year", 
                              conditions = conditions, 
                              re_formula = NULL),
          plot = FALSE)

fig_cap <- glue("Figure 3--figure supplement 3. A 2-level binomial logistic model *{deparse(f)}* reveals that all sequencing instrument models 
are associated with temporally increasing anti-conservative p value histograms, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Only GEO submissions utilizing 
single sequencing platform were used for model fitting. Lines denote best fit of linear model. Shaded area denotes 95% credible region.")

#+ fig.height=8, fig.cap=fig_cap
p$year + 
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(y = "Proportion of anti-conservative p value histograms",
       x = "Year") +
  facet_wrap(~ model, labeller = label_wrap_gen(width = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(here("figures/figure_3_figure_supplement_3.tiff"), height = 23, width = 18, dpi = 300, units = "cm")

#'
#+ FigS4
p <- data %>% 
  count(year, de_tool) %>% 
  add_count(year, name = "total", wt = n) %>% 
  mutate(Proportion = n / total) %>% 
  ggplot() +
  geom_line(aes(year, Proportion, color = de_tool), size = 1) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(2010, 2019, by = 3)) +
  labs(x = "Year") +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.8))

fig_cap <- glue("Figure 3--figure supplement 4. No single differential expression analysis tool dominates the field. 
                Y-axis shows the proportion of analysis platforms, 
                x-axis shows publication year of GEO submission, N = {prettyNum(sum(p$data$n), big.mark=',')}.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_3_figure_supplement_4.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#+ FigS5a
y_title <- "Prop. anti-cons."
f <- anticons ~ de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pa <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#'
#+ FigS5b
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_detool_all.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pb <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ FigS5c
f <- anticons ~ year + de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_year_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pc <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ FigS5d
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  ))
f <- anticons ~ organism + de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/anticons_organism_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pd <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ FigS5e
f <- anticons ~ de_tool + (1 | model)
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_detool__1_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pe <- p$de_tool + 
  labs(y = y_title) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

#' 
#+ FigS5f
f <- anticons ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons_detool__detool_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pf <- p$de_tool + 
  labs(y = y_title) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pfN <- prettyNum(summary(mod)$nobs, big.mark=',')
pfF <- deparse(f)

fig_cap <- glue('Figure 3--figure supplement 5. Binomial logistic models for proportion of anti-conservative p value histograms. 
                A, simple model *{paF}*, N = {paN}. 
                B, simple model *{pbF}* fitted on complete data, N = {pbN}. 
                C, model conditioned on year of GEO submission: *{pcF}*, N = {pcN}. 
                D, model conditioned on studied organism (human/mouse/other): *{pdF}*, N = {pdN}. 
                E, varying intercept model *{peF}* where "model" stands for sequencing instrument model, N = {peN}. 
                F, varying intercept and slope model *{pfF}*, N = {pfN}. Points denote best fit of linear model. Error bars, 95% credible interval.')

for (p in list(pa, pb, pc, pd, pe, pf)) {
  p$layers[[1]]$aes_params$size <- 1
}

#+ FigS5, fig.cap=fig_cap
(pa + pb + pc) / (pd + pe + pf) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
ggsave(here("figures/figure_3_figure_supplement_5.tiff"), height = 14, width = 18, dpi = 300, units = "cm")

#' 
#+ FigS6a
f <- pi0 ~ de_tool
family <- student()
data <- pvalues
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_detool_full_data.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pa <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
paN <- prettyNum(summary(mod)$nobs, big.mark=',')
paF <- deparse(f)

#' 
#+ FigS6b
f <- pi0 ~ year + de_tool
data <- pvalues_sample
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           file = here("results/models/pi0_year_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pb <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pbN <- prettyNum(summary(mod)$nobs, big.mark=',')
pbF <- deparse(f)

#' 
#+ FigS6c
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata) %>% 
  mutate(organism = case_when(
    tax_id == 9606 ~ "human",
    tax_id == 10090 ~ "mouse",
    TRUE ~ "other"
  ))
f <- pi0 ~ organism + de_tool
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           file = here("results/models/pi0_organism_detool.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pc <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pcN <- prettyNum(summary(mod)$nobs, big.mark=',')
pcF <- deparse(f)

#' 
#+ FigS6d
f <- pi0 ~ de_tool + (1 | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0_detool__1_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pd <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdN <- prettyNum(summary(mod)$nobs, big.mark=',')
pdF <- deparse(f)

#' 
#+ FigS6e
f <- pi0 ~ de_tool + (de_tool | model)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0_detool__detool_model.rds"))
p <- plot(conditional_effects(mod, 
                              effects = "de_tool"),
          plot = FALSE)
pe <- p$de_tool + 
  labs(y = expression(pi[0])) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
peN <- prettyNum(summary(mod)$nobs, big.mark=',')
peF <- deparse(f)

for (p in list(pa, pb, pc, pd, pe)) {
  p$layers[[1]]$aes_params$size <- 1
}
p <- (pa + pb + pc) / (pd + pe + plot_spacer()) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))
fig_cap <- glue("Figure 4--figure supplement 1. Robust (student's t likelihood) modeling of $\\pi_0$. 
                A, simple model *{paF}* fitted on complete data, N = {paN}. 
                B, model conditioned on year of GEO submission: *{pbF}*, N = {pbN}. 
                C, model conditioned on studied organism (human/mouse/other): *{pcF}*, N = {pcN}. 
                D, varying intercept model *{pdF}* where 'model' stands for sequencing instrument model, N = {pdN}. 
                E, varying intercept/slope model *{peF}*, N = {peN}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_4_figure_supplement_1.tiff"), height = 14, width = 18, dpi = 300, units = "cm")

#' 
#+ FigS7
f <- pi0 ~ (1 | model)
data <- data %>% 
  filter(!is.na(pi0), !is.na(model))
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__1_model.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_model[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_model) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = expression(pi[0])) +
  theme(axis.title.y = element_blank())
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 4--figure supplement 2. Modeling dependency of $\\pi_0$ on sequencing instrument model: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_4_figure_supplement_2.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

#' 
#+ FigS8
f <- pi0 ~ (1 | library_strategy)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__1_librarystrategy.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_strategy[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_strategy) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = expression(pi[0]), y = "Library strategy")
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 4--figure supplement 3. Modeling dependency of $\\pi_0$ on library strategy: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_4_figure_supplement_3.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ FigS9
f <- pi0 ~ (1 | library_selection)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores, 
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__1_libraryselection.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_selection[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_selection) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = expression(pi[0]), y = "Library selection")
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 4--figure supplement 4. Modeling dependency of $\\pi_0$ on library selection: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_4_figure_supplement_4.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+  FigS10
f <- pi0 ~ (1 | library_layout)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           iter = ifelse(is_ci(), 400, 4000),
           refresh = refresh,
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/pi0__1_librarylayout.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_layout[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_layout) %>%
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = expression(pi[0]), y = "Library layout")
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 4--figure supplement 5. Modeling dependency of $\\pi_0$ on library layout: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_4_figure_supplement_5.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ FigS11
data <- pvalues_sample %>% 
  inner_join(sequencing_metadata)
f <- anticons ~ (1 | model)
family <- bernoulli()
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__1_model.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_model[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_model) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative histograms") +
  theme(axis.title.y = element_blank())
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 3--figure supplement 6. Modeling dependency of proportion of anti-conservative histograms on sequencing platform: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_3_figure_supplement_6.tiff"), height = 12, width = 18, dpi = 300, units = "cm")

#' 
#+ FigS12
f <- anticons ~ (1 | library_strategy)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__1_librarystrategy.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_strategy[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_strategy) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative histograms", y = "Library strategy") +
  scale_x_continuous(limits = c(0, 1))
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 3--figure supplement 7. Modeling dependency of proportion of anti-conservative histograms on library strategy: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_3_figure_supplement_7.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ FigS13
f <- anticons ~ (1 | library_selection)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__1_libraryselection.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_selection[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_selection) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative histograms", y = "Library_selection") +
  scale_x_continuous(limits = c(0, 1))
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 3--figure supplement 8. Modeling dependency of proportion of anti-conservative histograms on library selection: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_3_figure_supplement_8.tiff"), height = 7, width = 10, dpi = 300, units = "cm")

#' 
#+ FigS14
f <- anticons ~ (1 | library_layout)
mod <- brm(formula = f, 
           data = data, 
           family = family, 
           chains = chains, 
           cores = cores,
           refresh = refresh,
           iter = ifelse(is_ci(), 400, 2000),
           control = list(adapt_delta = 0.99, max_treedepth = 12),
           file = here("results/models/anticons__1_librarylayout.rds"))
p <- mod %>%
  spread_draws(b_Intercept, r_library_layout[condition,]) %>%
  median_hdci(condition_mean = b_Intercept + r_library_layout) %>%
  mutate_at(vars(condition_mean, .lower,  .upper), inv_logit_scaled) %>% 
  ggplot(aes(y = condition, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  labs(x = "Proportion of anti-conservative histograms", y = "Library layout")
p$layers[[1]]$aes_params$size <- 1

fig_cap <- glue("Figure 3--figure supplement 9. Modeling dependency of proportion of anti-conservative histograms on library layout: *{deparse(f)}*, N = {prettyNum(summary(mod)$nobs, big.mark=',')}. Points denote best fit of linear model. Error bars denote 95% credible interval.")

#+ fig.cap=fig_cap
p
ggsave(here("figures/figure_3_figure_supplement_9.tiff"), height = 7, width = 10, dpi = 300, units = "cm")
