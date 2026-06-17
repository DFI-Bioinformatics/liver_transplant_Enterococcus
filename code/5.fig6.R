library(tidyverse)
library(scales)
library(cutpointr)
library(epiR)
library(patchwork)
library(caret)

# check relative abundance of E. faecium in clusters

source("code/0.pals.R")

u.meta <- read_csv("data/paired.umap_cluster.2025-03-12.csv") |> 
  mutate(cl = if_else(kcluster == 1, "Cluster D", "Cluster E"),
         kcluster = factor(kcluster)) |> 
  select(-c(paired, timepoint)) |> 
  distinct()

mpa_plt <- read_csv("data/mpa.june23.csv") |> 
  filter(!is.na(Kingdom)) |> 
  add_count(seq_id, wt = clean_relative_abundance, name = "totalAbd") |> 
  mutate(pctseqs = clean_relative_abundance/totalAbd)

mpa_plt |> 
  count(Species, wt = pctseqs, sort = T)

ef.abd <- mpa_plt |> 
  filter(Species == "Enterococcus_faecium") |> 
  count(seq_id, Species, wt = pctseqs, name = "pctseqs") |> 
  right_join(u.meta, by = c("seq_id" = "shotgun_seq_id")) |> 
  replace_na(list(pctseqs = 0)) 

cp.res <- cutpointr(ef.abd$pctseqs, ef.abd$cl, metric = youden) |> 
  mutate(label = paste0("Cutpoint: ", sprintf("%.3f", optimal_cutpoint),"\n",
                 "Accuracy: ", sprintf("%.3f", acc),"\n",
                 "Sensitivity: ", sprintf("%.3f", sensitivity),"\n",
                 "Specificity: ", sprintf("%.3f", specificity),"\n",
                 "AUC: ", sprintf("%.3f", AUC)))

ef.gg <- ef.abd |> 
  ggplot(aes(cl, pctseqs*100)) +
  geom_boxplot(aes( color = cl), outlier.shape = NA) +
  geom_jitter(aes(fill = cl), size = 2.5, 
              shape = 21, height = 0, width = 0.25, alpha = 0.65) +
  theme_bw() +
  labs(x = "Cluster", y = "Abund.",
       title = 'Optimal cutpoint for E. faecium abundance for Cluster E and D') +
  ggpubr::stat_compare_means() +
  scale_color_manual(values = c.pal) +
  scale_fill_manual(values = c.pal) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_hline(aes(yintercept = cp.res$optimal_cutpoint*100),
             color = "red", linetype = 2) +
  scale_y_continuous(breaks = c(0, 0.2, signif(cp.res$optimal_cutpoint, 3),
                                0.4, 0.6, 0.8, 1)*100) +
  annotate("text", y = 0.25, x = "Cluster E", 
           label = "optimal cutpoint", color = "red")

ef.gg


# validation cohort -----------------------------------------

val.tt_ratio <- read_csv("data/Fig6.validation_cohort.csv")

sublab <- val.tt_ratio |> 
  count(E_cl) |> 
  mutate(label = paste0(E_cl, " = ", n)) |> 
  pull(label) |> 
  paste0(collapse = "; ")

cm <- confusionMatrix(table(val.tt_ratio$tt_cl, val.tt_ratio$cl))

tt.epi <- epiR::epi.tests(table(val.tt_ratio$tt_cl, val.tt_ratio$cl), conf.level = 0.95)

tt.metrics <- tt.epi$detail |> 
  as_tibble() |> 
  filter(statistic %in% c("sp", "se")) |> 
  mutate(statistic = if_else(statistic == "sp", "Specificity", "Sensitivity")) |> 
  add_row(statistic = "Accuracy",
          est = cm$overall[1],
          lower = cm$overall[3],
          upper = cm$overall[4]) |> 
  mutate(lab = paste0(signif(est,3), " (",signif(lower,3),", ", signif(upper,3),")"))

tt.metrics <- tt.metrics |> 
  slice(c(3,1,2))

tt_ratio.metrics.lab <- paste0("Cutpoint=0.25\n",paste0(tt.metrics$statistic, "=", tt.metrics$lab,
       collapse = "\n"))

tt_ratio.gg <- val.tt_ratio |> 
  ggplot(aes(E_cl, `tyramine/tyrosine ratio`)) +
  geom_violin(aes( color = E_cl)) +
  geom_boxplot(aes( color = E_cl), outlier.shape = NA, width = 0.65, alpha = 0.65) +
  geom_jitter(aes(fill = E_cl), size = 2.5,
              shape = 21, height = 0, width = 0.25, alpha = 0.35) +
  theme_bw() +
  labs(x = "", y = "log10 transformed",
       title = 'tyramine/tyrosine ratio in validation cohort',
       caption = sublab) +
  ggpubr::stat_compare_means() +
  scale_color_manual(values = c.pal.new) +
  scale_fill_manual(values = c.pal.new) +
  scale_y_log10() +
  geom_hline(aes(yintercept = 0.25), color = "red", linetype = 2) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  annotate("text",
            x = "\u226526.6% E. faecium", y = 0.000001,
           label = tt_ratio.metrics.lab,
            hjust = 0.5, vjust = 0)

ef.gg +  tt_ratio.gg +
  plot_annotation(tag_levels = "A")

ggsave("figures/Fig6_E.faecium_validation.jitter.pdf", height = 5.5, width = 11)
