library(tidyverse)
library(forcats)
library(cutpointr)

# ROC curves of tyrosine, tyramine, and tyramine/tyrosine ratio to predict cluster E vs D

# load data ---------------------------------------------------------------

## load meta ----------
meta.lt <- read_csv("data/lt.meta_rawMetabolomics.csv")

## load imputed qual data ----

qual.imp <- read_csv("data/peri_post.qual.imput.2026-05-12.ori.csv")

## umap cluster ------------------------------------------------------------

u.meta <- read_csv("data/paired.umap_cluster.2025-03-12.csv") |> 
  mutate(cl = if_else(kcluster == 1, "Cluster D", "Cluster E"),
         kcluster = factor(kcluster)) |> 
  select(-c(paired, timepoint)) |> 
  distinct()

# assemble data -----------------------------------------------------------

nor_tyo <- qual.imp$tyrosine/max(qual.imp$tyrosine)*100

nor_tya <- qual.imp$tyramine/max(qual.imp$tyramine)*100

ready.df <- qual.imp |> 
  transmute(metabolomics_id = sample_id, 
            tyrosine, tyramine, `tyramine/tyrosine ratio` = nor_tya/nor_tyo) |> 
  gather("qual_metab", "value", -metabolomics_id) |> 
  left_join(meta.lt |> 
              select(metabolomics_id, shotgun_seq_id)) |> 
  left_join(u.meta |> 
              select(shotgun_seq_id, cl)) |> 
  mutate(cl_bi = if_else(cl == "Cluster D", 0, 1))

safe_cutpointr <- possibly(.f = cutpointr, otherwise = "Error")

cutpoints <-
  ready.df %>%
  group_by(qual_metab) %>%
  group_map(
    ~ safe_cutpointr(
      .,
      value,
      cl_bi,
      qual_metab,
      method = maximize_metric,
      metric = youden,
      pos_class = 1,
      neg_class = 0,
      boot_runs = 100,
      boot_stratify = TRUE,
      use_midpoints = TRUE,
      na.rm = TRUE
    ),
    .keep = TRUE
  )

ac_ci.df <- lapply(cutpoints, boot_ci, acc) |> 
  bind_rows() |> 
  spread(quantile, values) |> 
  mutate(acc_ci = paste0(" (",sprintf("%.2f",`0.025`),", ",sprintf("%.2f",`0.975`),")"))|> 
  select(subgroup, ends_with("ci"))

sens_ci.df <- lapply(cutpoints, boot_ci, sensitivity) |> 
  bind_rows() |> 
  spread(quantile, values) |> 
  mutate(sens_ci = paste0(" (",sprintf("%.2f",`0.025`),", ",sprintf("%.2f",`0.975`),")"))|> 
  select(subgroup, ends_with("ci"))

spec_ci.df <- lapply(cutpoints, boot_ci, specificity) |> 
  bind_rows() |> 
  spread(quantile, values) |> 
  mutate(spec_ci = paste0(" (",sprintf("%.2f",`0.025`),", ",sprintf("%.2f",`0.975`),")")) |> 
  select(subgroup, ends_with("ci"))

auc_ci.df <- lapply(cutpoints, boot_ci, AUC)|> 
  bind_rows() |> 
  spread(quantile, values) |> 
  mutate(auc_ci = paste0(" (",sprintf("%.2f",`0.025`),", ",sprintf("%.2f",`0.975`),")")) |> 
  select(subgroup, ends_with("ci"))

cutpoints_unnest <- cutpoints %>%
  map_df(as_tibble)

cp.df <- cutpoints_unnest %>%
  select(
    subgroup,
    direction,
    optimal_cutpoint,
    acc, sensitivity, specificity, AUC
  ) |> 
  left_join(ac_ci.df) |> 
  left_join(sens_ci.df) |> 
  left_join(spec_ci.df) |> 
  left_join(auc_ci.df) |> 
  mutate(label = paste0("Cutpoint: ", sprintf("%.2f", optimal_cutpoint),"\n",
                        "Accuracy: ", sprintf("%.2f", acc), acc_ci, "\n",
                        "Sensitivity: ", sprintf("%.2f", sensitivity),sens_ci,"\n",
                        "Specificity: ", sprintf("%.2f", specificity),spec_ci,"\n",
                        "AUC: ", sprintf("%.2f", AUC), auc_ci))

cp.gg <- cutpoints_unnest |> 
  unnest(roc_curve) |> 
  ggplot(aes(fpr,tpr)) +
  geom_line(aes(group = subgroup)) +
  facet_wrap(.~ subgroup, scales = "free") +
  geom_abline(intercept = 0, color = "red", linetype = 2) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_blank()) +
  geom_text(data = cp.df,
           aes(x = 1, y = 0.3, label = label),
           hjust = 1, vjust = 1) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(panel.grid = element_blank())

cp.gg

## box/violin ----------------------------------------------------------------

source("code/0.pals.R")

dis.gg <- ready.df |> 
  ggplot(aes(cl, value)) +
  geom_violin(aes( color = cl)) +
  geom_boxplot(aes( color = cl), outlier.shape = NA, width = 0.65, alpha = 0.65) +
  geom_jitter(aes( fill = cl), alpha = 0.65, shape = 21) + 
  theme_bw() +
  labs(x = "Cluster", y = "log10 transformed") +
  ggpubr::stat_compare_means() +
  scale_color_manual(values = c.pal) +
  scale_fill_manual(values = c.pal) +
  scale_y_log10() +
  facet_wrap(.~ qual_metab, scales = "free") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  geom_hline(data = cutpoints_unnest |> 
               mutate(qual_metab = subgroup), 
             aes(yintercept = optimal_cutpoint),
             color = "red", linetype = 2)
  
dis.gg

library(patchwork)

cp.gg / dis.gg +
  plot_annotation(tag_levels = "A")

ggsave("figures/Fig5_ttroc.box.pdf", height = 8.5, width = 12.5)
