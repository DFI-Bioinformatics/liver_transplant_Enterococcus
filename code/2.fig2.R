library(tidyverse)

# load meta and imputed data -----------------------------------------------

lt.sams <- read_csv("data/lt.meta_rawMetabolomics.csv")

qual.imp <- read_csv("data/peri_post.qual.imput.2026-05-12.ori.csv")

ba.cat <- read_csv("data/Bile Acid Identities Qual.csv") |> 
  filter(hmmf_panel == "BileAcid") |> 
  select(hmmf_panel, compound, Category) |> 
  mutate(Category = gsub("Secondsry", "Secondary", Category),
         Category = gsub(" BA","", Category)) |> 
  mutate(Category = factor(Category,
                           levels = c("Unconjugated Primary",
                                      "Conjugated Primary",
                                      "Unconjugated Secondary",
                                      "Conjugated Secondary")))|> 
  add_count(Category, name = "nCat")

ba.imp <- qual.imp |> 
  gather("compound", "value", -sample_id) |> 
  right_join(ba.cat |> 
               select(-c(hmmf_panel, nCat)))

# transformation test: zscore, minmax and rank -------------------------------

ba.cat.norm <- ba.imp %>%
  
  # STEP 1: Apply transformations to each compound individually
  # Grouping by compound ensures the scale normalizes for different ionization efficiencies
  group_by(compound) %>%
  mutate(
    # 1. Original Method: Min-Max Scaling (0-100 scale)
    scaled_minmax = (value - min(value, na.rm = TRUE)) / 
      (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) * 100,
    
    # 3. Reviewer Method B: Standardized Z-Score
    # (as.numeric prevents matrix formatting issues in downstream dplyr operations)
    scaled_zscore = as.numeric(scale(value)),
    
    # 4. Reviewer Method C: Percentile Rank (0-1 scale)
    # (Completely robust to magnitude outliers)
    scaled_pct_rank = percent_rank(value)
  ) %>%
  ungroup() %>%
  
  # STEP 2: Additive pooling of the transformed scores
  # Group by patient (sample_id) and bile acid class (category)
  group_by(sample_id, Category) %>%
  
  # Sum the individually transformed values to create Composite Scores
  summarize(
    composite_minmax = sum(scaled_minmax, na.rm = TRUE),
    composite_zscore = sum(scaled_zscore, na.rm = TRUE),
    composite_pct_rank = sum(scaled_pct_rank, na.rm = TRUE),
    
    .groups = "drop"
  )

ba.total.norm <- ba.imp |> 
  group_by(sample_id) %>%
  summarize(total_abundance = sum(value, na.rm = TRUE), .groups = "drop") %>%
  
  # STEP 2: Apply mathematical transformations
  # CRITICAL: Group by category so that z-scores and ranks are calculated 
  # relative to other samples WITHIN the same bile acid class, rather than across classes.
  mutate(
    # 1. Original Method: Min-Max Scaling (0-100 scale)
    scaled_minmax = (total_abundance - min(total_abundance, na.rm = TRUE)) / 
      (max(total_abundance, na.rm = TRUE) - min(total_abundance, na.rm = TRUE)) * 100,
    
    # 3. Reviewer Method B: Standardized Z-Score
    # (as.numeric prevents matrix formatting issues in downstream dplyr operations)
    scaled_zscore = as.numeric(scale(total_abundance)),
    
    # 4. Reviewer Method C: Percentile Rank (0-1 scale)
    # (Completely robust to magnitude outliers)
    scaled_pct_rank = percent_rank(total_abundance)
  ) %>%
  ungroup() |> 
  mutate(Category = "Total BA")

ba.trans.long.df <- ba.total.norm |> 
  gather("type", "trans_values", -c(sample_id, Category)) |> 
  bind_rows(ba.cat.norm |> 
              gather("type","trans_values",-c(sample_id, Category))) |> 
  mutate(new.cat = case_when(
    grepl("^scaled", type) ~ paste0(Category, " (normalized)"),
    grepl("^total", type) ~ paste0(Category, " (raw)"),
    TRUE ~ paste0(Category, "\n(normalized composite)")
  ),
         new.type = str_split(type, "_") |> sapply(function(x) x[length(x)]),
          new.cat = factor(new.cat, 
                           levels = c("Total BA (raw)",
                                      "Total BA (normalized)",
                                      "Unconjugated Primary\n(normalized composite)",
                                      "Conjugated Primary\n(normalized composite)",
                                      "Unconjugated Secondary\n(normalized composite)",
                                      "Conjugated Secondary\n(normalized composite)"))) 

ba.trans.long.df |> 
  count(new.type)

ba.trans.long.df |> 
  count(new.type, new.cat)

ba.trans.long.df |> 
  count(new.cat)

# load mixed effect model libraries ------------------------------------------

library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)

# Fig 2: pre vs post (combined) ----------------------------------------------

## panel A -----------------------------------------------------------------

pp.ba.df <- ba.trans.long.df |> 
  filter(new.type %in% c("abundance", "minmax")) |> 
  left_join(lt.sams |> 
              select(sample_id = metabolomics_id, 
                     study_id = studyid, timepoint, 
                     timepoint2, death))

pp.minmax.mix <- pp.ba.df |> 
  mutate(
    # Set "Early" as the reference level so the model calculates the change to "Late"
    timepoint2 = factor(timepoint2, levels = c("pre-transplant", "post-transplant")),
    study_id = factor(study_id)
  ) |> 
  group_by(new.type, new.cat) %>%
  # Collapse the rest of the data into a "list-column" for each taxon
  nest() %>%
  mutate(
    # Because lmerTest is loaded, this lmer() now calculates degrees of freedom!
    model = map(data, ~ lmer(trans_values ~ timepoint2 + (1 | study_id), data = .x)),
    tidied = map(model, ~ tidy(.x, conf.int = TRUE))
  ) %>%
  unnest(tidied) %>%
  filter(effect == "fixed") %>%
  filter(term != "(Intercept)") %>%
  # The p.value column will now successfully pull through
  select(new.type, new.cat, term, estimate, std.error, 
         statistic, p.value, conf.low, conf.high) |> 
  ungroup()

pp.minmax.mix.res <- pp.minmax.mix |> 
  group_by(new.type, new.cat) |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  
  # Optional: Reorder the columns so p.adj sits right next to the raw p.value
  relocate(p.adj, .after = p.value) %>%
  
  # Optional: Sort by the adjusted p-value to bring the most significant taxa to the top
  arrange(p.adj) |> 
  left_join(pp.ba.df |> 
              group_by(new.type, new.cat) |> 
              summarise(y.position = max(trans_values)*1.12)) |> 
  mutate(group1 = "pre-transplant", group2 = "post-transplant") |> 
  add_significance(p.col = "p.adj") |> 
  ungroup()

pp.gg <- pp.ba.df |> 
  mutate(timepoint2 = factor(timepoint2, 
                             levels = c("pre-transplant", "post-transplant"))) |> 
  ggplot(aes(timepoint2, trans_values)) +
  geom_violin(aes(color = timepoint2)) +
  geom_boxplot(aes(color = timepoint2), width = 0.65, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(aes(fill = timepoint2), shape = 21, size = 2.5, alpha = 0.65) +
  facet_wrap(. ~ new.cat, scales = "free", nrow = 4) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 9.5),
        axis.text.y=element_text(size = 9.5),
        legend.position="none",
        strip.text.x = element_text(angle=0, size = 11, face = "bold"),
        strip.background = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = c("#00C1CF", "#CF7500")) +
  scale_fill_manual(values = c("#00C1CF", "#CF7500")) +
  ggpubr::stat_pvalue_manual(pp.minmax.mix.res, 
                             label = "q-value: {scales::pvalue(p.adj)}{p.adj.signif}",
                             size = 3.5,
                             hide.ns = TRUE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs( x = "Timepoint", y = "Min-Max transformed values")

## panel B -----------------------------------------------------------------

ba.death.df <- pp.ba.df |>
  filter(timepoint == "Late") 

ba.death.res <- ba.death.df |> 
  group_by(new.type, new.cat) |> 
  wilcox_test(trans_values ~ death)|>
  mutate(p.adj = p.adjust(p, method = "BH")) |> 
  add_significance(p.col = "p.adj") |> 
  left_join(ba.death.df |> 
              group_by(new.type, new.cat) |> 
              summarise(y.position = max(trans_values)) ) |> 
  mutate(group1 = "Survival", group2 = "Death")

death.gg <- ba.death.df  |> 
  mutate(death = if_else(death == "Yes", "Death", "Survival")) |> 
  ggplot(aes(death, trans_values, color = death)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.45) +
  facet_wrap(. ~ new.cat, scales = "free", nrow = 4) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = d.pal,) +
  ggpubr::stat_pvalue_manual(data = ba.death.res,
                             label = "q-value: {scales::pvalue(p.adj)}{p.adj.signif}",
                             size = 3.5,
                             hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(axis.text.x=element_text(size = 9.5),
        axis.text.y=element_text(size = 9.5),
        legend.position="none",
        strip.text.x = element_text(angle=0, size = 11, face = "bold"),
        strip.background = element_blank(),
        axis.title = element_blank()) +
  scale_x_discrete(labels = c("Non-survivor", "Survivor"))

## combine -----------------------------------------------------------------

library(patchwork)

pp.gg + death.gg +
  plot_annotation(tag_levels = "A")

ggsave("figures/Fig2_minmax_qualPoolBA.ori.pdf", 
        height = 8.5, width = 11)
