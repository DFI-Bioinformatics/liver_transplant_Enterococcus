library(tidyverse)
library(vegan)
library(remotes)

source("https://raw.githubusercontent.com/yingeddi2008/DFIutility/master/getRdpPal.R")

# load meta ---------------------------------------------------------------

meta_all <- read_csv("data/3grp_A.sampeList.ori.csv")

source("code/0.pals.R")

# get metaphlan ----------------------------------------------------

mpa <- read_csv("data/mpa.june23.csv")

mpa.plt <- mpa |> 
  add_count(seq_id, 
            wt = clean_relative_abundance, 
            name = "totalAbd") %>%
  mutate(pctseqs = clean_relative_abundance/totalAbd,
         genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

# inverse simpson ---------------------------------------------------------

mpa.mat <- mpa %>% 
  select(seq_id, taxid, clean_relative_abundance) %>% 
  pivot_wider(names_from = taxid, 
              values_from = clean_relative_abundance, 
              values_fill = 0,
              values_fn = sum) %>% 
  column_to_rownames(var = "seq_id")

mpa.mat.trans <- round(mpa.mat/rowSums(mpa.mat)* 1000000)

invt <- diversity(mpa.mat.trans, index = "inv")

invt.df <- tibble(shotgun_seq_id = names(invt),
                  invSimp =  invt)

taxpal <- getRdpPal(mpa.plt)

mpa.gg <-  mpa.plt %>%
  group_by(seq_id, Kingdom,Phylum,Class,Order,Family, Genus, genLab) %>%
  summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() %>% # view
  left_join(meta_all |> 
              select(seq_id = shotgun_seq_id,
                     sample_id, invSimp, header)) |> 
  mutate(header = factor(header, levels =  c("pre-transplant (47)", 
                                             "early post-transplant (57)",
                                             "late post-transplant (86)",
                                             "Healthy Donors (36)"))) |> 
  ggplot(aes(x=reorder(sample_id, invSimp),y=pctseqs)) +
  geom_col(aes(fill=Genus),position ="fill") +
  scale_fill_manual(values=taxpal) +
  labs( y = "Metaphlan4 (v.Jun23) Bacteria Abundance",
        caption = paste0("Number of samples:", nrow(mpa.plt %>% count(seq_id)))) +
  theme_bw() +
  theme(# axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position="none",
    strip.text.x = element_text(angle=0, size = 9.5, face = "bold"),
    strip.background = element_blank()) +
  facet_grid(. ~ header, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0.001,0.001))

mpa.gg

# sequencing depth -------------------------------------------------------
## by Cluster D vs. E ----------------------------------------------------

u.meta <- read_csv("data/paired.umap_cluster.2025-03-12.csv") |> 
  mutate(cl = if_else(kcluster == 1, "Cluster D", "Cluster E"),
         kcluster = factor(kcluster)) |> 
  select(-c(paired, timepoint)) |> 
  distinct()

meta_all |> 
  inner_join(u.meta) |> 
  ggplot(aes(cl, count, fill = cl)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 21) +
  scale_fill_manual(values = c.pal) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  ggpubr::stat_compare_means() +
  labs(y = "Sequencing Depth")

ggsave("figures/SuppFig6_seqDepth.byCluster.pdf", height = 5.8, width = 5.8)

# bacteria abundance ------------------------------------------------------

bact.abd <- mpa.plt |> 
  filter(genLab %in% c("Enterococcus")) |> 
  count(seq_id, genLab, wt = pctseqs) |> 
  dplyr::rename(Taxon = genLab) |> 
  bind_rows(mpa.plt |> 
              filter(Family %in% c("Lachnospiraceae",
                                   "Bifidobacteriaceae",
                                   "Enterobacteriaceae",
                                   "Oscillospiraceae")) |> 
              count(seq_id, Family, wt = pctseqs) |> 
              dplyr::rename(Taxon = Family)) |>
  bind_rows(mpa.plt |> 
              filter(Phylum %in% c("Bacteroidota")) |> 
              count(seq_id, Phylum, wt = pctseqs) |> 
              dplyr::rename(Taxon = Phylum)) |> 
  spread(Taxon, n, fill = 0)

bact.df <- meta_all |>
  filter(timepoint3 != "Healthy Donors") |> 
  select(shotgun_seq_id, sample_id, timepoint3, study_id) |> 
  filter(shotgun_seq_id %in% mpa.plt$seq_id) |> 
  left_join(invt.df) |> 
  select(seq_id = shotgun_seq_id,
         `Inverse Simpson` = invSimp,
         timepoint3, study_id) |> 
  left_join(bact.abd) %>%
  replace(is.na(.),0) |> 
  gather("taxon", "value", -c(seq_id, timepoint3, study_id)) |> 
  mutate(taxon = factor(taxon,
                        levels = c("Inverse Simpson",
                                   "Enterococcus",
                                   "Enterobacteriaceae",
                                   "Lachnospiraceae",
                                   "Bifidobacteriaceae",
                                   "Bacteroidota",
                                   "Oscillospiraceae"))) 

# mixed effect model ------------------------------------------------------

library(lme4)
library(lmerTest)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)

taxa_results <- bact.df |> 
  mutate(
    # Set "Early" as the reference level so the model calculates the change to "Late"
    timepoint = factor(timepoint3, levels = c("early post-transplant","late post-transplant", "pre-transplant")),
    study_id = factor(study_id)
  ) |> 
  group_by(taxon) %>%
  # Collapse the rest of the data into a "list-column" for each taxon
  nest() %>%
  mutate(
    # Because lmerTest is loaded, this lmer() now calculates degrees of freedom!
    model = map(data, ~ lmer(value ~ timepoint3 + (1 | study_id), data = .x)),
    tidied = map(model, ~ tidy(.x, conf.int = TRUE))
  ) %>%
  unnest(tidied) %>%
  filter(effect == "fixed") %>%
  filter(term != "(Intercept)") %>%
  # The p.value column will now successfully pull through
  select(taxon, term, estimate, std.error, statistic, p.value, conf.low, conf.high)

final_stats <- taxa_results %>%
  ungroup() |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  
  # Optional: Reorder the columns so p.adj sits right next to the raw p.value
  relocate(p.adj, .after = p.value) %>%
  
  # Optional: Sort by the adjusted p-value to bring the most significant taxa to the top
  arrange(p.adj)

bact.pvals <- final_stats |> 
  add_significance(p.col = "p.adj") |> 
  mutate(group1 = "early post-transplant", 
         group2 = if_else(grepl("late post-transplant",term), "late post-transplant",
                          "pre-transplant")) |> 
  left_join(bact.df |> 
              group_by(taxon, timepoint3) |> 
              summarise(y.position = max(value)*1.12) |> 
              rename(group2 = timepoint3)) 

bact.gg <- bact.df |> 
  mutate(timepoint3 = factor(timepoint3, 
                             levels = c("pre-transplant",
                                        "early post-transplant",
                                        "late post-transplant"))) |> 
  ggplot(aes(timepoint3, value)) +
  geom_violin(aes(color = timepoint3)) +
  geom_boxplot(aes(color = timepoint3), width = 0.65, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(aes(fill = timepoint3), shape = 21, size = 2.5, alpha = 0.65) +
  facet_wrap(.~ taxon, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 9.5, angle = 45, hjust = 1),
        axis.text.y=element_text(size = 9.5),
        legend.position="none",
        strip.text.x = element_text(angle=0, size = 11, face = "bold"),
        strip.background = element_blank()) +
  scale_color_manual(values = t2.pal) +
  scale_fill_manual(values = t2.pal) +
  ggpubr::stat_pvalue_manual(bact.pvals, 
                             label = "q-value: {scales::pvalue(p.adj)}{p.adj.signif}",
                             size = 3.5,
                             hide.ns = TRUE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs( x = "", y = "")

bact.gg

# inverse simpson with paired ---------------------------------------------

bact.df |> 
  mutate(timepoint3 = factor(timepoint3, 
                             levels = c("pre-transplant",
                                        "early post-transplant",
                                        "late post-transplant"))) |> 
  filter(taxon == "Inverse Simpson") |> 
  ggplot(aes(timepoint3, value)) +
  geom_violin(aes(color = timepoint3)) +
  geom_boxplot(aes(color = timepoint3), 
               width = 0.55, alpha = 0.55, outlier.shape = NA) +
  geom_line(aes(group = study_id), linetype = 2, alpha = 0.6) +
  geom_jitter(aes(fill = timepoint3), shape = 21, size = 2.5, alpha = 0.65,
              width = 0.15) +
  facet_wrap(.~ taxon, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 9.5, angle = 45, hjust = 1),
        axis.text.y=element_text(size = 9.5),
        legend.position="none",
        strip.text.x = element_text(angle=0, size = 11, face = "bold"),
        strip.background = element_blank()) +
  scale_color_manual(values = t2.pal) +
  scale_fill_manual(values = t2.pal) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs( x = "", y = "")

ggsave("figures/SuppFig2_3grp.invSimp.ori.pdf", height = 6.35, width = 4.8)

# assemble Fig 1 A + B ----------------------------------------------------

library(patchwork)

(mpa.gg +
    labs(caption = "") +
    theme(strip.text.x = element_text(size = 11)))  / bact.gg +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 0.85))

ggsave("figures/Fig1_AB.pdf", height = 8.95, width = 13.5 )

# Fig 1 C: qual volcano plot -------------------------------------------------------

qual.imp <- read_csv("data/peri_post.qual.imput.2026-05-12.ori.csv")
lt.sams <- read_csv("data/lt.meta_rawMetabolomics.csv")

qual.imp.long <-  qual.imp |> 
  gather("compound","value", -sample_id) |> 
  as_tibble()

qual_v_df <- qual.imp.long |> 
  left_join(lt.sams |> 
              transmute(sample_id = sampleid, 
                        timepoint2, 
                        study_id = studyid))

qual_tp_df <- qual_v_df |> 
  mutate(
    timepoint2 = factor(timepoint2, levels = c("pre-transplant", "post-transplant")),
    study_id = factor(study_id)
  ) |> 
  group_by(compound) %>%
  # Collapse the rest of the data into a "list-column" for each taxon
  nest() %>%
  mutate(
    # Because lmerTest is loaded, this lmer() now calculates degrees of freedom!
    model = map(data, ~ lmer(value ~ timepoint2 + (1 | study_id), data = .x)),
    tidied = map(model, ~ tidy(.x, conf.int = TRUE))
  ) %>%
  unnest(tidied) %>%
  filter(effect == "fixed") %>%
  filter(term != "(Intercept)") %>%
  # The p.value column will now successfully pull through
  select(compound, term, estimate, std.error, statistic, p.value, conf.low, conf.high) |> 
  ungroup() |> 
  mutate(p.adj = p.adjust(p.value, method = "BH")) 

qual_tp_fc <- qual_v_df |> 
  group_by(compound, timepoint2) |> 
  summarise(ave = mean(value)) |> 
  mutate(ave = if_else(ave == 0, 0.0001, ave)) |> 
  spread(timepoint2, ave) |> 
  mutate(log2fc = log2(`post-transplant`/`pre-transplant`))

qual_tp_tot <- left_join(qual_tp_fc, qual_tp_df) %>%
  column_to_rownames(var = "compound")

#### p-adjusted -------
set.seed(123456)

xylims = ceiling(max(abs(qual_tp_tot$log2fc)))
plims <- ceiling(max(abs(-log10(qual_tp_tot$p.adj))))

pcut = 0.05
fccut = 1.5

library(EnhancedVolcano)

qual_tp_volcano <-
  EnhancedVolcano(qual_tp_tot,
                  lab = rownames(qual_tp_tot),
                  title = "",
                  y = "p.value",
                  x = "log2fc",
                  pCutoff = pcut,
                  FCcutoff = fccut,
                  pointSize = 3.4,
                  labSize = 3.5,
                  xlim = c(-xylims,xylims),
                  ylim = c(0, plims + 1),
                  col=c("#F0F0F0FF","#D2D2D2FF","#77AB43FF","#EA879CFF"),
                  colAlpha = 0.85,
                  legendPosition = "right",
                  legendLabels = c(bquote(p > .(pcut)*";" ~ Log[2] ~ FC < "\u00B1"*.(fccut)),
                                   bquote(p > .(pcut)*";" ~ Log[2] ~ FC >= "\u00B1"*.(fccut)),
                                   bquote(p <= .(pcut)*";" ~ Log[2] ~ FC < "\u00B1"*.(fccut)),
                                   bquote(p <= .(pcut)*";" ~ Log[2] ~ FC >= "\u00B1"*.(fccut))),
                  legendLabSize = 12,
                  legendIconSize = 5.0,
                  drawConnectors = T,
                  widthConnectors = 0.25,
                  maxoverlapsConnectors = Inf,
                  arrowheads = F,
                  gridlines.minor = F,
                  gridlines.major = F) +
  labs(subtitle = "",
       y = expression( -Log[10] ~ p)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.75, label = "post-transplant",
           size = 9, color = t.pal[2]) +
  annotate("rect", xmin = fccut + 0.05, xmax = Inf, ymin = -log(pcut, base = 10),
           ymax = Inf,
           alpha = .1, fill = t.pal[2]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.75, label = "pre-transplant",
           size = 9, color = t.pal[1]) +
  annotate("rect", xmin = -(fccut + 0.05), xmax = -Inf, ymin = -log(pcut, base = 10),
           ymax = Inf,
           alpha = .1, fill = t.pal[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

qual_tp_volcano

ggsave(plot = qual_tp_volcano,
       paste0("figures/Fig1_C.volcano.pvals.ori.pdf"),
       width = 11, height = 6.05, units = "in")
