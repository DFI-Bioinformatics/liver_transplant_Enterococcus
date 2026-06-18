library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(gt)

# load meta ---------------------------------------------------------------

lt.sams <- read_csv("data/lt.meta_rawMetabolomics.csv")

demo <- read_csv("data/demographic.clinical.abx.137sub.complete.csv")

## color palette -----

source("code/0.pals.R")

u.meta <- read_csv("data/paired.umap_cluster.2025-03-12.csv") |> 
  mutate(cl = if_else(kcluster == 1, "Cluster D", "Cluster E"),
         kcluster = factor(kcluster)) |> 
  select(-c(paired, timepoint)) |> 
  distinct()

cluster.pat <- lt.sams |> 
  left_join(u.meta) |> 
  select(studyid, kcluster) |> 
  group_by(studyid) |> 
  summarise(Cluster = if_else(any(kcluster == 2),
                              "Cluster E", "Cluster D"))

km.df <- lt.sams |> 
  distinct(studyid, icu_los, death_date, event_date, death) |> 
  left_join(demo |> 
              select(studyid, `Meld-Na`, `Clinical Status`,
                     Dialysis, `Antibiotics Exposure >= 5 Days`)) |> 
  # filter(trans_date == event_date)
  mutate(t2e = as.numeric(death_date - event_date),
         event = if_else(is.na(t2e) | t2e > 365, 0, 1)) |> 
  replace_na(list(t2e = 365)) |> 
  left_join(cluster.pat)

cluster <- km.df$Cluster
mortality <- factor(km.df$death, levels = c("Yes", "No"), labels = c("Death", "Survival"))

s.fit <- survfit(Surv(t2e, event) ~ Cluster, data = km.df) 

cs.fit <- coxph(Surv(t2e, event) ~ Cluster, data = km.df)

c.fit <- coxph(Surv(t2e, event) ~ Cluster + icu_los + `Meld-Na` + Dialysis + `Antibiotics Exposure >= 5 Days`, data = km.df)

hr <- cs.fit |> 
  tidy(conf.int = TRUE, exponentiate = TRUE) |> 
  filter(term == "ClusterCluster E") |> 
  select(estimate, starts_with("conf")) |> 
  unlist()

summary(c.fit)$coefficients |> 
  as.data.frame() |> 
  rownames_to_column() |> 
  left_join(exp(confint(c.fit)) |> 
              as.data.frame() |> 
              rownames_to_column()) |> 
  as_tibble() |> 
  select(rowname, `exp(coef)`, `Pr(>|z|)`,`2.5 %`, `97.5 %`) |> 
  write_csv("data/coxph_cluster.csv")


# figure 3 ----------------------------------------------------------------

## Panel B -----------------------------------------------------------------

lvl = "family"

masslin <- read_tsv(paste0("data/clusterMaaslin/",lvl,"/significant_results.tsv"))

m.gg <-  masslin |> 
  filter(abs(coef) >= 4, metadata == "cl") |> 
  mutate(group = if_else(coef > 0, "Cluster E", "Cluster D")) |> 
  ggplot(aes(coef, reorder(feature, coef), fill = group)) +
  geom_col() +
  scale_fill_manual(values = c.pal) +
  labs(x = "Coefficient", y = "Family", fill = "") +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "top")

m.gg


## Panel C: volcano -----------------------------------------------------------------

qual.imp <- read_csv("data/peri_post.qual.imput.2026-05-12.ori.csv")

qual.imp.long <-  qual.imp |> 
  gather("compound","value", -sample_id) |> 
  as_tibble()

qual_v_df <- qual.imp.long |> 
  left_join(u.meta |> 
              select(shotgun_seq_id, cl) |> 
              left_join(lt.sams |> 
                          select(shotgun_seq_id, 
                                 sample_id = metabolomics_id, 
                                 study_id = studyid)) )

qual_tp_df <- qual_v_df |> 
  mutate(
    cl = factor(cl, levels = c("Cluster D", "Cluster E")),
    study_id = factor(study_id)
  ) |> 
  group_by(compound) %>%
  # Collapse the rest of the data into a "list-column" for each taxon
  nest() %>%
  mutate(
    # Because lmerTest is loaded, this lmer() now calculates degrees of freedom!
    model = map(data, ~ lmer(value ~ cl + (1 | study_id), data = .x)),
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
  group_by(compound, cl) |> 
  summarise(ave = mean(value)) |> 
  mutate(ave = if_else(ave == 0, 0.0001, ave)) |> 
  spread(cl, ave) |> 
  mutate(log2fc = log2(`Cluster E`/`Cluster D`))

qual_tp_tot <- left_join(qual_tp_fc, qual_tp_df) %>%
  column_to_rownames(var = "compound")

#### p-adjusted -------
set.seed(123456)

xylims = ceiling(max(abs(qual_tp_tot$log2fc)))
plims <- ceiling(max(abs(log10(qual_tp_tot$p.adj))))

pcut = 0.05
fccut = 1.5

library(EnhancedVolcano)

qual_tp_volcano <-
  EnhancedVolcano(qual_tp_tot,
                  lab = rownames(qual_tp_tot),
                  title = "",
                  y = "p.adj",
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
       y = expression( -Log[10] ~ q)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.75, label = "Cluster E",
           size = 9, color = c.pal[2]) +
  annotate("rect", xmin = fccut + 0.05, xmax = Inf, ymin = -log(pcut, base = 10),
           ymax = Inf,
           alpha = .1, fill = c.pal[2]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.75, label = "Cluster D",
           size = 9, color = c.pal[1]) +
  annotate("rect", xmin = -(fccut + 0.05), xmax = -Inf, ymin = -log(pcut, base = 10),
           ymax = Inf,
           alpha = .1, fill = c.pal[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

qual_tp_volcano

## Panel A ----------------------------------------------

u.df <- u.meta |> 
  left_join(lt.sams |> 
              select(shotgun_seq_id, timepoint3))

paired.cluster <- u.meta |> 
  left_join(lt.sams |> 
              select(shotgun_seq_id, studyid, timepoint)) |> 
  add_count(studyid, name = "sample_n") |>  
  filter(sample_n == 2) |> 
  select(studyid, timepoint, kcluster) |> 
  spread(timepoint, kcluster) |> 
  mutate(cluster_pos = case_when(
    Early == 2 & Late == 2 ~ "Stayed in Cluster E",
    Early == 1 & Late == 1 ~ "Stayed in Cluster D",
    Early == 2 & Late == 1 ~ "Moved from Cluster E to D",
    TRUE ~ "Moved from Cluster D to E",
  ))

arow.df <- u.meta |> 
  left_join(lt.sams |> 
              select(shotgun_seq_id, studyid, timepoint)) |> 
  add_count(studyid, name = "sample_n") |>  
  filter(sample_n == 2)  |> 
  left_join(paired.cluster) |> 
  filter(cluster_pos == "Moved from Cluster D to E")

a.gg <- u.df |> 
  ggplot(aes(V1, V2)) +
  geom_point(aes(fill = timepoint3), size = 3.5, alpha = 0.85, shape = 21) +
  scale_fill_manual(values = t2.pal) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 6.5),
        legend.position="top") +
  geom_path(data = arow.df,
            aes(group = studyid, linetype = cluster_pos),
            arrow = arrow(length=unit(0.30,"cm"),
                          type = "closed"),
            alpha = 0.65, show.legend = FALSE) +
  labs(x = "UMAP1", y = "UMAP2", fill = "")

a.gg

u.df |> 
  ggplot(aes(V1, V2)) +
  geom_point(aes(fill = cl), size = 3.5, alpha = 0.85, shape = 21) +
  scale_fill_manual(values = c.pal) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 6.5),
        legend.position="top") +
  geom_path(data = arow.df,
            aes(group = studyid, linetype = cluster_pos),
            arrow = arrow(length=unit(0.30,"cm"),
                          type = "closed"),
            alpha = 0.65, show.legend = FALSE) +
  labs(x = "UMAP1", y = "UMAP2", fill = "")

fig3 <- (a.gg + m.gg ) /  qual_tp_volcano +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A")

ggsave("figures/Fig3_cluster.pdf", 
       fig3,
       height = 10, width = 12.5)

# figure 4 ----------------------------------------------------------------

gtbl <- table(cluster, mortality) |> 
  as.data.frame() |> 
  spread(mortality, Freq) |> 
  mutate(`Mortality Rate` = signif(Death/ (Death + Survival) * 100, digits = 2)) |> 
  gt()

names(c.pal) <- NULL

suv.gg <- ggsurvplot(s.fit, 
                     data = km.df,
                     pval = T,
                     pval.coord = c(5, .85),
                     pval.size = 6.5,
                     palette = c.pal,
                     xlab = "Time (days)",
                     conf.int = T,
                     risk.table = "nrisk_cumevents",
                     risk.table.pos = "in",
                     risk.table.fontsize = 5.5,
                     tables.y.text = FALSE,
                     tables.col = "strata",
                     legend.title = "",
                     legend.labs = c("Cluster D", "Cluster E"),
                     xlim = c(0, 365),
                     ylim = c(0.5, 1),
                     break.x.by = 60) 

# s.tbl <- suv.gg$table +
#   theme(legend.position = "none") +
#   labs(title = "", y = "")

s.plt <- suv.gg$plot +
  labs(color = "", fill = "")  + 
  annotate("text",x = 5, y = 0.75, label = paste0("HR = ",round(hr[["estimate"]],2),
                                                  "(",round(hr[["conf.low"]],2),",",
                                                  round(hr[["conf.high"]],2),")"), 
           size =6.5, hjust = 0) +
  inset_element(gtbl, 0.5, 0.5, 0.3, 0.3) 

s.plt

ggsave("figures/Fig4_KM.cluster.pdf", height = 6, width = 7.5)

# suppFig4 ----------------------------------------------------------------

source("https://raw.githubusercontent.com/yingeddi2008/DFIutility/master/getRdpPal.R")

mpa <- read_csv("data/mpa.june23.csv")

mpa.plt <- mpa |> 
  filter(seq_id %in% lt.sams$shotgun_seq_id) |> 
  add_count(seq_id, 
            wt = clean_relative_abundance, 
            name = "totalAbd") %>%
  mutate(pctseqs = clean_relative_abundance/totalAbd,
         genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

taxpal <- getRdpPal(mpa.plt)

mpa.gg <-  mpa.plt %>%
  group_by(seq_id, Kingdom,Phylum,Class,Order,Family, Genus, genLab) %>%
  summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() %>% # view
  left_join(u.meta |> 
              select(seq_id = shotgun_seq_id,
                     cl)) |> 
  ggplot(aes(x=seq_id ,y=pctseqs)) +
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
  facet_grid(. ~ cl, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0.001,0.001))

mpa.gg

ggsave("figures/SuppFig4_cluster.pdf",
       height = 4.95, width = 10.5)
