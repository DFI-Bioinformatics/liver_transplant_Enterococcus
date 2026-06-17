library(tidyverse)
library(umap)
library(vegan)

# load meta ---------------------------------------------------------------

lt.sams <- read_csv("data/lt.meta_rawMetabolomics.csv")

## color palette -----

source("code/0.pals.R")

# load metaphlan 4 --------------------------------------------------------

mpa <- read_csv("data/mpa.june23.csv")

mpa.plt <- mpa |> 
  add_count(seq_id, 
            wt = clean_relative_abundance, 
            name = "totalAbd") %>%
  mutate(pctseqs = clean_relative_abundance/totalAbd,
         genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

mpa_mat <- mpa.plt |> 
  filter(seq_id %in% lt.sams$shotgun_seq_id) |> 
  select(seq_id, taxid, pctseqs) |> 
  pivot_wider(names_from = taxid, 
              values_from = pctseqs,
              values_fn = sum,
              values_fill = 0) |> 
  column_to_rownames(var = "seq_id")

mpa_mat_trans <- round(mpa_mat*1000000)

# repeat umap 10 times ----------------------------------------------------

library(umap)

cus.umap.config <- umap.defaults
cus.umap.config$n_neighbors <- 50
cus.umap.config$min_dist <- 0.2

sds <- sample(1:100000, 12)

ucoord.list <- list()

for (s in sds){
  
  cus.umap.config$random_state <- s
  
  ucoord = umap(mpa_mat_trans, config = cus.umap.config)
  
  u.meta <- ucoord$layout |> 
    as.data.frame() %>% 
    rownames_to_column(var = "shotgun_seq_id") %>%
    as_tibble() |> 
    full_join(lt.sams |> 
                select(shotgun_seq_id, timepoint))
  
  km <- kmeans(ucoord$layout, centers = 2, nstart = 5)
  
  km_df = tibble(kcluster = km$cluster,seq_id = names(km$cluster))
  
  ucoord.list[[as.character(s)]] <- u.meta |> 
    left_join(km_df, by = c("shotgun_seq_id"="seq_id")) |> 
    mutate(seed = s)
  
}

ucoord.com <- ucoord.list |> 
  bind_rows() |> 
  mutate(seed = factor(paste0("seed:", seed))) 

u.meta <- read_csv("data/paired.umap_cluster.2025-03-12.csv") |> 
  mutate(cl = if_else(kcluster == 1, "Cluster D", "Cluster E"),
         kcluster = factor(kcluster)) |> 
  select(-c(paired, timepoint)) |> 
  distinct()

ucoord.com |> 
  left_join(u.meta |> 
              select(-c(V1,V2,kcluster))) |> 
  ggplot(aes(V1, V2, fill = cl)) +
  geom_point(alpha = 0.65, shape = 21) +
  theme_bw() +
  scale_fill_manual(values = c( "#9344d5", "#36a02d")) +
  facet_wrap("seed", scales = "free", ncol = 3) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(fill = "")

ggsave("figures/SuppFig3_umap12.byCluster.pdf", height = 11, width = 9)

## enterotype 3k -----------------------------------------------------

fk3 <- read_tsv("data/EnterotypeAssignment_fuzzyClustering_k3_tpt.tsv") |> 
  mutate(sampleid = gsub("\\.","_", sample_id)) |> 
  left_join(lt.sams |> 
              select(sampleid, timepoint))

fk3 |> 
  select(sampleid, Enterotype_Dysbiosis_Score) |> 
  left_join(lt.sams |> 
              select(shotgun_seq_id, sampleid)) |> 
  left_join(u.meta) |> 
  ggplot(aes(cl, Enterotype_Dysbiosis_Score, fill = cl )) +
  geom_boxplot(outlier.shape = NA, alpha = 0.65) +
  geom_jitter(shape = 21, alpha = 0.85) +
  scale_fill_manual(values = c(  "#9344d5", "#36a02d")) +
  theme_bw() +
  ggpubr::stat_compare_means() +
  theme(legend.position = "none") +
  labs(y = "Enterotype Dysbiosis Score (EDS)", x = "")

ggsave("figures/SuppFig5_enterotype3k_cluster.EDS.bwplt.byCluster.pdf", 
       height = 5.5, width = 6)
