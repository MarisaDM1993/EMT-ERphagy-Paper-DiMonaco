# =============================================================================
# Gaussian mixture clustering of whole-cell-lysate proteomics, and GO
# enrichment of each cluster.
#
# Paper   : EMT Induces FAM134A-mediated ER-phagy to Regulate Procollagen
#           Secretion. Raw proteomics: ProteomeXchange/PRIDE PXD073014.
#
# Produces: Fig 4C       - CC enrichment dot plots for clusters of interest
#           Supp Fig 4A  - standardised expression profile of all 15 clusters
#           Supp Table 2 - sheets ClusterData, GOenrich_full,
#                          GOenrich_simplified
#
# Input   : data_from_mass_spec_facility/CTRL.ConA_vs_CTRL.CTRL.csv
#           Facility output from DIA-NN. Only the per-sample intensity
#           columns are used; the facility's log2FC and p-value columns are
#           dropped (see columns_to_remove). The log2FC sign error corrected
#           in proteomics_differential_expression_and_GO.R therefore does
#           not affect anything in this script.
#
# Writes  : rawdata_long.csv, output_clustered_data.csv,
#           go_enrichment_output/{full,simplified}/cluster<N>_GO_<ONT>.csv
#
# Run     : top to bottom, from the directory holding
#           data_from_mass_spec_facility/. Independent of the differential
#           expression script; neither reads the other's output.
# =============================================================================

library(tidyverse)
library(mclust)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

################################################################################
# ADD DATA FROM MS FACILITY
################################################################################

facility_rawdata <- read_csv(file="data_from_mass_spec_facility/CTRL.ConA_vs_CTRL.CTRL.csv")

# store the names of raw value columns and imputation columns
raw_cols <- grep("^[A-Za-z]+_[A-Za-z]+\\.\\d+$", colnames(facility_rawdata), value = TRUE)
imp_cols <- paste0("Imputation_", raw_cols)

# remove columsn which are not needed later (store id columns) 
columns_to_remove <- c(raw_cols, imp_cols, 
                       "-log10 adj.P.Value", "-log10 P.Value",
                       "Comparison", "log2FC", "Significance")

id_cols <- facility_rawdata |>
  select(-all_of(columns_to_remove))

# tidy up the raw value data (pivot_longer)

raw_long <- facility_rawdata |>
  dplyr::select(Protein.Group, all_of(raw_cols)) |>
  pivot_longer(cols=all_of(raw_cols), names_to="Sample", values_to="Intensity")

# tidy up the imputation data (pivot_longer)
imp_long <- facility_rawdata |>
  dplyr::select(Protein.Group, all_of(imp_cols)) |>
  pivot_longer(cols=all_of(imp_cols), names_to="Sample", values_to="Imputation") |>
  mutate(Sample=str_remove(Sample, "^Imputation_"))

# join up raw_long and imp_long
long_data <- left_join(raw_long, imp_long, by=c("Protein.Group", "Sample"))

# create the finalised raw dataset
tidy_rawdata <- long_data |>
  left_join(id_cols, by="Protein.Group") |>
  separate(Sample, into=c("Treatment", "Replicate"), sep = "\\.", convert = TRUE) |>
  drop_na(Intensity) |>
  mutate(Imputation=if_else(is.na(Imputation), 0, Imputation))

tidy_rawdata |>head()

# write to disk
write_csv(tidy_rawdata, file="rawdata_long.csv")

# remove unnecessary objects
remove(tidy_rawdata, facility_rawdata, imp_long, raw_long, long_data, imp_cols, raw_cols)

################################################################################
#Import raw data and prepare it for clustering

df_long <- read_csv("rawdata_long.csv")

# summarise intensity across replicates by median
df_median <- df_long |>
  group_by(Protein.Group, Treatment) |>
  summarise(median_intensity=median(Intensity), .groups = "drop")

# change to wide format
df_median_wide <- df_median |>
  pivot_wider(
    names_from = Treatment,
    values_from = c(median_intensity)
  )

# z-score standardisation
df_scaled_wide <- df_median_wide |>
  column_to_rownames("Protein.Group") |> # extract numerical matrix
  as.matrix() |>
  t() |> # transpose
  scale(center=TRUE, scale=TRUE) |> # z-score standardisation 
  t() |> # re-transpose
  as.data.frame() |> 
  rownames_to_column("Protein.Group") # return rownames to columns

################################################################################
# perform clustering

# Fixes the Mclust initialisation. This matters for more than numerical
# reproducibility: GMM cluster labels are arbitrary, and the manuscript
# refers to specific clusters by number (Cluster 8, Cluster 11). A different
# seed, or a different mclust version, may permute the labels even when the
# partition is unchanged.
set.seed(123)

# model fit
gmm_model <- Mclust(df_scaled_wide[-1], G=15, modelNames = "EII")
gmm_model |> summary()

# add clustering data and post. prob. to original df
df_clustered_wide <- df_scaled_wide |>
  mutate(Cluster.GMM=gmm_model$classification,
         Post.Prob=apply(gmm_model$z, 1, max))

# change back to long format
df_clustered <- df_clustered_wide |>
  pivot_longer(
    cols = -c(Protein.Group, Cluster.GMM, Post.Prob),
    names_to = "Treatment",
    values_to = "Std.Intensity") |>
  left_join(df_long |>
              dplyr::select(Protein.Group, Protein.Names, Genes) |> # add id columns
              distinct(Protein.Group, .keep_all = TRUE),
            by = "Protein.Group") |>
  mutate(Treatment=factor(Treatment, levels=c("CTRL_CTRL", "CTRL_ConA", "TGFb_CTRL", "TGFb_ConA")))

# create labels for the plot

cluster_counts <- df_clustered |>
  mutate(Protein.Group = as.character(Protein.Group)) |>
  distinct(Protein.Group, Cluster.GMM) |>
  group_by(Cluster.GMM) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(label = paste0("Cluster ", Cluster.GMM, " || n = ", n))

labels_named <- setNames(cluster_counts$label, cluster_counts$Cluster.GMM)

# calculate mean profile for each cluster
cluster_means <- df_clustered |>
  group_by(Cluster.GMM, Treatment) |>
  summarise(
    mean_intensity = mean(Std.Intensity, na.rm = TRUE),
    .groups = "drop"
  )

# plot the clusters (with post. probability)
cluster_plot <- ggplot(df_clustered, aes(x = Treatment, y = Std.Intensity, group = Protein.Group)) +
  ggrastr::rasterise(geom_line(alpha = 0.1, aes(colour=Post.Prob)), dpi=600) +
  geom_line(data = cluster_means, aes(x = Treatment, y = mean_intensity, group = Cluster.GMM),
            colour = "black", linewidth = 0.75) +
  facet_wrap(~ Cluster.GMM, nrow = 3, ncol = 5,labeller = labeller(Cluster.GMM = labels_named)) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(-1.5, 2, 1.5), limits=c(-2,2)) +
  scale_colour_gradient(low="navy", high="red", name = "Posterior\nProb.",
                        limits=c(min(df_clustered$Post.Prob), 1)) +
  labs(title = "Cluster profiles",
       x = "Treatment",
       y = "Std intensity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster_plot 

# write clustered data to disk
write_csv(df_clustered, "output_clustered_data.csv")

# re-import data if needed 
# df_clustered <- read_csv("output_clustered_data.csv") |>
  # mutate(Treatment=factor(Treatment, levels=c("CTRL_CTRL", "CTRL_ConA", "TGFb_CTRL", "TGFb_ConA")))


################################################################################
# proceed with the GO enrichment analysis

dir.create("go_enrichment_output/full", showWarnings = FALSE, recursive = TRUE)
dir.create("go_enrichment_output/simplified", showWarnings = FALSE, recursive = TRUE)

universe_genes <- df_clustered |>
  filter(!grepl("cRAP", Protein.Group),!is.na(Genes), Genes != "") |>
  pull(Genes) |>
  strsplit(";") |>
  unlist() |>
  unique()

for (cluster_no in 1:15) {
  
  cluster_genes <- df_clustered |>
    filter(Cluster.GMM == cluster_no) |>
    pull(Genes) |>
    strsplit(";") |>
    unlist() |>
    unique()
  
  go_results_full <- enrichGO(gene = cluster_genes, universe = universe_genes,
                     OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")
  
  # skip if empty
  if (is.null(go_results_full) || nrow(as.data.frame(go_results_full)) == 0) next
  
  # save results
  go_full <- as.data.frame(go_results_full)
  
  # write to disk
  for (ont in c("BP", "CC", "MF")) {
    
    write.csv(
      subset(go_full, ONTOLOGY == ont),
      paste0(
        "go_enrichment_output/full/cluster",
        cluster_no, "_GO_", ont, ".csv"
      ),
      row.names = FALSE
    )
  }
  
# simplify and save the reduced term list
  go_results_simp <- simplify(
    go_results_full,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )
  
  go_simp <- as.data.frame(go_results_simp)
  
  for (ont in c("BP", "CC", "MF")) {
    
    write.csv(
      subset(go_simp, ONTOLOGY == ont),
      paste0(
        "go_enrichment_output/simplified/cluster",
        cluster_no, "_GO_", ont, ".csv"
      ),
      row.names = FALSE
    )
  }
}

################################################################################
# re-import the data as single dfs

combine_results <- function(folder_path) {
  
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  map_dfr(files, \(f) {
    read_csv(
      f,
      col_types = cols(.default = col_character()),
      show_col_types = FALSE
    ) |>
      mutate(source_file = basename(f))
  })
}

df_full <- combine_results("go_enrichment_output/full") |>
  mutate(
    across(c(RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, Count), as.numeric),
    cluster_no = str_extract(source_file, "(?<=cluster)\\d+") |>
      as.integer()
    )

df_simp <- combine_results("go_enrichment_output/simplified") |>
  mutate(
    across(c(RichFactor, FoldEnrichment, zScore, pvalue, p.adjust, qvalue, Count), as.numeric),
    cluster_no = str_extract(source_file, "(?<=cluster)\\d+") |>
      as.integer()
  )

# write_csv(df_full, file="go_enrichment_output/GO_full.csv")
# write_csv(df_simp, file="go_enrichment_output/GO_simplified.csv")

# plot_cluster_go_dot, version 1 of 3. SUPERSEDED: redefined below before any
# call is made, so this definition never runs. Kept for the record.
plot_cluster_go_dot <- function(cluster_id, data, ontology, aspectratio, top){
  
  data |>
    filter(`ONTOLOGY` %in% ontology) |>
    filter(cluster_no == cluster_id) |>
    group_by(`ONTOLOGY`) |>
    slice_min(p.adjust, n = top) |> # filter the top n most significant terms
    ungroup() |>
    mutate(Description = forcats::fct_reorder(Description, FoldEnrichment)) |>
    
    ggplot(aes(x = FoldEnrichment, y = Description)) +
    geom_point(shape = 21, stroke = 0.5, fill="red", 
               aes(size = -log10(p.adjust))) +
    facet_wrap(~ONTOLOGY, scales = "free") +
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title = paste0("Cluster ", cluster_id),
      x = "Fold enrichment",
      y = NULL
    ) +
    theme(strip.text = element_text(size = 5),
          axis.text.x = element_text(size=6),
          aspect.ratio = aspectratio)
}


# plot_cluster_go_dot, version 2 of 3: GO ID on the y axis. This one serves the
# calls immediately below it.
plot_cluster_go_dot <- function(cluster_id, data, ontology, aspectratio, top){
  
  data |>
    filter(`ONTOLOGY` %in% ontology) |>
    filter(cluster_no == cluster_id) |>
    group_by(`ONTOLOGY`) |>
    slice_min(p.adjust, n = top) |> # filter the top n most significant terms
    ungroup() |>
    mutate(ID = forcats::fct_reorder(ID, FoldEnrichment)) |>
    
    ggplot(aes(x = FoldEnrichment, y = ID)) +
    geom_point(shape = 21, stroke = 0.5, 
               aes(size = -log10(p.adjust), fill=Count)) +
    facet_wrap(~ONTOLOGY, scales = "free") +
    scale_fill_continuous(palette = c("#FEE0D2", "#FC9272", "#DE2D26"))+
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title = paste0("Cluster ", cluster_id),
      x = "Fold enrichment",
      y = NULL
    ) +
    theme(strip.text = element_text(size = 5),
          axis.text.x = element_text(size=6),
          aspect.ratio = aspectratio)
}

# plotting clusters
c8_cc <- plot_cluster_go_dot(8, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c9_cc <- plot_cluster_go_dot(9, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c14_cc <- plot_cluster_go_dot(14, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c15_cc <- plot_cluster_go_dot(15, df_simp, c("CC"), aspectratio=1/0.7,top=20)

gridExtra::grid.arrange(c8_cc, c9_cc, c14_cc, c15_cc, ncol=2)


# plot_cluster_go_dot, version 3 of 3: GO term Description on the y axis. This
# is the version whose output appears in Fig 4C, since the calls below reassign
# c8_cc/c9_cc/c14_cc/c15_cc.

df_simp

plot_cluster_go_dot <- function(cluster_id, data, ontology, aspectratio, top){
  
  data |>
    filter(`ONTOLOGY` %in% ontology) |>
    filter(cluster_no == cluster_id) |>
    group_by(`ONTOLOGY`) |>
    slice_min(p.adjust, n = top) |> # filter the top n most significant terms
    ungroup() |>
    mutate(Description = forcats::fct_reorder(Description, FoldEnrichment)) |>
    
    ggplot(aes(x = FoldEnrichment, y = Description)) +
    geom_point(shape = 21, stroke = 0.5, 
               aes(size = -log10(p.adjust), fill=Count)) +
    facet_wrap(~ONTOLOGY, scales = "free") +
    scale_fill_continuous(palette = c("#FEE0D2", "#FC9272", "#DE2D26"))+
    scale_size_continuous(range = c(1, 6)) +
    labs(
      title = paste0("Cluster ", cluster_id),
      x = "Fold enrichment",
      y = NULL
    ) +
    theme(strip.text = element_text(size = 5),
          axis.text.x = element_text(size=6),
          aspect.ratio = aspectratio)
}

c8_cc <- plot_cluster_go_dot(8, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c9_cc <- plot_cluster_go_dot(9, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c14_cc <- plot_cluster_go_dot(14, df_simp, c("CC"), aspectratio=1/0.7, top=20)
c15_cc <- plot_cluster_go_dot(15, df_simp, c("CC"), aspectratio=1/0.7,top=20)

gridExtra::grid.arrange(c8_cc, c9_cc, c14_cc, c15_cc, ncol=2)

# =============================================================================
# Record package versions alongside the outputs.
# =============================================================================
writeLines(capture.output(sessionInfo()), "sessionInfo_clustering_and_GO.txt")
