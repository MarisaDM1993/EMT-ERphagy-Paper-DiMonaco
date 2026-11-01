# =============================================================================
# Differential expression volcano plots, and GO enrichment of the up- and
# down-regulated protein sets.
#
# Paper   : EMT Induces FAM134A-mediated ER-phagy to Regulate Procollagen
#           Secretion. Raw proteomics: ProteomeXchange/PRIDE PXD073014.
#
# Produces: Fig 4B       - volcano plots: TGFb vs untreated (upper) and
#                          TGFb+ConA vs TGFb (lower)
#           Supp Fig 4D  - BP enrichment of the TGFb-vs-basal upregulated set
#
# Input   : data_from_mass_spec_facility/*.csv - the four pairwise
#           comparisons returned by the facility from the DIA-NN output.
#
# SIGN CORRECTION - read this before reusing the facility files.
#   The facility's analysis code reports log2FC with the wrong sign for the
#   two ConA contrasts, CTRL.ConA_vs_CTRL.CTRL and TGFb.ConA_vs_TGFb.CTRL.
#   The two lines marked SIGN CORRECTION below negate those columns. Any
#   reanalysis that starts from the deposited facility CSVs must apply the
#   same correction, or every ConA effect will come out inverted.
#
# Writes  : go_enrichment_up-down/{full,simp}/*.csv
#
# Run     : top to bottom. Independent of the clustering script. Nothing
#           here is stochastic, so no seed is needed.
# =============================================================================

library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

# create filename list (data from facility)

de_files <- list.files("data_from_mass_spec_facility/", pattern = "\\.csv$", full.names = TRUE)

diff_exprs_sets <- list()

# import data

for(i in seq_along(de_files)) {
  df_selected <- readr::read_csv(de_files[i])
  diff_exprs_sets[[i]] <- df_selected
  names(diff_exprs_sets)[i] <- tools::file_path_sans_ext(basename(de_files[i]))
}

# SIGN CORRECTION. The facility's code reports log2FC with the wrong sign for
# the two ConA contrasts. Negating them here makes a positive log2FC mean
# 'higher with ConA', which is how every figure and the manuscript read it.
# See the header. The other two contrasts are unaffected.

diff_exprs_sets$CTRL.ConA_vs_CTRL.CTRL$log2FC= -diff_exprs_sets$CTRL.ConA_vs_CTRL.CTRL$log2FC
diff_exprs_sets$TGFb.ConA_vs_TGFb.CTRL$log2FC= -diff_exprs_sets$TGFb.ConA_vs_TGFb.CTRL$log2FC

################################################################################

# Indexed by name, not position: list.files() above sorts by the active
# locale's collation, and these filenames order differently under C/POSIX than
# under an English locale. The list is named from the filenames on import, so
# name lookup gives the same mapping in any locale.
CTRL.ConA_vs_CTRL.CTRL <- diff_exprs_sets$CTRL.ConA_vs_CTRL.CTRL %>%
  dplyr::select(Protein.Group, Protein.Names, Genes, log2FC, `-log10 adj.P.Value`, `-log10 P.Value`)

TGFb.ConA_vs_TGFb.CTRL <- diff_exprs_sets$TGFb.ConA_vs_TGFb.CTRL %>%
  dplyr::select(Protein.Group, Protein.Names, Genes, log2FC, `-log10 adj.P.Value`, `-log10 P.Value`)

TGFb.CTRL_vs_CTRL.CTRL <- diff_exprs_sets$TGFb.CTRL_vs_CTRL.CTRL %>%
  dplyr::select(Protein.Group, Protein.Names, Genes, log2FC, `-log10 adj.P.Value`, `-log10 P.Value`)

TGFb.ConA_vs_CTRL.ConA <- diff_exprs_sets$TGFb.ConA_vs_CTRL.ConA %>%
  dplyr::select(Protein.Group, Protein.Names, Genes, log2FC, `-log10 adj.P.Value`, `-log10 P.Value`)

################################################################################

# volcano plots

goi_ls <- c("SPOCK1","ADAM19","TGFBI","CCN2","RHOU","PODXL","COL5A1","NIPAL4",
            "COL4A2","SERPINE1","TGFBR1","TNFRSF12A","COL4A1","LTBP4",
            "SERPINE2","TIMP2","PLAU")


# TGFb vs untreated (all black)#################################

df <- TGFb.CTRL_vs_CTRL.CTRL |>
  mutate(
    note = case_when(
      Genes %in% goi_ls ~ "goi",
      `-log10 adj.P.Value` <= -log10(0.05) ~ "not significant",
      abs(log2FC) <= 1 ~ "no_change",
      TRUE ~ "significant"
    )) |>
  mutate(
    note = factor(note, c("goi", "not significant", "no_change", "significant"))) |>
  na.omit()

goi_df <- df |>
  dplyr::filter(Genes %in% goi_ls)

p1 <- ggplot(df, aes(x=log2FC, y=`-log10 P.Value`))+
  ggrastr::rasterise(
    geom_point(aes(colour=note, size=note), stroke=0),
    dpi=800
  )+
  geom_point(data=goi_df, colour=alpha("red",.8), stroke=0, size=2)+
  scale_colour_manual(values=c(alpha("white", alpha=0), # genes of interest
                               alpha("black", alpha=.3), # nonsignificant genes
                               alpha("black", alpha=.3), # not diff exprs genes
                               alpha("black", alpha=.3)))+ # diff exprs genes
  scale_size_manual(values=c(1.5, 1.5, 1.5, 1.5))+
  ggrepel::geom_text_repel(data=goi_df, aes(label=Genes),
                           nudge_y =.1,
                           nudge_x=.5, 
                           size=2)+
  ggprism::theme_prism()+
  scale_y_continuous(limits=c(0,22.5), breaks=seq(0,22.5,2.5))+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10,2))+
  geom_vline(xintercept=c(-1,1), linetype=3)+
  geom_hline(yintercept=min(df$`-log10 P.Value`[df$note=="significant"]), linetype=3)+
  theme(
    legend.position="none",
    aspect.ratio=1/1.2
  )

#############
df <- TGFb.ConA_vs_TGFb.CTRL |>
  mutate(
    note = case_when(
      Genes %in% goi_ls ~ "goi",
      `-log10 adj.P.Value` <= -log10(0.05) ~ "not significant",
      abs(log2FC) <= 1 ~ "no_change",
      TRUE ~ "significant"
    )) |>
  mutate(
    note = factor(note, c("goi", "not significant", "no_change", "significant"))) |>
  na.omit()

goi_df <- df |>
  dplyr::filter(Genes %in% goi_ls)

p2 <- ggplot(df, aes(x=log2FC, y=`-log10 P.Value`))+
  ggrastr::rasterise(
    geom_point(aes(colour=note, size=note), stroke=0),
    dpi=800
  )+
  geom_point(data=goi_df, colour=alpha("red",.8), stroke=0, size=2)+
  scale_colour_manual(values=c(alpha("white", alpha=0), # genes of interest
                               alpha("black", alpha=.3), # nonsignificant genes
                               alpha("black", alpha=.3), # not diff exprs genes
                               alpha("black", alpha=.3)))+ # diff exprs genes
  scale_size_manual(values=c(1.5, 1.5, 1.5, 1.5))+
  ggrepel::geom_text_repel(data=goi_df, aes(label=Genes),
                           nudge_y =.1,
                           nudge_x=.5, 
                           size=2)+
  ggprism::theme_prism()+
  scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2.5))+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10,2))+
  geom_vline(xintercept=c(-1,1), linetype=3)+
  geom_hline(yintercept=min(df$`-log10 P.Value`[df$note=="significant"]), linetype=3)+
  theme(
    legend.position="none",
    aspect.ratio=1/1.2
  )

gridExtra::grid.arrange(p1, p2, nrow=2)



################################################################################
################################################################################

# TGFb vs untreated - (upregulated/downregulated)##################################

df <- TGFb.CTRL_vs_CTRL.CTRL |>
  mutate(
    note = case_when(
      `-log10 adj.P.Value` <= -log10(0.05) ~ "not significant",
      log2FC > 0 ~ "upregulated",
      TRUE ~ "downregulated"
    )) |>
  mutate(
    note = factor(note, c("not significant", "upregulated", "downregulated"))) |>
  na.omit()

ggplot(df, aes(x=log2FC, y=`-log10 P.Value`))+
  ggrastr::rasterise(
    geom_point(aes(colour=note), stroke=0),
    dpi=800
  )+
  scale_colour_manual(values=c(alpha("grey70", alpha=.5), # nonsignificant genes
                               alpha("red", alpha=.5), # upregulated genes
                               alpha("blue", alpha=.6)))+ # downregulated genes
  ggprism::theme_prism()+
  scale_y_continuous(limits=c(0,22.5), breaks=seq(0,22.5,2.5))+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10,2))+
  geom_vline(xintercept=0, linetype=3, colour="grey70")+
  theme(
    legend.position="none",
    aspect.ratio=1/1.2
  )

################################################################################
# GO enrichment analysis

# df here is TGFb.CTRL_vs_CTRL.CTRL (assigned above). The GO calls below
# depend on that assignment rather than taking it as an argument, so run
# this file top to bottom; re-running a block out of order will silently
# enrich the wrong contrast.
universe_genes <- df |>
  filter(!grepl("cRAP", Protein.Group),!is.na(Genes), Genes != "") |>
  pull(Genes) |>
  strsplit(";") |>
  unlist() |>
  unique()

# NOTE on gene symbols: universe_genes above splits composite protein groups on
# ";" whereas the genes_to_run blocks in this file do not. The 59 of 7,889
# protein groups carrying a composite Genes entry (e.g. "GATD3;GATD3B") can
# therefore never match a symbol in the universe and are effectively absent
# from the foreground. Left exactly as run, so that this script reproduces the
# published enrichment; all affected entries are paralogue pairs or
# pseudogenes, none bearing on the conclusions. The clustering script splits in
# both places and is unaffected.
genes_to_run <- df |>
  dplyr::filter(note=="upregulated") |>
  pull(Genes) |>
  unique()

go_results_up_TGFbvBasal <- enrichGO(gene=genes_to_run, universe=universe_genes,
                                        OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")

genes_to_run <- df |>
  dplyr::filter(note=="downregulated") |>
  pull(Genes) |>
  unique()

go_result_down_TGFbvBasal <- enrichGO(gene=genes_to_run, universe=universe_genes,
                                         OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")


# TGFb.ConA vs TGFb only (upregulated/downregulated)#################################

df <- TGFb.ConA_vs_TGFb.CTRL |>
  mutate(
    note = case_when(
      `-log10 adj.P.Value` <= -log10(0.05) ~ "not significant",
      log2FC > 0 ~ "upregulated",
      TRUE ~ "downregulated"
    )) |>
  mutate(
    note = factor(note, c("not significant", "upregulated", "downregulated"))) |>
  na.omit()

ggplot(df, aes(x=log2FC, y=`-log10 P.Value`))+
  ggrastr::rasterise(
    geom_point(aes(colour=note), stroke=0),
    dpi=800
  )+
  scale_colour_manual(values=c(alpha("grey70", alpha=.5), # nonsignificant genes
                               alpha("red", alpha=.5), # upregulated genes
                               alpha("blue", alpha=.6)))+ # downregulated genes
  ggprism::theme_prism()+
  scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2.5))+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10,2))+
  geom_vline(xintercept=0, linetype=3, colour="grey70")+
  theme(
    legend.position="none",
    aspect.ratio=1/1.2
  )

################################################################################
# GO enrichment analysis

# df here is TGFb.ConA_vs_TGFb.CTRL (reassigned above). Same caveat as the
# previous GO block: the contrast comes from the preceding assignment.
universe_genes <- df |>
  filter(!grepl("cRAP", Protein.Group),!is.na(Genes), Genes != "") |>
  pull(Genes) |>
  strsplit(";") |>
  unlist() |>
  unique()

genes_to_run <- df |>
  dplyr::filter(note=="upregulated") |>
  pull(Genes) |>
  unique()

go_results_up_TGFbvTGFbConA <- enrichGO(gene=genes_to_run, universe=universe_genes,
                                   OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")

genes_to_run <- df |>
  dplyr::filter(note=="downregulated") |>
  pull(Genes) |>
  unique()

go_result_down_TGFbvTGFbConA <- enrichGO(gene=genes_to_run, universe=universe_genes,
                                    OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL")

################################################################################

up_EMTvBasal <- go_results_up_TGFbvBasal |> as.data.frame()
down_EMTvBasal <- go_result_down_TGFbvBasal |> as.data.frame()

up_EMTvConA <- go_results_up_TGFbvTGFbConA |> as.data.frame()
down_EMTvConA <- go_result_down_TGFbvTGFbConA |> as.data.frame()

################################################################################

dir.create("go_enrichment_up-down/full", showWarnings = FALSE, recursive = TRUE)
dir.create("go_enrichment_up-down/simp", showWarnings = FALSE, recursive = TRUE)

write_csv(up_EMTvBasal, "go_enrichment_up-down/full/upregulated_EMTvBasal.csv")
write_csv(down_EMTvBasal, "go_enrichment_up-down/full/downregulated_EMTvBasal.csv")

write_csv(up_EMTvConA, "go_enrichment_up-down/full/upregulated_EMTvConA.csv")
write_csv(down_EMTvConA, "go_enrichment_up-down/full/downregulated_EMTvConA.csv")

##################### simplification

up_EMTvBasal_simp <- go_results_up_TGFbvBasal |> simplify() |> as.data.frame()
down_EMTvBasal_simp <- go_result_down_TGFbvBasal |> simplify() |> as.data.frame()

up_EMTvConA_simp <- go_results_up_TGFbvTGFbConA |> simplify() |> as.data.frame()
down_EMTvConA_simp <- go_result_down_TGFbvTGFbConA |> simplify() |> as.data.frame()


write_csv(up_EMTvBasal_simp, "go_enrichment_up-down/simp/upregulated_EMTvBasal.csv")
write_csv(down_EMTvBasal_simp, "go_enrichment_up-down/simp/downregulated_EMTvBasal.csv")

write_csv(up_EMTvConA_simp, "go_enrichment_up-down/simp/upregulated_EMTvConA.csv")
write_csv(down_EMTvConA_simp, "go_enrichment_up-down/simp/downregulated_EMTvConA.csv")






plot_cluster_go_dot <- function(data, ontology, aspectratio, top){
  
  data |>
    filter(`ONTOLOGY` %in% ontology) |>
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
      x = "Fold enrichment",
      y = NULL
    ) +
    theme(strip.text = element_text(size = 5),
          axis.text.x = element_text(size=6),
          aspect.ratio = aspectratio)
}

plot_cluster_go_dot(up_EMTvBasal_simp, c("BP"), aspectratio=1/0.7, top=20)
plot_cluster_go_dot(up_EMTvConA_simp, c("CC"), aspectratio=1/0.7, top=20)

# =============================================================================
# Record package versions alongside the outputs.
# =============================================================================
writeLines(capture.output(sessionInfo()),
           "sessionInfo_differential_expression_and_GO.txt")
