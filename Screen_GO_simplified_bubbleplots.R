# =============================================================================
# GO enrichment of the genome-wide CRISPRa ER-phagy screen, and the bubble
# plots of the enriched terms.
#
# Paper   : EMT Induces FAM134A-mediated ER-phagy to Regulate Procollagen
#           Secretion. Screen data: GEO GSE315931.
#
# Produces: Fig 2A - GO Biological Process (upper) and Molecular Function
#                    (lower) enrichment of genes whose overexpression
#                    ENHANCES ER-phagy
#           Fig 2B - the same for genes whose overexpression INHIBITS it
#           Supp Table 1 - sheets GO_BP_* and GO_MF_*, and symbol_map.csv
#
# Input   : either of two files, detected from the extension -
#             Supp_Table1-CRISPRa_starved_screen.xlsx  (as published)
#             CRISPRa_starved_genetable_collapsed_namesfixed.txt
#               (ScreenProcessing's collapsed gene table, the original input)
#           Both carry the same numbers and give byte-identical output. See
#           "INPUT" below.
#
# Writes  : <output_dir>/
#             GO_BP_simplified_bubbleplot.{pdf,png}      <- Fig 2A/2B upper
#             GO_{Enhanced,Inhibited}_{BP,MF}_simplified_results.csv
#             supp_tables/GO_{Enhanced,Inhibited}_{BP,MF}.csv
#             supp_tables/symbol_map.csv
#             sessionInfo_screen_GO.txt
#
# Run     : top to bottom. Set input_file and output_dir immediately below;
#           nothing else needs editing. Independent of the two proteomics
#           scripts in this directory.
#
# Note    : three parameters below change which terms appear and are flagged
#           PARAMETER where they are set - maxGSSize, the phenotype cutoff,
#           and n_top. The values here are the ones used for the published
#           figure.
# =============================================================================


# ---- Configuration ----------------------------------------------------------
# Edit these two lines only. Paths may be absolute, or relative to the working
# directory; the script does not call setwd().

# input_file may be EITHER the published Supplementary Table 1 (.xlsx) or
# ScreenProcessing's gene table (.txt); the script works out which from the
# extension. Table S1 is the default because it is self-contained: everything
# the script needs is published in it.

input_file <- "Supp_Table1-CRISPRa_starved_screen.xlsx"
# input_file <- "CRISPRa_starved_result/CRISPRa_starved_genetable_collapsed_namesfixed.txt"

output_dir <- "screen_GO_output"


# ---- Requirements -----------------------------------------------------------
# Required: clusterProfiler, org.Hs.eg.db, AnnotationDbi, ggplot2, patchwork,
#           stringr
# Optional: readxl (only to read Supplementary Table 1), ragg (PNG export),
#           systemfonts (font check). Without the last two the script still
#           runs; see "GRAPHICS DEVICES" below.

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# AnnotationDbi::select/keys/mapIds are called with their namespace throughout,
# because dplyr::select masks select() if the tidyverse is also loaded.

dir.create(file.path(output_dir, "supp_tables"), showWarnings = FALSE,
           recursive = TRUE)


################################################################################
# INPUT
################################################################################
# The gene table is ScreenProcessing output. Each row is one gene; columns are
# named <comparison>-<replicate>-<statistic>. The comparison used here is
# EnIn = Enhanced vs Inhibited, i.e. the sorted top 25% against the bottom 25%
# of reporter signal.
#
# "namesfixed": the hCRISPRa-v2 library contains genes whose symbols Excel
# silently converts to dates (SEPT7 -> Sep.07, MARCH1 -> Mar.01, and so on).
# Those entries were restored to their library symbols before this script ran.
# If you regenerate the table from ScreenProcessing, do not open it in Excel.
#
# The library pairs every real gene with a matched negative-control
# "pseudogene" built from non-targeting sgRNAs. These 18,915 pseudo rows are
# dropped here; they are the screen's null distribution, not candidate hits.

# The two accepted inputs hold the same numbers under different column names,
# so each is read into the same two-column frame: gene, phenotype.

read_screen_table <- function(path) {

  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    # --- Supplementary Table 1, as published ---------------------------------
    # Sheet Gene_level_results holds the 18,915 library genes and
    # Negative_controls their 18,915 matched control pseudogenes, with
    # ScreenProcessing's columns renamed for publication.
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Reading Supplementary Table 1 needs the readxl package; ",
           "install it, or point input_file at the .txt gene table instead.")

    phenotype_col <- "Enh vs Inh | mean Rep1-2 | phenotype"
    genes  <- as.data.frame(readxl::read_excel(path, sheet = "Gene_level_results"))
    pseudo <- as.data.frame(readxl::read_excel(path, sheet = "Negative_controls"))
    stopifnot(phenotype_col %in% names(genes),
              phenotype_col %in% names(pseudo))

    message(sprintf("Read Supplementary Table 1: %d genes + %d control pseudogenes.",
                    nrow(genes), nrow(pseudo)))
    rbind(
      data.frame(gene = genes[["Library symbol"]],
                 phenotype = genes[[phenotype_col]], stringsAsFactors = FALSE),
      data.frame(gene = pseudo[["Negative-control pseudogene"]],
                 phenotype = pseudo[[phenotype_col]], stringsAsFactors = FALSE))

  } else {
    # --- ScreenProcessing gene table -----------------------------------------
    phenotype_col <- "EnIn-ave_Rep1_Rep2-average phenotype of strongest 3"
    raw <- read.delim(path, check.names = FALSE)
    names(raw)[1] <- "gene"        # ScreenProcessing leaves this header blank
    stopifnot(phenotype_col %in% names(raw))

    message(sprintf("Read ScreenProcessing gene table: %d rows.", nrow(raw)))
    data.frame(gene = raw$gene, phenotype = raw[[phenotype_col]],
               stringsAsFactors = FALSE)
  }
}

stopifnot(file.exists(input_file))
screen <- read_screen_table(input_file)
screen <- screen[!grepl("^pseudo", screen$gene), ]


################################################################################
# LIBRARY SYMBOL -> CURRENT HGNC SYMBOL
################################################################################
# hCRISPRa-v2 was built in 2016, so some library symbols are now aliases.
# enrichGO silently discards symbols it cannot map, which would drop those
# genes from the enrichment without warning. Each library symbol is therefore
# resolved to a current symbol first:
#
#   already current          -> kept unchanged
#   maps to exactly one current symbol via an alias -> updated
#   maps to several          -> left alone (ambiguous; no safe choice)
#   no match at all          -> left alone, and reported at the end
#
# symbol_map.csv records the decision for every gene, so the mapping is
# auditable rather than implicit.

current_symbols <- AnnotationDbi::keys(org.Hs.eg.db, "SYMBOL")
non_current <- setdiff(screen$gene, current_symbols)

aliases <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = non_current,
                        keytype = "ALIAS", columns = "SYMBOL"))
aliases <- aliases[!is.na(aliases$SYMBOL), ]

# keep only aliases that resolve to a single current symbol
alias_n <- table(aliases$ALIAS)
unique_aliases <- aliases[aliases$ALIAS %in% names(alias_n)[alias_n == 1], ]

symbol_map <- data.frame(library_symbol = screen$gene, stringsAsFactors = FALSE)
symbol_map$current_symbol <- ifelse(
  symbol_map$library_symbol %in% current_symbols,
  symbol_map$library_symbol,
  unique_aliases$SYMBOL[match(symbol_map$library_symbol, unique_aliases$ALIAS)])

symbol_map$status <- ifelse(
  symbol_map$library_symbol %in% current_symbols, "current",
  ifelse(!is.na(symbol_map$current_symbol), "updated (unique alias)",
  ifelse(symbol_map$library_symbol %in% names(alias_n),
         "ambiguous alias - not updated", "not found")))

symbol_map$entrez <- NA_character_
mapped <- !is.na(symbol_map$current_symbol)
symbol_map$entrez[mapped] <- unname(suppressMessages(
  AnnotationDbi::mapIds(org.Hs.eg.db, symbol_map$current_symbol[mapped],
                        "ENTREZID", "SYMBOL", multiVals = "first")))

write.csv(symbol_map, file.path(output_dir, "supp_tables", "symbol_map.csv"),
          row.names = FALSE)
cat("Symbol mapping:\n"); print(table(symbol_map$status))


################################################################################
# HIT LISTS
################################################################################
# PARAMETER (phenotype_cutoff). ScreenProcessing's "average phenotype of
# strongest 3" averages the three strongest-effect sgRNAs per gene, averaged
# again over the two biological replicates. Positive = overexpression enhances
# ER-phagy, negative = inhibits it. +/- 0.7 is the cutoff used for the
# published figure; relaxing it adds weaker hits and broadens the enrichment.

phenotype_cutoff <- 0.7          # PARAMETER

phenotype <- screen$phenotype

to_current <- function(g) {
  s <- symbol_map$current_symbol[match(g, symbol_map$library_symbol)]
  ifelse(is.na(s), g, s)
}

hit_lists <- list(
  Enhanced  = to_current(screen$gene[phenotype >=  phenotype_cutoff]),
  Inhibited = to_current(screen$gene[phenotype <= -phenotype_cutoff]))

cat(sprintf("Hits: Enhanced %d, Inhibited %d\n",
            length(hit_lists$Enhanced), length(hit_lists$Inhibited)))


################################################################################
# GO ENRICHMENT
################################################################################
# PARAMETER (maxGSSize). Raised from the clusterProfiler default of 500 to
# 1000 so that "MAPK cascade" (777 annotated genes) is tested at all. Raising
# it further lets very broad parent terms (>1000 genes) dominate the result.
#
# qvalueCutoff = 1 disables the separate q-value filter; significance is
# governed by pvalueCutoff = 0.05 applied to the BH-adjusted p-value.
#
# simplify() removes terms whose semantic similarity to a more significant
# term exceeds 0.7, keeping the more significant one.

go_min_size <- 10
go_max_size <- 1000          # PARAMETER - see above
go_p_cutoff <- 0.05
simplify_cutoff <- 0.7

run_go <- function(genes, ont) {
  ego <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = ont,
                  pvalueCutoff = go_p_cutoff, qvalueCutoff = 1,
                  minGSSize = go_min_size, maxGSSize = go_max_size)
  list(full   = ego,
       simple = simplify(ego, cutoff = simplify_cutoff, by = "p.adjust",
                         select_fun = min))
}


################################################################################
# TERM SELECTION FOR THE FIGURE
################################################################################
# PARAMETER (n_top). The figure shows the n_top most significant
# non-redundant terms, plus any significant MAPK/ERK term even if simplify()
# merged it into a relative - those terms are the subject of Supp Fig 2 and
# would otherwise be invisible in the main figure.
#
# EMT-related terms are only emphasised in bold, never force-included: every
# EMT term shown reached the top list on its own p-value. This matters for the
# claim in the text, so it is deliberately not a force-include.

n_top <- 15                                             # PARAMETER
mapk_pattern <- "MAPK|ERK1|JNK|p38|Ras protein signal"  # bolded AND included
emt_pattern  <- "differentiation|epithelial cell migration|epithelium migration"
                                                        # bolded only

pick_terms <- function(res) {
  simple <- as.data.frame(res$simple)
  top <- head(simple[order(simple$p.adjust), ], n_top)

  full <- as.data.frame(res$full)
  mapk <- full[grepl(mapk_pattern, full$Description) &
               full$p.adjust < go_p_cutoff, ]

  out <- rbind(top, mapk[!mapk$ID %in% top$ID, ])
  out$MAPK <- grepl(mapk_pattern, out$Description)
  out$Bold <- out$MAPK | grepl(emt_pattern, out$Description)
  out
}


################################################################################
# FIGURE GEOMETRY AND GRAPHICS DEVICES
################################################################################
# The figure is composed on a notional 20-inch canvas and scaled down to the
# journal's full page width. Every font size, point size, line width and
# margin is multiplied by k, so the proportions are fixed and only the output
# size changes. Text is floored at 5.5 pt, the journal's minimum.

fig_width_in <- 183 / 25.4      # Developmental Cell full page width, 183 mm
k <- fig_width_in / 20
pt <- function(x, min_pt = 5.5) max(x * k, min_pt)

panel_colours <- list(Enhanced  = c("#dfe6f1", "#1f2f5c"),   # pale -> navy
                      Inhibited = c("#f7ddd3", "#9e1b16"))   # pale -> red
ontology_names <- c(BP = "GO Biological Process",
                    MF = "GO Molecular Function")
stripe_fill <- "grey60"
stripe_alpha <- 0.18

# Arial Narrow is a system font. Where it is missing, fall back to the device
# default rather than failing: the layout is unchanged, only the typeface.
# The check itself must never stop the analysis: enumerating system fonts can
# fail or exhaust memory on some systems, in which case the font is left as it
# is and the graphics device substitutes silently.
label_font <- "Arial Narrow"
font_available <- tryCatch(
  requireNamespace("systemfonts", quietly = TRUE) &&
    label_font %in% systemfonts::system_fonts()$family,
  error = function(e) NA)
if (isFALSE(font_available)) {
  message("Arial Narrow not installed; using the default font.")
  label_font <- ""
} else if (is.na(font_available)) {
  message("Could not enumerate system fonts; leaving the font to the device.")
}

# The published PDF was written with the macOS Quartz device, which embeds
# system fonts correctly. On other platforms cairo_pdf is used instead; glyph
# spacing may differ very slightly from the published figure.
pdf_device <- if (capabilities("aqua")) {
  function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...)
} else {
  grDevices::cairo_pdf
}
png_device <- if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png
} else {
  message("ragg not installed; PNG will use the default device.")
  NULL
}

# ggplot2 3.5.0 moved inside-panel legends to legend.position.inside.
legend_inside <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.03),
        legend.justification = c(1, 0))
} else {
  theme(legend.position = c(0.98, 0.03),
        legend.justification = c(1, 0))
}


################################################################################
# PANEL CONSTRUCTION
################################################################################
# Each panel is a label column plus a bubble plot, drawn as two plots that
# share the same y scale and the same alternating grey row bands, so that each
# GO term reads across to its own dot. Returns both halves plus the row count,
# which the caller needs to set the figure height.

bubble <- function(df, set_name, ont, tag) {
  df$ratio <- sapply(strsplit(df$GeneRatio, "/"),
                     function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df$lab <- str_wrap(df$Description, 40)
  df <- df[order(df$ratio, -df$p.adjust), ]
  df$row <- seq_len(nrow(df))
  df$face <- ifelse(df$Bold, "bold", "plain")

  n <- nrow(df)
  stripes <- data.frame(row = seq(2, n, by = 2))
  y_scale <- scale_y_continuous(limits = c(0.4, n + 0.6), expand = c(0, 0))
  band <- geom_rect(data = stripes,
                    aes(xmin = -Inf, xmax = Inf,
                        ymin = row - 0.5, ymax = row + 0.5),
                    inherit.aes = FALSE,
                    fill = stripe_fill, alpha = stripe_alpha)

  labels <- ggplot(df, aes(x = 1, y = row)) +
    band +
    geom_text(aes(label = lab, fontface = face), hjust = 1, nudge_x = -0.015,
              size = pt(20) / .pt, lineheight = 0.75,
              family = label_font, colour = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    y_scale +
    labs(tag = tag) +
    theme_void() +
    theme(plot.tag = element_text(face = "bold", size = pt(26)),
          plot.tag.position = c(0, 1),
          plot.margin = margin(30 * k, 0, 5 * k, 5 * k))

  dots <- ggplot(df, aes(ratio, row, size = Count,
                         colour = -log10(p.adjust))) +
    band +
    geom_point() +
    y_scale +
    scale_colour_gradient(low = panel_colours[[set_name]][1],
                          high = panel_colours[[set_name]][2],
                          name = expression(-log[10]*"(adj. "*italic(p)*")")) +
    scale_size_continuous(range = c(2, 7) * k, name = "Count",
                          breaks = function(x) unique(round(pretty(x)))) +
    guides(colour = guide_colourbar(order = 1),
           size = guide_legend(order = 2)) +
    labs(x = "Gene ratio", y = NULL,
         title = paste(set_name, "-", ontology_names[[ont]])) +
    theme_bw(base_size = 14 * k) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(30 * k, 10 * k, 5 * k, 0),
          axis.text.x = element_text(size = pt(14)),
          axis.title.x = element_text(size = pt(17)),
          legend.title = element_text(size = pt(15)),
          legend.text = element_text(size = pt(13)),
          legend.key.size = unit(pt(15), "pt"),
          legend.margin = margin(4, 4, 4, 4),
          plot.title = element_text(face = "bold", hjust = 0.5,
                                    size = pt(18)),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black",
                                               linewidth = 0.8 * k),
          legend.box = "vertical", legend.box.just = "left",
          panel.border = element_rect(colour = "black",
                                      linewidth = 0.8 * k)) +
    legend_inside

  list(labels = labels, dots = dots, n = n)
}


################################################################################
# RUN
################################################################################
# Both ontologies are enriched and both sets of tables are written, because
# Fig 2A and 2B each show BP above and MF below. Only the BP composite is
# assembled here as a single file; the MF panels were laid out separately.

for (ont in c("BP", "MF")) {
  panels <- list()

  for (set_name in names(hit_lists)) {
    res <- run_go(hit_lists[[set_name]], ont)
    figure_terms <- pick_terms(res)

    # one table per set and ontology, flagging which terms survived simplify()
    # and which reached the figure, so the figure can be traced to the table
    full <- as.data.frame(res$full)
    full$NonRedundant <- full$ID %in% as.data.frame(res$simple)$ID
    full$InFigure <- full$ID %in% figure_terms$ID

    write.csv(full,
              file.path(output_dir, "supp_tables",
                        sprintf("GO_%s_%s.csv", set_name, ont)),
              row.names = FALSE)
    write.csv(as.data.frame(res$simple),
              file.path(output_dir,
                        sprintf("GO_%s_%s_simplified_results.csv",
                                set_name, ont)),
              row.names = FALSE)

    cat(sprintf("%-9s %s: %5d terms, %4d non-redundant, %2d plotted (%d MAPK/ERK)\n",
                set_name, ont, nrow(full), sum(full$NonRedundant),
                nrow(figure_terms), sum(figure_terms$MAPK)))

    panels[[set_name]] <- bubble(figure_terms, set_name, ont,
                                 if (set_name == "Enhanced") "A" else "B")
  }

  if (ont != "BP") next   # MF tables are written above; only BP is composited

  composite <- wrap_plots(panels$Enhanced$labels, panels$Enhanced$dots,
                          panels$Inhibited$labels, panels$Inhibited$dots,
                          nrow = 1, widths = c(1.12, 1, 1.12, 1))

  # height scales with the longest panel so rows keep a constant pitch
  fig_height_in <- (1.6 + 0.54 * max(sapply(panels, `[[`, "n"))) * k

  suppressWarnings({
    ggsave(file.path(output_dir,
                     sprintf("GO_%s_simplified_bubbleplot.pdf", ont)),
           composite, width = fig_width_in, height = fig_height_in,
           device = pdf_device)
    ggsave(file.path(output_dir,
                     sprintf("GO_%s_simplified_bubbleplot.png", ont)),
           composite, width = fig_width_in, height = fig_height_in,
           dpi = 450, device = png_device)
  })
}


################################################################################
# UNMAPPED HITS
################################################################################
# Hit genes with no current HGNC symbol are silently dropped by enrichGO.
# They are listed here so that the exclusion is on the record.

for (set_name in names(hit_lists)) {
  unmapped <- setdiff(hit_lists[[set_name]], current_symbols)
  cat(sprintf("%-9s hits: %d, unmapped: %d%s\n",
              set_name, length(hit_lists[[set_name]]), length(unmapped),
              if (length(unmapped)) paste0(" (", paste(unmapped, collapse = ", "), ")") else ""))
}

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo_screen_GO.txt"))
