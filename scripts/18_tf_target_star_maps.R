#!/usr/bin/env Rscript

# Generate TF-centred star maps from the exact Top-2,000 Enrichr overlap table.
#
# Each spoke connects one selected transcription-factor target set to one gene
# found in the corresponding Enrichr overlap. Spokes do NOT establish TF
# activity, direct binding, regulatory direction, or causality.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

INPUT_FILE <- "results/tf_target_spider_maps/top2000_tf_target_evidence_table.csv"
OUTPUT_DIR <- "results/tf_target_spider_maps/R_star_maps"
COMBINED_PDF <- file.path(OUTPUT_DIR, "top2000_TF_target_star_maps_R.pdf")

FIGURE_WIDTH <- 13.33
FIGURE_HEIGHT <- 7.50
PNG_DPI <- 300

# Canonical regulator symbols used for display and deduplication.
# RNF2 is the official gene symbol; RING1B is the familiar protein/regulator
# label already used in the professor-facing figures.
REGULATOR_ALIASES <- c("RNF2" = "RING1B")

canonical_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  aliased <- unname(REGULATOR_ALIASES[x])
  replace <- !is.na(aliased)
  x[replace] <- aliased[replace]
  x
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

direction_palette <- function(direction_name) {
  if (identical(direction_name, "downregulated")) {
    list(hub = "#0072B2", target = "#B9DDF0", edge = "#5D8FA8")
  } else if (identical(direction_name, "upregulated")) {
    list(hub = "#D55E00", target = "#F4C9A8", edge = "#B87A4D")
  } else {
    stop("Unsupported direction: ", direction_name, call. = FALSE)
  }
}

required_packages <- c(
  "readr", "dplyr", "tibble", "ggplot2", "igraph", "ggraph", "ggrepel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input evidence table not found: ", INPUT_FILE,
    "\nRun scripts/04_threshold_deg_tf_sensitivity.py and ",
    "scripts/17_tf_target_deg_spider_maps.py first.",
    call. = FALSE
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

evidence <- read_csv(INPUT_FILE, show_col_types = FALSE)

required_columns <- c(
  "regulator", "direction", "overlap_gene", "selected_for_figure"
)
missing_columns <- setdiff(required_columns, names(evidence))
if (length(missing_columns) > 0L) {
  stop(
    "Evidence table is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Accept logical, numeric, or textual TRUE values without silently treating
# arbitrary non-empty strings as selected.
is_selected <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x == 1)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

evidence <- evidence %>%
  mutate(
    regulator_display = canonical_symbol(regulator),
    overlap_gene_display = canonical_symbol(overlap_gene),
    direction = tolower(trimws(as.character(direction))),
    selected_flag = is_selected(selected_for_figure)
  )

selected_pairs <- evidence %>%
  filter(selected_flag) %>%
  distinct(direction, regulator_display) %>%
  arrange(factor(direction, levels = c("downregulated", "upregulated")),
          regulator_display)

if (nrow(selected_pairs) == 0L) {
  stop("No rows are marked selected_for_figure in the evidence table.",
       call. = FALSE)
}

plot_tf_star <- function(tf_name, direction_name) {
  tf_name <- canonical_symbol(tf_name)
  direction_name <- tolower(trimws(direction_name))
  colours <- direction_palette(direction_name)

  selected <- evidence %>%
    filter(
      regulator_display == tf_name,
      direction == direction_name,
      selected_flag
    ) %>%
    filter(
      !is.na(overlap_gene_display),
      overlap_gene_display != "",
      overlap_gene_display != tf_name
    ) %>%
    distinct(overlap_gene_display, .keep_all = TRUE) %>%
    arrange(overlap_gene_display)

  if (nrow(selected) == 0L) {
    stop(
      "No peripheral overlap genes found for ", tf_name,
      " (", direction_name, ").",
      call. = FALSE
    )
  }

  edges <- selected %>%
    transmute(
      from = tf_name,
      to = overlap_gene_display,
      edge_type = "Exact Enrichr overlap"
    )

  node_names <- unique(c(tf_name, edges$to))
  nodes <- tibble(
    name = node_names,
    node_type = if_else(
      node_names == tf_name,
      "TF target set",
      "Overlap gene"
    )
  )

  graph <- graph_from_data_frame(
    d = edges,
    directed = FALSE,
    vertices = nodes
  )

  centre_index <- which(V(graph)$name == tf_name)
  if (length(centre_index) != 1L) {
    stop("Could not identify a unique central node for ", tf_name,
         call. = FALSE)
  }

  star_order <- c(centre_index, setdiff(seq_len(vcount(graph)), centre_index))
  coordinates <- layout_as_star(
    graph,
    center = centre_index,
    order = star_order
  )
  # Add radial space so the labelled peripheral nodes do not crowd the hub.
  coordinates <- coordinates * 1.22

  plot_data <- create_layout(
    graph,
    layout = "manual",
    x = coordinates[, 1],
    y = coordinates[, 2]
  )

  # Decorative spider-web guides. These polygons improve readability and
  # visual resemblance to a spider map, but deliberately carry no biological
  # meaning. Real evidence remains encoded only by the TF-to-gene spokes.
  peripheral_coordinates <- as_tibble(plot_data) %>%
    filter(node_type == "Overlap gene") %>%
    mutate(angle = atan2(y, x)) %>%
    arrange(angle)

  web_guides <- bind_rows(lapply(
    c(0.34, 0.64, 1.00),
    function(ring_fraction) {
      ring <- peripheral_coordinates %>%
        transmute(
          x = x * ring_fraction,
          y = y * ring_fraction,
          ring = factor(ring_fraction)
        )
      bind_rows(ring, slice(ring, 1L))
    }
  ))

  libraries <- character(0)
  if ("Gene_set" %in% names(selected)) {
    libraries <- sort(unique(as.character(selected$Gene_set)))
    libraries <- libraries[!is.na(libraries) & libraries != ""]
  }
  library_line <- if (length(libraries) > 0L) {
    paste0(" | Source: ", paste(libraries, collapse = " + "))
  } else {
    ""
  }

  ggraph(plot_data) +
    geom_path(
      data = web_guides,
      aes(x = x, y = y, group = ring),
      inherit.aes = FALSE,
      colour = "#D2D7DC",
      linewidth = 0.72,
      linetype = "dotted",
      alpha = 0.95,
      lineend = "round"
    ) +
    geom_edge_link(
      colour = colours$edge,
      linewidth = 1.35,
      alpha = 0.88,
      lineend = "round"
    ) +
    # Draw labels as the nodes themselves. This keeps every gene symbol inside
    # its node and more closely matches the professor's hub-and-spoke example.
    geom_node_label(
      aes(label = if_else(node_type == "Overlap gene", name, NA_character_)),
      fill = colours$target,
      fontface = "plain",
      size = 4.6,
      colour = "#333333",
      linewidth = 0.65,
      label.padding = grid::unit(0.42, "lines"),
      label.r = grid::unit(0.65, "lines"),
      repel = FALSE,
      na.rm = TRUE
    ) +
    geom_node_label(
      aes(label = if_else(node_type == "TF target set", name, NA_character_)),
      fill = colours$hub,
      colour = "white",
      fontface = "bold",
      size = 6.2,
      linewidth = 0.85,
      label.padding = grid::unit(0.68, "lines"),
      label.r = grid::unit(0.9, "lines"),
      repel = FALSE,
      na.rm = TRUE
    ) +
    labs(
      title = paste0(tf_name, "-centered target-set overlap spider map"),
      subtitle = paste0(
        "Top-2,000 ", direction_name,
        " candidate-gene input | ", nrow(selected),
        " peripheral overlap genes", library_line
      ),
      caption = paste0(
        "Solid spokes represent exact membership in the reported Enrichr ",
        "target-set overlap; pale dotted polygons are layout guides only. ",
        "The map does not establish TF activity, direct binding, regulatory ",
        "direction, gene-gene interaction, or causality."
      )
    ) +
    coord_equal(clip = "off") +
    theme_graph(base_size = 16, base_family = "sans") +
    theme(
      plot.title = element_text(
        face = "bold", size = 22, hjust = 0.5, colour = "#142B52"
      ),
      plot.subtitle = element_text(
        size = 15, hjust = 0.5, colour = "#3D4A5C",
        margin = margin(b = 12)
      ),
      plot.caption = element_text(
        size = 11.5, colour = "#555555", hjust = 0,
        margin = margin(t = 12)
      ),
      plot.margin = margin(24, 34, 24, 34),
      legend.position = "none"
    )
}

plots <- vector("list", nrow(selected_pairs))
manifest <- vector("list", nrow(selected_pairs))

for (i in seq_len(nrow(selected_pairs))) {
  direction_name <- selected_pairs$direction[[i]]
  tf_name <- selected_pairs$regulator_display[[i]]
  message("Generating: ", direction_name, " / ", tf_name)

  p <- plot_tf_star(tf_name, direction_name)
  plots[[i]] <- p

  stem <- safe_filename(paste(direction_name, tf_name, "R_star_map", sep = "_"))
  png_path <- file.path(OUTPUT_DIR, paste0(stem, ".png"))
  pdf_path <- file.path(OUTPUT_DIR, paste0(stem, ".pdf"))

  ggsave(
    png_path, p,
    width = FIGURE_WIDTH, height = FIGURE_HEIGHT,
    dpi = PNG_DPI, bg = "white"
  )
  ggsave(
    pdf_path, p,
    width = FIGURE_WIDTH, height = FIGURE_HEIGHT,
    device = cairo_pdf, bg = "white"
  )

  peripheral_count <- evidence %>%
    filter(
      regulator_display == tf_name,
      direction == direction_name,
      selected_flag,
      !is.na(overlap_gene_display),
      overlap_gene_display != "",
      overlap_gene_display != tf_name
    ) %>%
    distinct(overlap_gene_display) %>%
    nrow()

  manifest[[i]] <- tibble(
    direction = direction_name,
    regulator = tf_name,
    peripheral_overlap_genes_shown = peripheral_count,
    png_file = png_path,
    pdf_file = pdf_path
  )
}

# One page per TF, convenient for review and supplementary submission.
grDevices::cairo_pdf(
  COMBINED_PDF,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  onefile = TRUE,
  family = "sans"
)
for (p in plots) print(p)
dev.off()

manifest_table <- bind_rows(manifest)
manifest_path <- file.path(OUTPUT_DIR, "top2000_R_star_map_manifest.csv")
write_csv(manifest_table, manifest_path)

message("\nCompleted successfully.")
message("Generated ", length(plots), " TF-centred star maps.")
message("Individual figures: ", OUTPUT_DIR)
message("Combined PDF: ", COMBINED_PDF)
message("Manifest: ", manifest_path)
