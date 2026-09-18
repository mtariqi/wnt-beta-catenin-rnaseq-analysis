#!/usr/bin/env Rscript

# Evidence-based TF target-set interaction maps
#
# Layer 1: TF -> gene spokes are exact Enrichr overlap membership from the
#          Top-2,000 evidence table.
# Layer 2: gene -> gene edges are mouse STRING v12.0 functional associations
#          among those same peripheral genes. No extra proteins are added.
#
# IMPORTANT: STRING functional associations can reflect several evidence
# channels and do not by themselves establish direct physical binding,
# regulatory direction, TF activity, or causality.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(igraph)
  library(curl)
  library(scales)
})

# ----------------------------- Configuration -----------------------------

INPUT_FILE <- "results/tf_target_spider_maps/top2000_tf_target_evidence_table.csv"
OUTPUT_DIR <- "results/tf_target_spider_maps/STRING_interaction_maps"
COMBINED_PDF <- file.path(
  OUTPUT_DIR,
  "top2000_TF_target_STRING_interaction_maps.pdf"
)

# Version-specific endpoint is used for reproducibility.
STRING_API_BASE <- "https://version-12-0.string-db.org/api"
STRING_SPECIES <- 10090L                    # Mus musculus
STRING_REQUIRED_SCORE <- 400L               # medium confidence or greater
STRING_NETWORK_TYPE <- "functional"         # functional, not physical-only
CALLER_IDENTITY <- "wnt-beta-catenin-rnaseq-analysis"
API_WAIT_SECONDS <- 1.05                     # STRING requests >=1 s spacing

FIGURE_WIDTH <- 13.33
FIGURE_HEIGHT <- 7.50
PNG_DPI <- 300
LAYOUT_SEED <- 20260916L

REGULATOR_ALIASES <- c("RNF2" = "RING1B")

# Contextual examples reported by Doumpas et al., EMBO J (2019),
# doi:10.15252/embj.201798873. These flags are annotations only: the paper used
# human HEK293T cells and does not validate STRING edges or the present mouse
# RNA-seq contrast. GADD45A is the current gene-symbol rendering of the paper's
# "GADD45" example.
PAPER_REFERENCE <- "Doumpas et al., EMBO J 2019, doi:10.15252/embj.201798873"
PAPER_CANONICAL_TCF_EXAMPLES <- c("AXIN2", "LEF1", "DKK1", "NKD1", "SP5")
PAPER_TCF_INDEPENDENT_EXAMPLES <- c("HOXC4", "ZNF503", "GADD45A")

# ------------------------------- Utilities -------------------------------

canonical_regulator <- function(x) {
  x <- toupper(trimws(as.character(x)))
  alias <- unname(REGULATOR_ALIASES[x])
  replace <- !is.na(alias)
  x[replace] <- alias[replace]
  x
}

clean_gene_symbol <- function(x) {
  toupper(trimws(as.character(x)))
}

paper_context_class <- function(x) {
  x <- clean_gene_symbol(x)
  case_when(
    x %in% PAPER_CANONICAL_TCF_EXAMPLES ~
      "Paper: canonical/TCF-dependent example",
    x %in% PAPER_TCF_INDEPENDENT_EXAMPLES ~
      "Paper: TCF-independent example",
    TRUE ~ "No paper example flag"
  )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

is_selected <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x == 1)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

direction_palette <- function(direction_name) {
  if (identical(direction_name, "downregulated")) {
    list(hub = "#0072B2", target = "#B9DDF0", spoke = "#4E8BA8")
  } else if (identical(direction_name, "upregulated")) {
    list(hub = "#D55E00", target = "#F4C9A8", spoke = "#B87849")
  } else {
    stop("Unsupported direction: ", direction_name, call. = FALSE)
  }
}

urlencode_form <- function(values) {
  stopifnot(!is.null(names(values)))
  encoded <- Map(
    function(key, value) {
      paste0(curl_escape(key), "=", curl_escape(as.character(value)))
    },
    names(values),
    values
  )
  paste(unlist(encoded, use.names = FALSE), collapse = "&")
}

string_post_tsv <- function(method, parameters, allow_not_found = FALSE) {
  endpoint <- paste(STRING_API_BASE, "tsv", method, sep = "/")
  handle <- new_handle()
  handle_setheaders(
    handle,
    "Content-Type" = "application/x-www-form-urlencoded",
    "Accept" = "text/tab-separated-values"
  )
  handle_setopt(
    handle,
    postfields = urlencode_form(parameters),
    useragent = paste0(CALLER_IDENTITY, "/1.0")
  )

  response <- tryCatch(
    curl_fetch_memory(endpoint, handle = handle),
    error = function(e) {
      stop("STRING API request failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (response$status_code == 404L && allow_not_found) {
    Sys.sleep(API_WAIT_SECONDS)
    return(tibble())
  }
  if (response$status_code < 200L || response$status_code >= 300L) {
    body <- rawToChar(response$content)
    stop(
      "STRING API returned HTTP ", response$status_code,
      if (nzchar(body)) paste0(": ", substr(body, 1, 300)) else "",
      call. = FALSE
    )
  }

  body <- rawToChar(response$content)
  Sys.sleep(API_WAIT_SECONDS)
  if (!nzchar(trimws(body))) return(tibble())

  read_tsv(I(body), show_col_types = FALSE, progress = FALSE)
}

map_string_identifiers <- function(symbols, chunk_size = 80L) {
  symbols <- sort(unique(clean_gene_symbol(symbols)))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0L) return(tibble())

  chunks <- split(symbols, ceiling(seq_along(symbols) / chunk_size))
  mapped <- lapply(chunks, function(chunk) {
    message("Mapping ", length(chunk), " identifiers to STRING IDs...")
    string_post_tsv(
      "get_string_ids",
      list(
        identifiers = paste(chunk, collapse = "\r"),
        species = STRING_SPECIES,
        echo_query = 1,
        caller_identity = CALLER_IDENTITY
      ),
      allow_not_found = TRUE
    )
  })

  bind_rows(mapped) %>%
    mutate(query_symbol = clean_gene_symbol(queryItem)) %>%
    arrange(query_symbol, queryIndex) %>%
    distinct(query_symbol, .keep_all = TRUE)
}

query_string_network <- function(string_ids) {
  string_ids <- unique(as.character(string_ids))
  string_ids <- string_ids[!is.na(string_ids) & nzchar(string_ids)]
  if (length(string_ids) < 2L) return(tibble())

  string_post_tsv(
    "network",
    list(
      identifiers = paste(string_ids, collapse = "\r"),
      species = STRING_SPECIES,
      required_score = STRING_REQUIRED_SCORE,
      network_type = STRING_NETWORK_TYPE,
      add_nodes = 0,
      caller_identity = CALLER_IDENTITY
    ),
    allow_not_found = TRUE
  )
}

edge_coordinates <- function(edges, node_coordinates) {
  # Preserve the expected plotting columns even when STRING reports no
  # qualifying interactions. A zero-column tibble makes ggplot2 fail while
  # evaluating x_from/y_from in an otherwise valid empty layer.
  if (nrow(edges) == 0L) {
    return(
      edges %>%
        mutate(
          x_from = numeric(),
          y_from = numeric(),
          x_to = numeric(),
          y_to = numeric()
        )
    )
  }
  edges %>%
    left_join(
      node_coordinates %>% select(from = name, x_from = x, y_from = y),
      by = "from"
    ) %>%
    left_join(
      node_coordinates %>% select(to = name, x_to = x, y_to = y),
      by = "to"
    )
}

# ------------------------------- Input data ------------------------------

required_packages <- c(
  "readr", "dplyr", "tibble", "ggplot2", "igraph", "curl", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
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

evidence <- evidence %>%
  mutate(
    regulator_original = clean_gene_symbol(regulator),
    regulator_display = canonical_regulator(regulator),
    overlap_gene_display = clean_gene_symbol(overlap_gene),
    direction = tolower(trimws(as.character(direction))),
    selected_flag = is_selected(selected_for_figure)
  )

selected_pairs <- evidence %>%
  filter(selected_flag) %>%
  group_by(direction, regulator_display) %>%
  summarise(regulator_original = first(regulator_original), .groups = "drop") %>%
  arrange(
    factor(direction, levels = c("downregulated", "upregulated")),
    regulator_display
  )

if (nrow(selected_pairs) == 0L) {
  stop("No rows are marked selected_for_figure.", call. = FALSE)
}

# Map every peripheral symbol once; retain unresolved symbols in the output.
all_peripheral_symbols <- evidence %>%
  filter(selected_flag) %>%
  pull(overlap_gene_display) %>%
  unique()

string_mapping <- map_string_identifiers(all_peripheral_symbols)
mapping_audit <- tibble(query_symbol = sort(unique(all_peripheral_symbols))) %>%
  left_join(string_mapping, by = "query_symbol") %>%
  mutate(mapped_to_string = !is.na(stringId))

write_csv(
  mapping_audit,
  file.path(OUTPUT_DIR, "STRING_identifier_mapping_audit.csv")
)

# ---------------------------- Map construction ---------------------------

build_tf_interaction_map <- function(
    tf_name,
    regulator_original,
    direction_name,
    pair_index) {

  tf_name <- canonical_regulator(tf_name)
  regulator_original <- clean_gene_symbol(regulator_original)
  direction_name <- tolower(direction_name)
  colours <- direction_palette(direction_name)

  self_symbols <- unique(c(tf_name, regulator_original))
  if (tf_name == "RING1B") self_symbols <- unique(c(self_symbols, "RNF2"))

  selected <- evidence %>%
    filter(
      regulator_display == tf_name,
      direction == direction_name,
      selected_flag,
      !is.na(overlap_gene_display),
      nzchar(overlap_gene_display),
      !overlap_gene_display %in% self_symbols
    ) %>%
    distinct(overlap_gene_display, .keep_all = TRUE) %>%
    arrange(overlap_gene_display)

  if (nrow(selected) == 0L) {
    stop("No peripheral genes for ", tf_name, " / ", direction_name,
         call. = FALSE)
  }

  peripheral_genes <- selected$overlap_gene_display
  mapped_for_pair <- mapping_audit %>%
    filter(query_symbol %in% peripheral_genes, mapped_to_string)

  message(
    "Querying STRING: ", direction_name, " / ", tf_name,
    " (", nrow(mapped_for_pair), "/", length(peripheral_genes), " mapped)"
  )

  string_raw <- query_string_network(mapped_for_pair$stringId)

  id_to_symbol <- mapped_for_pair %>%
    select(stringId, query_symbol)

  if (nrow(string_raw) > 0L) {
    string_edges <- string_raw %>%
      left_join(
        id_to_symbol %>% rename(stringId_A = stringId, from = query_symbol),
        by = "stringId_A"
      ) %>%
      left_join(
        id_to_symbol %>% rename(stringId_B = stringId, to = query_symbol),
        by = "stringId_B"
      ) %>%
      mutate(
        from = coalesce(from, clean_gene_symbol(preferredName_A)),
        to = coalesce(to, clean_gene_symbol(preferredName_B)),
        edge_key = if_else(from < to, paste(from, to, sep = "--"),
                           paste(to, from, sep = "--"))
      ) %>%
      filter(from %in% peripheral_genes, to %in% peripheral_genes, from != to) %>%
      arrange(desc(score)) %>%
      distinct(edge_key, .keep_all = TRUE) %>%
      select(-edge_key)
  } else {
    string_edges <- tibble(
      from = character(), to = character(), score = numeric()
    )
  }

  overlap_edges <- tibble(
    from = tf_name,
    to = peripheral_genes,
    edge_type = "Exact Enrichr overlap",
    layout_weight = 1.0
  )

  graph_edges <- bind_rows(
    overlap_edges %>% select(from, to, edge_type, layout_weight),
    string_edges %>% transmute(
      from, to,
      edge_type = "STRING functional association",
      layout_weight = 0.6 + 1.8 * score
    )
  )

  node_names <- unique(c(tf_name, peripheral_genes))
  nodes <- tibble(
    name = node_names,
    node_type = if_else(name == tf_name, "TF target set", "Overlap gene"),
    paper_context = if_else(
      name == tf_name,
      "TF hub",
      paper_context_class(name)
    )
  )

  graph <- graph_from_data_frame(graph_edges, directed = FALSE, vertices = nodes)
  set.seed(LAYOUT_SEED + pair_index)
  coordinates <- layout_with_fr(
    graph,
    weights = E(graph)$layout_weight,
    niter = 1500,
    grid = "nogrid"
  )

  tf_index <- which(V(graph)$name == tf_name)
  coordinates[, 1] <- coordinates[, 1] - coordinates[tf_index, 1]
  coordinates[, 2] <- coordinates[, 2] - coordinates[tf_index, 2]
  max_radius <- max(sqrt(rowSums(coordinates^2)))
  if (is.finite(max_radius) && max_radius > 0) coordinates <- coordinates / max_radius

  node_coordinates <- nodes %>%
    mutate(
      x = coordinates[match(name, V(graph)$name), 1],
      y = coordinates[match(name, V(graph)$name), 2]
    )

  overlap_coordinates <- edge_coordinates(overlap_edges, node_coordinates)
  string_coordinates <- edge_coordinates(string_edges, node_coordinates)

  subtitle <- paste0(
    "Top-2,000 ", direction_name, " candidate-gene input | ",
    length(peripheral_genes), " overlap genes | ",
    nrow(string_edges), " STRING edges (score >= ",
    STRING_REQUIRED_SCORE / 1000, ")"
  )

  p <- ggplot() +
    geom_segment(
      data = string_coordinates,
      aes(
        x = x_from, y = y_from, xend = x_to, yend = y_to,
        linewidth = score
      ),
      colour = "#8B8F94",
      alpha = 0.72,
      lineend = "round"
    ) +
    geom_segment(
      data = overlap_coordinates,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      colour = colours$spoke,
      linewidth = 1.12,
      alpha = 0.90,
      lineend = "round"
    ) +
    geom_label(
      data = node_coordinates %>% filter(node_type == "Overlap gene"),
      aes(x = x, y = y, label = name, colour = paper_context),
      fill = colours$target,
      family = "sans",
      fontface = "plain",
      size = 4.3,
      linewidth = 0.60,
      label.padding = grid::unit(0.35, "lines"),
      label.r = grid::unit(0.55, "lines")
    ) +
    geom_label(
      data = node_coordinates %>% filter(node_type == "TF target set"),
      aes(x = x, y = y, label = name),
      fill = colours$hub,
      colour = "white",
      family = "sans",
      fontface = "bold",
      size = 6.2,
      linewidth = 0.80,
      label.padding = grid::unit(0.62, "lines"),
      label.r = grid::unit(0.8, "lines")
    ) +
    scale_linewidth_continuous(
      name = "STRING combined score",
      range = c(0.55, 2.6),
      limits = c(STRING_REQUIRED_SCORE / 1000, 1),
      breaks = c(0.4, 0.7, 0.9),
      labels = number_format(accuracy = 0.1)
    ) +
    scale_colour_manual(
      name = "Doumpas et al. context",
      values = c(
        "No paper example flag" = "#30343A",
        "Paper: canonical/TCF-dependent example" = "#1B7837",
        "Paper: TCF-independent example" = "#762A83"
      ),
      breaks = c(
        "Paper: canonical/TCF-dependent example",
        "Paper: TCF-independent example"
      ),
      drop = TRUE
    ) +
    coord_equal(clip = "off", xlim = c(-1.25, 1.25), ylim = c(-1.18, 1.18)) +
    labs(
      title = paste0(tf_name, " target-overlap and STRING interaction map"),
      subtitle = subtitle,
      caption = paste0(
        "Colored spokes: exact Enrichr target-set overlap. Grey edges: ",
        "STRING v12.0 mouse functional associations among displayed genes; ",
        "width encodes combined score. Associations do not establish direct ",
        "binding, regulatory direction, TF activity, or causality. Green or ",
        "purple labels mark contextual examples from Doumpas et al. (human ",
        "HEK293T); they are not validation in this dataset."
      )
    ) +
    theme_void(base_family = "sans", base_size = 15) +
    theme(
      plot.title = element_text(
        face = "bold", size = 21, hjust = 0.5, colour = "#142B52"
      ),
      plot.subtitle = element_text(
        size = 14, hjust = 0.5, colour = "#3D4A5C",
        margin = margin(b = 12)
      ),
      plot.caption = element_text(
        size = 10.5, colour = "#555555", hjust = 0,
        margin = margin(t = 12)
      ),
      plot.margin = margin(24, 36, 26, 36),
      legend.position = if (nrow(string_edges) > 0L) "right" else "none",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )

  list(
    plot = p,
    string_edges = string_edges,
    selected = selected,
    node_coordinates = node_coordinates,
    mapped_count = nrow(mapped_for_pair),
    peripheral_count = length(peripheral_genes)
  )
}

# --------------------------- Generate all maps ----------------------------

plots <- vector("list", nrow(selected_pairs))
manifest <- vector("list", nrow(selected_pairs))
all_string_edges <- vector("list", nrow(selected_pairs))

for (i in seq_len(nrow(selected_pairs))) {
  pair <- selected_pairs[i, ]
  result <- build_tf_interaction_map(
    pair$regulator_display,
    pair$regulator_original,
    pair$direction,
    i
  )

  tf_name <- pair$regulator_display
  direction_name <- pair$direction
  p <- result$plot
  plots[[i]] <- p

  stem <- safe_filename(paste(
    direction_name, tf_name, "STRING_interaction_map", sep = "_"
  ))
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

  if (nrow(result$string_edges) > 0L) {
    all_string_edges[[i]] <- result$string_edges %>%
      mutate(
        regulator = tf_name,
        direction = direction_name,
        .before = 1
      )
  } else {
    all_string_edges[[i]] <- tibble()
  }

  manifest[[i]] <- tibble(
    direction = direction_name,
    regulator = tf_name,
    peripheral_overlap_genes = result$peripheral_count,
    genes_mapped_to_STRING = result$mapped_count,
    plotted_STRING_edges = nrow(result$string_edges),
    required_score = STRING_REQUIRED_SCORE / 1000,
    network_type = STRING_NETWORK_TYPE,
    png_file = png_path,
    pdf_file = pdf_path
  )
}

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
string_edge_audit <- bind_rows(all_string_edges)
if (ncol(string_edge_audit) == 0L) {
  string_edge_audit <- tibble(
    regulator = character(), direction = character(),
    from = character(), to = character(), score = numeric()
  )
}

paper_context_audit <- bind_rows(lapply(seq_along(plots), function(i) {
  pair <- selected_pairs[i, ]
  evidence %>%
    filter(
      selected_flag,
      regulator_display == pair$regulator_display,
      direction == pair$direction
    ) %>%
    transmute(
      direction,
      regulator = regulator_display,
      gene = overlap_gene_display,
      paper_context = paper_context_class(overlap_gene_display),
      paper_reference = PAPER_REFERENCE
    ) %>%
    distinct()
}))

write_csv(
  manifest_table,
  file.path(OUTPUT_DIR, "top2000_STRING_interaction_map_manifest.csv")
)
write_csv(
  string_edge_audit,
  file.path(OUTPUT_DIR, "top2000_STRING_interaction_edge_audit.csv")
)
write_csv(
  paper_context_audit,
  file.path(OUTPUT_DIR, "top2000_Doumpas2019_context_annotation.csv")
)
write_csv(
  tibble(
    parameter = c(
      "STRING_API_BASE", "species", "required_score", "network_type",
      "add_nodes", "caller_identity", "input_file", "paper_reference",
      "paper_annotation_scope"
    ),
    value = c(
      STRING_API_BASE, STRING_SPECIES, STRING_REQUIRED_SCORE,
      STRING_NETWORK_TYPE, 0, CALLER_IDENTITY, INPUT_FILE, PAPER_REFERENCE,
      paste(
        "Context-only flags for selected examples; human HEK293T paper,",
        "not validation of mouse RNA-seq or STRING edges"
      )
    )
  ),
  file.path(OUTPUT_DIR, "STRING_query_parameters.csv")
)

message("\nCompleted successfully.")
message("Generated ", length(plots), " evidence-based interaction maps.")
message("Output directory: ", OUTPUT_DIR)
message("Combined PDF: ", COMBINED_PDF)
message(
  "Audited STRING edges: ",
  file.path(OUTPUT_DIR, "top2000_STRING_interaction_edge_audit.csv")
)
