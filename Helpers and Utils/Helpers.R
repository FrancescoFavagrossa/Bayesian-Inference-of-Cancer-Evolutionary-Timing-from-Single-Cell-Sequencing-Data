default_mutation_key_cols <- function() {
  c("chr", "from", "ref", "alt")
}

get_single_value <- function(x, label) {
  if (length(x) != 1 || is.na(x)) {
    stop("Expected one non-missing value for ", label, ".", call. = FALSE)
  }

  x[[1]]
}

get_mrca_time <- function(phylo_forest, ids, label) {
  if (length(ids) == 0) {
    stop("No cell ids available for MRCA calculation: ", label, call. = FALSE)
  }

  coalescent <- phylo_forest$get_coalescent_cells(ids)
  get_single_value(coalescent[["birth_time"]], label)
}

subtract_interval <- function(regions, chr, begin, end) {
  chr <- as.character(chr)
  subtract_begin <- as.integer(begin)
  subtract_end <- as.integer(end)

  if (subtract_begin > subtract_end) {
    stop("Invalid interval: begin is greater than end.", call. = FALSE)
  }

  regions <- tibble::as_tibble(regions)
  require_columns(regions, c("chr", "begin", "end"), "regions")

  regions[["chr"]] <- as.character(regions[["chr"]])
  regions[["begin"]] <- as.integer(regions[["begin"]])
  regions[["end"]] <- as.integer(regions[["end"]])

  out <- vector("list", max(1L, nrow(regions) * 2L))
  out_idx <- 1L

  for (i in seq_len(nrow(regions))) {
    r_chr <- regions$chr[[i]]
    r_begin <- regions$begin[[i]]
    r_end <- regions$end[[i]]

    if (r_chr != chr || subtract_end < r_begin || subtract_begin > r_end) {
      out[[out_idx]] <- tibble::tibble(chr = r_chr, begin = r_begin, end = r_end)
      out_idx <- out_idx + 1L
      next
    }

    if (subtract_begin > r_begin) {
      out[[out_idx]] <- tibble::tibble(chr = r_chr, begin = r_begin, end = subtract_begin - 1L)
      out_idx <- out_idx + 1L
    }

    if (subtract_end < r_end) {
      out[[out_idx]] <- tibble::tibble(chr = r_chr, begin = subtract_end + 1L, end = r_end)
      out_idx <- out_idx + 1L
    }
  }

  out <- out[seq_len(out_idx - 1L)]

  if (length(out) == 0) {
    return(tibble::tibble(chr = character(), begin = integer(), end = integer()))
  }

  dplyr::bind_rows(out)
}

get_node_descendants <- function(node, parents) {
  if (length(node) != 1 || is.na(node)) {
    stop("node must be a single non-missing index.", call. = FALSE)
  }

  node <- as.integer(node)
  descendants <- node
  children <- which(!is.na(parents) & parents == node)

  for (child in children) {
    descendants <- c(descendants, get_node_descendants(child, parents))
  }

  descendants
}

make_common_normal_regions <- function(
    cna_regions,
    sampled_clone_order,
    genome_regions = NULL,
    phylo_forest = NULL,
    cell_label_map = NULL,
    assign_total_normal_length = TRUE,
    envir = parent.frame()
) {
  cna_regions <- tibble::as_tibble(cna_regions)
  require_columns(cna_regions, c("target_label", "chr", "begin", "end"), "cna_regions")

  if (is.null(genome_regions)) {
    if (is.null(phylo_forest)) {
      phylo_forest <- lookup_object("phylo_forest", envir)
    }

    if (is.null(cell_label_map)) {
      cell_label_map <- lookup_object("cell_label_map", envir)
    }

    if (is.null(phylo_forest) || is.null(cell_label_map)) {
      stop(
        "Either pass genome_regions or call make_common_normal_regions() after ",
        "phylo_forest and cell_label_map exist.",
        call. = FALSE
      )
    }

    common_regions <- infer_common_normal_regions_from_simulation(
      phylo_forest = phylo_forest,
      cell_label_map = cell_label_map,
      sampled_clone_order = sampled_clone_order
    )
  } else {
    genome_regions <- tibble::as_tibble(genome_regions)
    require_columns(genome_regions, c("chr", "begin", "end"), "genome_regions")

    common_regions <- purrr::map_dfr(sampled_clone_order, function(label) {
      regions <- genome_regions

      clone_cnas <- cna_regions |>
        dplyr::filter(.data$target_label == label)

      if (nrow(clone_cnas) > 0) {
        for (i in seq_len(nrow(clone_cnas))) {
          regions <- subtract_interval(
            regions = regions,
            chr = clone_cnas$chr[[i]],
            begin = clone_cnas$begin[[i]],
            end = clone_cnas$end[[i]]
          )
        }
      }

      regions |>
        dplyr::mutate(
          cell_label = label,
          chr = as.character(.data$chr),
          begin = as.integer(.data$begin),
          end = as.integer(.data$end),
          length = pmax(0L, .data$end - .data$begin + 1L),
          .before = 1
        )
    })
  }

  total_normal_length <- common_regions |>
    dplyr::group_by(.data$cell_label) |>
    dplyr::summarise(total_11_length = sum(.data$length), .groups = "drop") |>
    dplyr::right_join(
      tibble::tibble(cell_label = sampled_clone_order),
      by = "cell_label"
    ) |>
    dplyr::mutate(total_11_length = tidyr::replace_na(.data$total_11_length, 0L)) |>
    dplyr::select("cell_label", "total_11_length")

  if (isTRUE(assign_total_normal_length)) {
    assign("total_normal_length", total_normal_length, envir = envir)
  }

  common_regions
}

infer_common_normal_regions_from_simulation <- function(
    phylo_forest,
    cell_label_map,
    sampled_clone_order,
    normal_major = 1L,
    normal_minor = 1L
) {
  frag <- phylo_forest$get_cell_allelic_fragmentation()
  frag <- data.table::as.data.table(frag)
  label_map <- data.table::as.data.table(cell_label_map)

  require_columns(
    frag,
    c("cell_id", "chr", "begin", "end", "major", "minor"),
    "phylo_forest$get_cell_allelic_fragmentation()"
  )
  require_columns(label_map, c("cell_id", "cell_label"), "cell_label_map")

  label_map <- unique(label_map[, c("cell_id", "cell_label"), with = FALSE])
  frag <- merge(frag, label_map, by = "cell_id", all = FALSE)
  frag <- frag[frag[["cell_label"]] %in% sampled_clone_order]
  frag <- frag[frag[["major"]] == normal_major & frag[["minor"]] == normal_minor]

  if (nrow(frag) == 0) {
    stop("No simulated 1:1 regions found for sampled clone labels.", call. = FALSE)
  }

  frag[, chr := as.character(chr)]
  frag[, begin := as.integer(begin)]
  frag[, end := as.integer(end)]

  clone_n <- unique(frag[, .(cell_id, cell_label)])
  clone_n <- clone_n[, .(n_cells = data.table::uniqueN(cell_id)), by = cell_label]

  normal_regions <- unique(frag[, .(cell_label, cell_id, chr, begin, end)])
  normal_regions <- normal_regions[
    ,
    .(n_cells = data.table::uniqueN(cell_id)),
    by = .(cell_label, chr, begin, end)
  ]
  normal_regions <- merge(normal_regions, clone_n, by = "cell_label", suffixes = c("_region", "_clone"))
  normal_regions <- normal_regions[normal_regions[["n_cells_region"]] == normal_regions[["n_cells_clone"]]]
  normal_regions <- normal_regions[, .(cell_label, chr, begin, end)]
  normal_regions[, length := pmax(0L, end - begin + 1L)]

  out <- tibble::as_tibble(normal_regions) |>
    dplyr::arrange(.data$cell_label, .data$chr, .data$begin, .data$end)

  missing_labels <- setdiff(sampled_clone_order, unique(out$cell_label))
  if (length(missing_labels) > 0) {
    stop(
      "No common simulated 1:1 regions found for: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  out
}

get_diploid_regions <- function(common_normal_regions, region_labels) {
  if (length(region_labels) == 0) {
    stop("region_labels must contain at least one clone label.", call. = FALSE)
  }

  regions <- tibble::as_tibble(common_normal_regions)
  require_columns(regions, c("cell_label", "chr", "begin", "end"), "common_normal_regions")

  regions <- regions |>
    dplyr::filter(.data$cell_label %in% region_labels) |>
    dplyr::mutate(
      chr = as.character(.data$chr),
      begin = as.integer(.data$begin),
      end = as.integer(.data$end)
    )

  if (length(region_labels) == 1) {
    return(
      regions |>
        dplyr::distinct(.data$chr, .data$begin, .data$end)
    )
  }

  missing_labels <- setdiff(region_labels, unique(regions$cell_label))
  if (length(missing_labels) > 0) {
    stop(
      "No normal regions available for: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  common <- regions |>
    dplyr::filter(.data$cell_label == region_labels[[1]]) |>
    dplyr::distinct(.data$chr, .data$begin, .data$end)

  for (label in region_labels[-1]) {
    label_regions <- regions |>
      dplyr::filter(.data$cell_label == label) |>
      dplyr::distinct(.data$chr, .data$begin, .data$end)

    common <- dplyr::inner_join(
      common,
      label_regions,
      by = "chr",
      suffix = c(".common", ".label")
    ) |>
      dplyr::transmute(
        chr = .data$chr,
        begin = pmax(.data$begin.common, .data$begin.label),
        end = pmin(.data$end.common, .data$end.label)
      ) |>
      dplyr::filter(.data$begin <= .data$end) |>
      dplyr::distinct(.data$chr, .data$begin, .data$end)

    if (nrow(common) == 0) {
      break
    }
  }

  common
}

build_mutation_clean <- function(
    mutation_row,
    cell_label_map,
    sampled_clone_order,
    mutation_key_cols = default_mutation_key_cols()
) {
  dt <- data.table::as.data.table(mutation_row)
  map <- data.table::as.data.table(cell_label_map)

  require_columns(dt, c("cell_id", "type", "class", mutation_key_cols, "allele"), "mutation_row")
  require_columns(map, c("cell_id", "cell_label"), "cell_label_map")

  map <- unique(map[, c("cell_id", "cell_label"), with = FALSE])
  dt <- dt[dt[["type"]] == "SNV" & (is.na(dt[["class"]]) | dt[["class"]] != "pre-neoplastic")]
  dt <- dt[dt[["cell_id"]] %in% map[["cell_id"]]]
  dt <- merge(dt, map, by = "cell_id", all = FALSE)
  dt <- dt[dt[["cell_label"]] %in% sampled_clone_order]

  if ("chr" %in% mutation_key_cols) {
    dt[, chr := as.character(chr)]
  }

  cols <- c("cell_label", "cell_id", mutation_key_cols, "allele")
  tibble::as_tibble(unique(dt[, cols, with = FALSE]))
}

build_mutation_presence <- function(
    mut_clean,
    mutation_key_cols = default_mutation_key_cols()
) {
  dt <- data.table::as.data.table(mut_clean)
  require_columns(dt, c("cell_label", "cell_id", mutation_key_cols), "mut_clean")

  if ("chr" %in% mutation_key_cols) {
    dt[, chr := as.character(chr)]
  }

  cols <- c("cell_label", "cell_id", mutation_key_cols)
  dt <- unique(dt[, cols, with = FALSE])
  out <- dt[, .(n_cells = data.table::uniqueN(cell_id)), by = c(mutation_key_cols, "cell_label")]

  tibble::as_tibble(out)
}

summarise_mutation_presence <- function(
    mutation_presence,
    required_labels,
    compatible_labels = character(0),
    excluded_labels = character(0),
    mutation_key_cols = default_mutation_key_cols()
) {
  dt <- data.table::as.data.table(mutation_presence)
  require_columns(dt, c(mutation_key_cols, "cell_label", "n_cells"), "mutation_presence")

  if ("chr" %in% mutation_key_cols) {
    dt[, chr := as.character(chr)]
  }

  if (nrow(dt) == 0) {
    out <- unique(dt[, mutation_key_cols, with = FALSE])
    out[, `:=`(n_req = integer(), n_compatible = integer(), n_exc = integer())]
    return(tibble::as_tibble(out))
  }

  out <- dt[, .(
    n_req = sum(n_cells[cell_label %in% required_labels]),
    n_compatible = sum(n_cells[cell_label %in% compatible_labels]),
    n_exc = sum(n_cells[cell_label %in% excluded_labels])
  ), by = mutation_key_cols]

  tibble::as_tibble(out)
}

#' Count clonal SNVs on a clone/node branch (sibling-independent)
#'
#' Returns the SNV loci that are present in EVERY sampled cell of `clone_label`.
#' `clone_label` may be a single residual clone label or the full descendant
#' label set for a node. The latter is needed for MRCA-node branch counts: a
#' mutation on the branch before a child MRCA must be clonal across all sampled
#' descendants of that child, not only across the residual cells with that label.
#'
#' The count subtracts loci also fixed/clonal in `parent_label`, and restricts
#' loci to `region_label`'s own 1:1 (diploid) regions. This keeps later
#' descendant/sibling CNAs from shrinking the ancestor branch's normal length.
#'
#' Unlike the previous definition it does NOT exclude mutations seen in sibling
#' clones, and it subtracts the parent's *fixed* set rather than requiring strict
#' absence from the parent. Mutations on the shared pre-divergence stem (present
#' in a brother, or subclonal in the parent) are therefore retained, which is the
#' fix for the brother-clone under-count.
count_m <- function(
    mutation_presence,
    common_normal_regions,
    clone_label,
    parent_label = NULL,
    region_label = clone_label[[1]],
    n_cells_per_clone,
    mutation_key_cols = default_mutation_key_cols(),
    clone_map = lookup_clone_map(parent.frame())
) {
  if (length(clone_label) == 0 || anyNA(clone_label)) {
    stop("clone_label must contain at least one non-missing clone label.", call. = FALSE)
  }

  if (length(region_label) != 1 || is.na(region_label)) {
    stop("region_label must be a single non-missing clone label.", call. = FALSE)
  }

  clone_label <- unique(as.character(clone_label))
  region_label <- as.character(region_label)

  mp <- data.table::as.data.table(mutation_presence)
  require_columns(mp, c(mutation_key_cols, "cell_label", "n_cells"), "mutation_presence")

  if ("chr" %in% mutation_key_cols) {
    mp[, chr := as.character(chr)]
  }

  # loci present in EVERY sampled cell of `labels` == clonal(labels)
  clonal_loci <- function(labels) {
    labels <- unique(as.character(labels))
    n_total <- sum_clone_cells(n_cells_per_clone, labels, clone_map)

    if (is.na(n_total) || n_total <= 0) {
      return(unique(mp[0, mutation_key_cols, with = FALSE]))
    }

    presence <- data.table::as.data.table(
      summarise_mutation_presence(
        mutation_presence = tibble::as_tibble(mp),
        required_labels = labels,
        mutation_key_cols = mutation_key_cols
      )
    )

    unique(presence[n_req == n_total, mutation_key_cols, with = FALSE])
  }

  clone_clonal <- clonal_loci(clone_label)

  # clonal(clone) \ clonal(parent): drop only loci FIXED in the parent,
  # regardless of any sibling clones
  if (!is.null(parent_label) && length(parent_label) == 1 && !is.na(parent_label)) {
    parent_clonal <- clonal_loci(parent_label)

    if (nrow(clone_clonal) > 0 && nrow(parent_clonal) > 0) {
      clone_clonal <- data.table::fsetdiff(clone_clonal, parent_clonal)
    }
  }

  # restrict to the branch clone's own diploid (1:1) regions
  regions <- get_diploid_regions(common_normal_regions, region_label)

  if (nrow(regions) == 0) {
    stop("No normal regions available for: ", region_label, call. = FALSE)
  }

  clone_clonal <- unique(
    filter_mutations_in_regions(clone_clonal, regions)[, mutation_key_cols, with = FALSE]
  )

  tibble::as_tibble(clone_clonal)
}

count_alpha_beta <- function(
    mut_clean,
    cna_regions,
    cna_name,
    target_label,
    n_cells_per_clone,
    required_labels = target_label,
    excluded_labels = character(0),
    compatible_labels = character(0),
    redistribute_non_amp = TRUE,
    mutation_presence = NULL,
    mutation_key_cols = default_mutation_key_cols(),
    clone_map = lookup_clone_map(parent.frame())
) {
  if (length(required_labels) == 0) {
    stop("required_labels must contain at least one clone label.", call. = FALSE)
  }

  required_n <- sum_clone_cells(n_cells_per_clone, required_labels, clone_map)

  if (is.na(required_n) || required_n <= 0) {
    stop("Could not determine required cell count for: ", paste(required_labels, collapse = ", "), call. = FALSE)
  }

  region <- tibble::as_tibble(cna_regions)
  require_columns(region, c("cna_name", "target_label", "chr", "begin", "end"), "cna_regions")
  region <- region |>
    dplyr::filter(.data$cna_name == !!cna_name, .data$target_label == !!target_label)

  if (nrow(region) != 1) {
    stop(
      "Expected exactly one CNA region for ",
      cna_name,
      " in ",
      target_label,
      "; found ",
      nrow(region),
      call. = FALSE
    )
  }

  candidates <- candidate_mutations_in_regions(
    mut_clean = mut_clean,
    regions = region,
    labels = required_labels,
    mutation_key_cols = mutation_key_cols
  )

  if (is.null(mutation_presence)) {
    mutation_presence <- build_mutation_presence(mut_clean, mutation_key_cols)
  }

  presence <- summarise_mutation_presence(
    mutation_presence = mutation_presence,
    required_labels = required_labels,
    compatible_labels = compatible_labels,
    excluded_labels = excluded_labels,
    mutation_key_cols = mutation_key_cols
  )

  clonal <- merge(
    data.table::as.data.table(candidates),
    data.table::as.data.table(presence),
    by = mutation_key_cols
  )

  clonal <- clonal[clonal[["n_req"]] == required_n & clonal[["n_exc"]] == 0]

  if (nrow(clonal) == 0) {
    classified <- clonal[, c(mutation_key_cols, "n_req", "n_compatible", "n_exc"), with = FALSE]
    classified[, class := character()]
  } else {
    mut_dt <- data.table::as.data.table(mut_clean)

    if ("chr" %in% mutation_key_cols) {
      mut_dt[, chr := as.character(chr)]
    }

    allele_dt <- mut_dt[mut_dt[["cell_label"]] %in% required_labels]
    allele_dt <- unique(allele_dt[, c(mutation_key_cols, "allele"), with = FALSE])
    allele_dt <- allele_dt[clonal[, mutation_key_cols, with = FALSE], on = mutation_key_cols, nomatch = 0]

    allele_summary <- allele_dt[, {
      alleles <- sort(unique(as.character(allele)))
      .(allele_n = length(alleles), allele_str = paste(alleles, collapse = ","))
    }, by = mutation_key_cols]

    classified <- merge(clonal, allele_summary, by = mutation_key_cols)
    classified[, class := data.table::fifelse(
      allele_n > 1,
      "alpha",
      data.table::fifelse(
        allele_str %in% c("0", "2"),
        "beta",
        data.table::fifelse(allele_str == "1", "non_amp", "other")
      )
    )]
  }

  n_alpha <- sum(classified[["class"]] == "alpha")
  n_beta <- sum(classified[["class"]] == "beta" & classified[["n_compatible"]] == 0)
  n_non <- sum(classified[["class"]] == "non_amp" & classified[["n_compatible"]] == 0)

  den <- n_alpha + n_beta

  if (isTRUE(redistribute_non_amp) && den > 0) {
    n_alpha_noamp <- n_non * n_alpha / den
    n_beta_noamp <- n_non * n_beta / den
  } else {
    n_alpha_noamp <- 0
    n_beta_noamp <- 0
  }

  tibble::tibble(
    cna_name = cna_name,
    target_label = target_label,
    N_alpha_amplified = n_alpha,
    N_beta_amplified = n_beta,
    N_non_amp = n_non,
    N_alpha_noamp = n_alpha_noamp,
    N_beta_noamp = n_beta_noamp,
    N_alpha_total = n_alpha + n_alpha_noamp,
    N_beta_total = n_beta + n_beta_noamp,
    n_clonal_loci = nrow(clonal)
  )
}

candidate_mutations_in_regions <- function(
    mut_clean,
    regions,
    labels,
    mutation_key_cols = default_mutation_key_cols()
) {
  mut_dt <- data.table::as.data.table(mut_clean)
  regions <- tibble::as_tibble(regions)

  require_columns(mut_dt, c("cell_label", mutation_key_cols), "mut_clean")
  require_columns(regions, c("chr", "begin", "end"), "regions")

  if (!"from" %in% mutation_key_cols) {
    stop("mutation_key_cols must include 'from' for interval filtering.", call. = FALSE)
  }

  if ("chr" %in% mutation_key_cols) {
    mut_dt[, chr := as.character(chr)]
  }

  candidate_list <- vector("list", nrow(regions))
  mut_dt <- mut_dt[mut_dt[["cell_label"]] %in% labels]

  for (i in seq_len(nrow(regions))) {
    candidate_list[[i]] <- mut_dt[
      mut_dt[["chr"]] == as.character(regions$chr[[i]]) &
        mut_dt[["from"]] >= regions$begin[[i]] &
        mut_dt[["from"]] <= regions$end[[i]],
      mutation_key_cols,
      with = FALSE
    ]
  }

  candidates <- data.table::rbindlist(candidate_list, use.names = TRUE, fill = TRUE)
  candidates <- unique(candidates)

  if (nrow(candidates) == 0) {
    empty <- mut_dt[0, mutation_key_cols, with = FALSE]
    return(tibble::as_tibble(empty))
  }

  tibble::as_tibble(candidates)
}

filter_mutations_in_regions <- function(dt, regions) {
  dt <- data.table::as.data.table(dt)
  regions <- tibble::as_tibble(regions)

  if (nrow(dt) == 0 || nrow(regions) == 0) {
    return(dt[0])
  }

  require_columns(dt, c("chr", "from"), "mutations")
  require_columns(regions, c("chr", "begin", "end"), "regions")

  keep <- rep(FALSE, nrow(dt))
  chr <- as.character(dt[["chr"]])
  from <- dt[["from"]]

  for (i in seq_len(nrow(regions))) {
    keep <- keep | (
      chr == as.character(regions$chr[[i]]) &
        from >= regions$begin[[i]] &
        from <= regions$end[[i]]
    )
  }

  dt[keep]
}

sum_clone_cells <- function(
    n_cells_per_clone,
    labels,
    clone_map = NULL
) {
  if (length(labels) == 0) {
    return(0L)
  }

  counts <- normalise_clone_counts(n_cells_per_clone, clone_map)
  matched <- counts[counts[["cell_label"]] %in% labels, , drop = FALSE]

  if (nrow(matched) == 0) {
    return(NA_integer_)
  }

  sum(matched[["n_cells_total"]])
}

normalise_clone_counts <- function(n_cells_per_clone, clone_map = NULL) {
  counts <- tibble::as_tibble(n_cells_per_clone)
  require_columns(counts, "n_cells_total", "n_cells_per_clone")

  if (!"cell_label" %in% names(counts)) {
    if (!is.null(clone_map) && all(c("mutant", "cell_label") %in% names(clone_map)) && "mutant" %in% names(counts)) {
      counts <- dplyr::left_join(
        counts,
        tibble::as_tibble(clone_map)[, c("mutant", "cell_label")],
        by = "mutant"
      )
    } else if ("mutant" %in% names(counts)) {
      mutant_labels <- as.character(counts[["mutant"]])
      counts[["cell_label"]] <- ifelse(
        startsWith(mutant_labels, "Clone_"),
        mutant_labels,
        paste0("Clone_", mutant_labels)
      )
    }
  }

  require_columns(counts, "cell_label", "n_cells_per_clone")
  counts[["n_cells_total"]] <- as.integer(counts[["n_cells_total"]])

  counts
}

lookup_clone_map <- function(envir) {
  if (exists("clones", envir = envir, inherits = TRUE)) {
    clone_map <- get("clones", envir = envir, inherits = TRUE)

    if (all(c("mutant", "cell_label") %in% names(clone_map))) {
      return(clone_map)
    }
  }

  NULL
}

lookup_object <- function(name, envir) {
  if (exists(name, envir = envir, inherits = TRUE)) {
    return(get(name, envir = envir, inherits = TRUE))
  }

  NULL
}

require_columns <- function(x, cols, label) {
  missing_cols <- setdiff(cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      "Missing ",
      label,
      " columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

make_posterior_ci_plot <- function(
    posterior_df,
    plot_title,
    x_label,
    x_limits = NULL
) {
  posterior_ci <- posterior_df %>%
    group_by(index, label) %>%
    summarise(
      q025 = quantile(value, 0.025),
      q50  = quantile(value, 0.50),
      q975 = quantile(value, 0.975),
      true_value = first(true_value),
      .groups = "drop"
    )
  
  p <- ggplot(posterior_df, aes(x = value)) +
    geom_rect(
      data = posterior_ci,
      aes(xmin = q025, xmax = q975, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#FFF2A1",
      alpha = 0.45
    ) +
    geom_histogram(
      bins = 40,
      fill = "lightblue",
      color = "black",
      linewidth = 0.25
    ) +
    geom_vline(
      data = posterior_ci,
      aes(xintercept = true_value),
      inherit.aes = FALSE,
      color = "#D7301F",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      data = posterior_ci,
      aes(xintercept = q50),
      inherit.aes = FALSE,
      color = "#1F2937",
      linewidth = 0.9
    ) + 
    facet_wrap(~ label, scales = "free_y", ncol = 2) +
    labs(
      title = plot_title,
      subtitle = "Shaded regions represent 95% credible intervals; dashed red lines indicate true values; solid black lines denote posterior medians",
      x = x_label,
      y = "Posterior frequency"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "grey35", hjust = 0.5),
      plot.margin = margin(8, 30, 8, 30),
      strip.background = element_rect(fill = "grey90", color = "grey35"),
      strip.text = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88"),
      axis.title = element_text(size = 14),
      axis.text = element_text(color = "grey25")
    )
  
  if (!is.null(x_limits)) {
    p <- p + coord_cartesian(xlim = x_limits)
  }
  
  p
}
